//! `tarotx label-by-gene`: label stranded long-read RNA alignments by the annotated gene they
//! most likely come from: how clear the choice is (GeneAssignment_Status) and how the gene is
//! supported (GeneAssignment_SupportLevel). Formerly the separate tool rnalabel.
//!
//! Assumes stranded reads (see `common::INPUT_READS`): only genes on the alignment's strand can be
//! assigned, and opposite-strand genes are reported as antisense.
//!
//! - `annot`: the GFF3 annotation (genes, merged exons, annotated introns);
//! - `blocks`: an alignment's introns and aligned blocks, walked from its cs tag or CIGAR;
//! - `label`: the classification of one alignment (the rules in docs/label-by-gene.md);
//! - `run`: the PAF and TSV runners, and the summary;
//! - `parquet_in`: the Parquet runner;
//! - `table`: the output table (TSV or Parquet).

mod annot;
mod blocks;
mod label;
mod parquet_in;
mod run;
mod table;

use std::time::Instant;

use anyhow::{bail, Result};
use clap::ValueEnum;

use annot::{Annotation, Options};
use blocks::Walk;
use label::{GeneLabel, Params};
use run::{Ctx, TableCols};
use table::Sink;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum InputFormat {
    /// From the file name: .paf, .tsv/.txt (each may end in .gz), .parquet
    Auto,
    Paf,
    Tsv,
    Parquet,
}

#[derive(clap::Args, Debug)]
#[command(
    about = "Assign each alignment to the annotated gene it most likely comes from",
    long_about = "Label each stranded long-read RNA alignment by the annotated gene it most likely comes from: \
                  how clear that choice is (GeneAssignment_Status) and how the gene is supported \
                  (GeneAssignment_SupportLevel).",
    after_help = crate::common::INPUT_READS
)]
pub struct Args {
    /// Alignments: a minimap2 PAF (with cs or cg tags), or a TSV (with a header line) or Parquet table
    /// with a cs or CIGAR column.
    /// Text may be gzipped; "-" reads stdin with --format
    #[arg(short = 'a', long, help_heading = "Inputs")]
    alignments: String,

    /// GENCODE-style GFF3 (may be gzipped). Only gene and exon lines are used, so a slim GFF works
    #[arg(short = 'g', long, help_heading = "Inputs")]
    gff: String,

    /// Input format
    #[arg(long, value_enum, default_value_t = InputFormat::Auto, help_heading = "Inputs")]
    format: InputFormat,

    /// The per-alignment label table (the same columns for every input): .parquet, .tsv, .tsv.gz,
    /// or "-" for TSV on stdout
    #[arg(short = 'o', long, help_heading = "Outputs")]
    output: String,

    /// Also write the status × support-level counts (the summary printed at the end) as a TSV
    #[arg(long, help_heading = "Outputs")]
    summary: Option<String>,

    /// Share of the read's aligned bp on a gene's exons that makes its support exonic (without a
    /// known junction)
    #[arg(long, default_value_t = 0.2, help_heading = "Thresholds")]
    min_exonic_frac: f64,

    /// How far ahead in exonic fraction the top gene must be, at equal matched introns, for
    /// best_of_multiple rather than tied
    #[arg(long, default_value_t = 0.1, help_heading = "Thresholds")]
    tie_margin: f64,

    /// Max bp from an aligned block to the nearest same-strand gene for nearby_gene
    #[arg(long, default_value_t = 5000, help_heading = "Thresholds")]
    near_max_dist: u64,

    /// Leave out genes tagged readthrough_gene
    #[arg(long, help_heading = "Genes")]
    skip_readthrough: bool,

    /// Use only these gene types (comma-separated; default: all)
    #[arg(long, value_delimiter = ',', help_heading = "Genes")]
    gene_types: Vec<String>,

    /// Read ID column, written as query_name
    #[arg(long, default_value = "query_name", help_heading = "Table input columns (TSV, Parquet)")]
    id_col: String,

    /// tp column (P, S, I or i), written as tp_tag; empty when the table has no such column
    #[arg(long, default_value = "tp_tag_A", help_heading = "Table input columns (TSV, Parquet)")]
    tp_col: String,

    /// Reference name column
    #[arg(long, default_value = "target_name_A", help_heading = "Table input columns (TSV, Parquet)")]
    chrom_col: String,

    /// Strand column (+/-; anything else is unmapped)
    #[arg(long, default_value = "strand_A", help_heading = "Table input columns (TSV, Parquet)")]
    strand_col: String,

    /// 0-based start column
    #[arg(long, default_value = "target_start_A", help_heading = "Table input columns (TSV, Parquet)")]
    start_col: String,

    /// End (exclusive) column
    #[arg(long, default_value = "target_end_A", help_heading = "Table input columns (TSV, Parquet)")]
    end_col: String,

    /// minimap2 cs-tag column; the introns are walked from it
    #[arg(long, default_value = "cs_A", help_heading = "Table input columns (TSV, Parquet)")]
    cs_col: String,

    /// Walk the introns from this CIGAR column instead of --cs-col
    #[arg(long, conflicts_with = "cs_col", help_heading = "Table input columns (TSV, Parquet)")]
    cigar_col: Option<String>,
}

fn detect_format(path: &str) -> Result<InputFormat> {
    let p = path.to_ascii_lowercase();
    let p = p.strip_suffix(".gz").unwrap_or(&p);
    Ok(match p.rsplit('.').next() {
        Some("sam" | "bam" | "cram") => {
            bail!("SAM/BAM/CRAM is not read; convert it to PAF first with `maligno sam2paf -U` (keeps unmapped reads)")
        }
        Some("paf") => InputFormat::Paf,
        Some("tsv" | "txt" | "tab") => InputFormat::Tsv,
        Some("parquet") => InputFormat::Parquet,
        _ => bail!("cannot tell the format of '{path}' from its name; set --format paf, tsv or parquet"),
    })
}

pub fn run(args: Args) -> Result<()> {
    let format = match args.format {
        InputFormat::Auto => detect_format(&args.alignments)?,
        f => f,
    };
    if args.summary.as_deref().is_some_and(|s| s == args.output) {
        bail!("--summary and -o are the same file");
    }
    if format == InputFormat::Parquet && args.alignments == "-" {
        bail!("Parquet input can't be read from stdin; give a file");
    }
    let gene_types: Vec<String> = args.gene_types.iter().filter(|d| !d.is_empty()).cloned().collect();

    let t0 = Instant::now();
    let ann = Annotation::load(&args.gff, &Options { gene_types: &gene_types, skip_readthrough: args.skip_readthrough })?;
    eprintln!("loaded the annotation in {:.1} s", t0.elapsed().as_secs_f64());

    let params = Params { min_exonic_frac: args.min_exonic_frac, tie_margin: args.tie_margin, near_max_dist: args.near_max_dist };
    let ctx = Ctx { ann: &ann, params: &params };
    let out = Sink::open(&args.output, GeneLabel::fields())?;
    let (introns, walk) = match &args.cigar_col {
        Some(c) => (c.clone(), Walk::Cigar),
        None => (args.cs_col.clone(), Walk::Cs),
    };
    let cols = TableCols {
        id: args.id_col.clone(),
        tp: args.tp_col.clone(),
        chrom: args.chrom_col.clone(),
        strand: args.strand_col.clone(),
        start: args.start_col.clone(),
        end: args.end_col.clone(),
        introns,
        walk,
    };
    let summary = match format {
        InputFormat::Tsv => run::run_tsv(&args.alignments, &cols, &ctx, out)?,
        InputFormat::Parquet => parquet_in::run_parquet(&args.alignments, &cols, &ctx, out)?,
        InputFormat::Paf | InputFormat::Auto => run::run_paf(&args.alignments, &ctx, out)?,
    };
    summary.print(&ann);
    if let Some(p) = &args.summary {
        summary.write_tsv(p)?;
    }
    eprintln!("done in {:.1} s", t0.elapsed().as_secs_f64());
    Ok(())
}
