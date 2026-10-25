//! `tarotx label-by-gene`: label stranded long-read RNA alignments by the annotated gene they
//! most likely come from: how clear the choice is (GeneAssignment_Status) and how the gene is
//! supported (GeneAssignment_SupportLevel). Formerly the separate tool rnalabel.
//!
//! Assumes stranded reads (see the README's "Input reads" section): only genes on the alignment's
//! strand can be assigned, and opposite-strand genes are reported as antisense.
//!
//! - `annot`: the GFF3 annotation (genes, merged exons, annotated introns);
//! - `blocks`: an alignment's introns and aligned blocks, walked from its cs tag or CIGAR;
//! - `label`: the classification of one alignment (the rules in docs/label-by-gene.md);
//! - `run`: labelling every input alignment, and the summary.
//!
//! The input, the leading columns and the output table are shared with every subcommand
//! (`crate::common`).

mod annot;
mod blocks;
mod label;
mod run;

use std::time::Instant;

use anyhow::Result;

use crate::common::args::{check_outputs, InputArgs, TableColumnArgs};
use crate::common::input::Walk;
use annot::{Annotation, Options};
use label::Params;

#[derive(clap::Args, Debug)]
#[command(
    about = "Assign each alignment to its most likely annotated gene",
    long_about = None
)]
pub struct Args {
    #[command(flatten)]
    input: InputArgs,

    /// GENCODE GFF3 (may be gzipped). Its gene and exon lines are used
    #[arg(short = 'g', long, help_heading = "Inputs")]
    gff: String,

    /// Table of inferred gene labels per alignment (.parquet, .tsv, .tsv.gz, or "-" for TSV on stdout)
    #[arg(short = 'o', long, help_heading = "Outputs")]
    output: String,

    /// Also write the counts of each status and support level as a TSV
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

    #[command(flatten)]
    table: TableColumnArgs,

    /// minimap2 cs-tag column the introns are walked from (PAF input reads its cs:Z or cg:Z tag)
    #[arg(long, default_value = "cs_A", help_heading = "Table input columns (TSV, Parquet)")]
    cs_col: String,

    /// Walk the introns from this CIGAR column instead of --cs-col
    #[arg(long, conflicts_with = "cs_col", help_heading = "Table input columns (TSV, Parquet)")]
    cigar_col: Option<String>,
}

pub fn run(args: Args) -> Result<()> {
    let format = args.input.format()?;
    check_outputs(&args.output, args.summary.as_deref())?;
    let gene_types: Vec<String> = args.gene_types.iter().filter(|d| !d.is_empty()).cloned().collect();

    let t0 = Instant::now();
    let ann = Annotation::load(&args.gff, &Options { gene_types: &gene_types, skip_readthrough: args.skip_readthrough })?;
    eprintln!("loaded the annotation in {:.1} s", t0.elapsed().as_secs_f64());

    let params = Params { min_exonic_frac: args.min_exonic_frac, tie_margin: args.tie_margin, near_max_dist: args.near_max_dist };
    let introns = match &args.cigar_col {
        Some(c) => (c.clone(), Walk::Cigar),
        None => (args.cs_col.clone(), Walk::Cs),
    };
    let cols = args.table.cols(Some(introns));
    let summary = run::run(&args.input.alignments, format, &cols, &ann, &params, &args.output)?;
    summary.print(&ann);
    if let Some(p) = &args.summary {
        summary.write_tsv(p)?;
    }
    eprintln!("done in {:.1} s", t0.elapsed().as_secs_f64());
    Ok(())
}
