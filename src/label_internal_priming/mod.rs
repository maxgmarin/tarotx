//! `tarotx label-internal-priming`: label long-read RNA-seq alignments whose 3' end sits next to
//! a same-strand genomic poly(A) tract, and flag likely internal priming (poly(A)-adjacent and,
//! when an annotation is given, not near an annotated transcript end). Formerly the separate tool
//! iprime.
//!
//! Assumes stranded, trimmed reads (see `common::INPUT_READS`): the strand says which end of the
//! alignment is the 3' end, and trimming makes that end the point where the transcript ended, which
//! is what the poly(A) tracts and annotated ends are measured from.
//!
//! - `points`: the reference points: polyscan poly(A) tracts and annotated transcript ends;
//! - `label`: the labels of one alignment's 3' end, and the summary;
//! - `run`: the PAF and TSV runners, and the tagged PAF output;
//! - `parquet_in`: the Parquet runner;
//! - `table`: the output table (TSV or Parquet).
//!
//! Every input alignment gets one row, in input order (every PAF line, every table row); nothing
//! is selected among a read's alignments. SAM/BAM is not read: convert it to PAF first.

mod label;
mod parquet_in;
mod points;
mod run;
mod table;

use std::time::Instant;

use anyhow::{bail, Result};
use clap::ValueEnum;

use label::{Params, Refs};
use points::{Tracts, TxEnds};
use run::TableCols;
use table::OutSpec;

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
    about = "Flag alignments whose 3' end sits next to a genomic poly(A) tract (likely internal priming)",
    long_about = "Label each alignment's 3' end by its distance to the nearest same-strand genomic poly(A) tract \
                  edge and, optionally, to the nearest annotated transcript end; flag likely internal priming.",
    after_help = crate::common::INPUT_READS
)]
pub struct Args {
    /// Alignments: a minimap2 PAF, or a TSV (with a header line) or Parquet table. Every alignment
    /// gets one row, in input order. Text may be gzipped; "-" reads stdin with --format
    #[arg(short = 'a', long, help_heading = "Inputs")]
    alignments: String,

    /// polyscan poly(A) tracts: BED6, or 5-column bedtools-merge output (may be gzipped)
    #[arg(short = 'p', long, help_heading = "Inputs")]
    polya_tracts: String,

    /// GENCODE-style GFF3 annotation (may be gzipped). Without it, the flag is poly(A) adjacency alone
    #[arg(short = 'g', long, help_heading = "Inputs")]
    gff: Option<String>,

    /// Input format
    #[arg(long, value_enum, default_value_t = InputFormat::Auto, help_heading = "Inputs")]
    format: InputFormat,

    /// Per-alignment table with every label and its supporting values: .parquet, .tsv, .tsv.gz,
    /// or "-" for TSV on stdout
    #[arg(short = 'o', long, help_heading = "Outputs")]
    output: String,

    /// PAF input only: every input line, with pd/td/ip tags appended to the mapped ones
    /// (.paf, or .paf.gz for gzip)
    #[arg(long, help_heading = "Outputs")]
    tagged_out: Option<String>,

    /// Write only rows where internal_priming_flag is False to the table
    #[arg(long, help_heading = "Outputs")]
    exclude_flagged: bool,

    /// Also write the counts (the summary printed at the end) as a TSV: one row per combination of
    /// the label columns, and the unmapped alignments
    #[arg(long, help_heading = "Outputs")]
    summary: Option<String>,

    /// Max bp from the 3' end to either edge of a poly(A) tract
    #[arg(long, default_value_t = 3, help_heading = "Thresholds")]
    polya_max_dist: i64,

    /// Max bp from the 3' end to an annotated transcript end
    #[arg(long, default_value_t = 50, help_heading = "Thresholds")]
    tes_max_dist: i64,

    /// Use only transcripts tagged 'basic' as annotated ends
    #[arg(long, help_heading = "Annotated transcript ends")]
    gencode_basic_only: bool,

    /// Keep transcripts tagged mRNA_end_NF (incomplete 3' end)
    #[arg(long, help_heading = "Annotated transcript ends")]
    keep_end_nf: bool,

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

    /// Read length column. With all four of --query-len-col, --query-start-col, --query-end-col and
    /// --mapq-col, a table gets the same record columns as PAF, read_3p_clip included
    #[arg(long, default_value = "query_len", help_heading = "Table input columns (TSV, Parquet)")]
    query_len_col: String,

    /// Query start column: 0-based, on the read as sequenced (as in PAF)
    #[arg(long, default_value = "query_start_A", help_heading = "Table input columns (TSV, Parquet)")]
    query_start_col: String,

    /// Query end (exclusive) column, on the read as sequenced
    #[arg(long, default_value = "query_end_A", help_heading = "Table input columns (TSV, Parquet)")]
    query_end_col: String,

    /// Mapping quality column
    #[arg(long, default_value = "mapq_A", help_heading = "Table input columns (TSV, Parquet)")]
    mapq_col: String,
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

/// --tagged-out writes the input PAF with tags, so it needs PAF input, a PAF file name, and a path
/// of its own.
fn check_tagged_out(tagged: Option<&str>, format: InputFormat, table: &str) -> Result<()> {
    let Some(p) = tagged else { return Ok(()) };
    if format != InputFormat::Paf {
        bail!("--tagged-out needs PAF input");
    }
    if ![".paf", ".paf.gz"].iter().any(|e| p.to_ascii_lowercase().ends_with(e)) {
        bail!("--tagged-out '{p}': the input is PAF, so this output is PAF too; name it .paf or .paf.gz");
    }
    if p == table {
        bail!("--tagged-out and -o are the same file '{p}'");
    }
    Ok(())
}

pub fn run(args: Args) -> Result<()> {
    let format = match args.format {
        InputFormat::Auto => detect_format(&args.alignments)?,
        f => f,
    };
    let tagged = args.tagged_out.as_deref();
    check_tagged_out(tagged, format, &args.output)?;
    if let Some(p) = args.summary.as_deref() {
        if p == args.output || tagged == Some(p) {
            bail!("--summary and another output are the same file '{p}'");
        }
    }
    if format == InputFormat::Parquet && args.alignments == "-" {
        bail!("Parquet input can't be read from stdin; give a file");
    }
    let t0 = Instant::now();
    // Both files are usually large and gzipped; decode them side by side.
    let (tracts, tx_ends) = std::thread::scope(|s| {
        let gff = args.gff.as_ref().map(|g| s.spawn(|| TxEnds::load(g, args.gencode_basic_only, !args.keep_end_nf)));
        let tracts = Tracts::load(&args.polya_tracts);
        let tx_ends = gff.map(|h| h.join().expect("annotation loader panicked")).transpose();
        (tracts, tx_ends)
    });
    let (tracts, tx_ends) = (tracts?, tx_ends?);
    let refs = Refs {
        tracts,
        tx_ends,
        params: Params { polya_max_dist: args.polya_max_dist, tes_max_dist: args.tes_max_dist },
    };
    eprintln!("loaded references in {:.1} s", t0.elapsed().as_secs_f64());

    let out = OutSpec { path: &args.output, exclude_flagged: args.exclude_flagged };
    let cols = TableCols {
        id: args.id_col.clone(),
        tp: args.tp_col.clone(),
        chrom: args.chrom_col.clone(),
        strand: args.strand_col.clone(),
        start: args.start_col.clone(),
        end: args.end_col.clone(),
        query_len: args.query_len_col.clone(),
        query_start: args.query_start_col.clone(),
        query_end: args.query_end_col.clone(),
        mapq: args.mapq_col.clone(),
    };
    let summary = match format {
        InputFormat::Tsv => run::run_tsv(&args.alignments, &cols, &refs, &out)?,
        InputFormat::Parquet => parquet_in::run_parquet(&args.alignments, &cols, &refs, &out)?,
        InputFormat::Paf | InputFormat::Auto => run::run_paf(&args.alignments, tagged, &refs, &out)?,
    };
    summary.print(&refs);
    if let Some(p) = &args.summary {
        summary.write_tsv(p, &refs)?;
    }
    eprintln!("done in {:.1} s", t0.elapsed().as_secs_f64());
    Ok(())
}
