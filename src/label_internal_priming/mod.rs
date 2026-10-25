//! `tarotx label-internal-priming`: label long-read RNA-seq alignments whose 3' end sits next to
//! a same-strand genomic poly(A) tract, and flag likely internal priming (poly(A)-adjacent and,
//! when an annotation is given, not near an annotated transcript end). Formerly the separate tool
//! iprime.
//!
//! Assumes stranded, trimmed reads (see the README's "Input reads" section): the strand says which
//! end of the alignment is the 3' end, and trimming makes that end the point where the transcript
//! ended, which is what the poly(A) tracts and annotated ends are measured from.
//!
//! - `points`: the reference points: polyscan poly(A) tracts and annotated transcript ends;
//! - `label`: the labels of one alignment's 3' end, and the summary;
//! - `run`: labelling every input alignment (and the tagged-PAF code kept for later use).
//!
//! The input, the leading columns and the output table are shared with every subcommand
//! (`crate::common`).

mod label;
mod points;
mod run;

use std::time::Instant;

use anyhow::Result;

use crate::common::args::{check_outputs, InputArgs, TableColumnArgs};
use label::{Params, Refs};
use points::{Tracts, TxEnds};

#[derive(clap::Args, Debug)]
#[command(
    about = "Flag reads with alignments that are likely due to RT internal priming artifacts",
    long_about = None
)]
pub struct Args {
    #[command(flatten)]
    input: InputArgs,

    /// poly(A) tracts from polyscan --mode tracts: BED6, or BED8 with --extended (may be gzipped)
    #[arg(short = 'p', long, help_heading = "Inputs")]
    polya_tracts: String,

    /// GENCODE GFF3 (may be gzipped). Its transcript lines give the annotated 3' ends. Without it,
    /// the flag is poly(A) adjacency alone
    #[arg(short = 'g', long, help_heading = "Inputs")]
    gff: Option<String>,

    /// Table of inferred internal-priming labels per alignment (.parquet, .tsv, .tsv.gz, or "-" for TSV on
    /// stdout)
    #[arg(short = 'o', long, help_heading = "Outputs")]
    output: String,

    /// Also write the counts of flagged, poly(A)-adjacent and near-annotated-end alignments as a TSV
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

    #[command(flatten)]
    table: TableColumnArgs,
}

pub fn run(args: Args) -> Result<()> {
    let format = args.input.format()?;
    check_outputs(&args.output, args.summary.as_deref())?;
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

    let cols = args.table.cols(None);
    let summary = run::run(&args.input.alignments, format, &cols, &refs, &args.output)?;
    summary.print(&refs);
    if let Some(p) = &args.summary {
        summary.write_tsv(p, &refs)?;
    }
    eprintln!("done in {:.1} s", t0.elapsed().as_secs_f64());
    Ok(())
}
