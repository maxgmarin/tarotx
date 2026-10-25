//! Code shared by every subcommand.

pub mod io;
pub mod strand;

/// The input reads every subcommand assumes, shown at the end of each --help. The README's "Input
/// reads" section explains what relies on each assumption.
pub const INPUT_READS: &str = "Input reads: tarotx is designed for stranded long-read RNA-seq reads, oriented 5'->3' \
along their transcript, whose adapters and poly(A) tails have been trimmed (e.g. PacBio Iso-Seq/Kinnex reads \
after lima and isoseq refine). The alignment strand is taken as the transcript strand, and the alignment's end \
as the transcript's end. These assumptions are not checked: unstranded or untrimmed reads are labelled anyway, \
and their labels can be wrong. Such reads are not supported yet.";

use parquet::file::metadata::KeyValue;

/// Key-value metadata for every Parquet table tarotx writes: its version, and the subcommand
/// that wrote the table (e.g. "label-by-gene").
pub fn parquet_metadata(command: &str) -> Vec<KeyValue> {
    vec![
        KeyValue::new("tarotx_version".to_string(), env!("CARGO_PKG_VERSION").to_string()),
        KeyValue::new("tarotx_command".to_string(), command.to_string()),
    ]
}
