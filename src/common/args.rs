//! Command-line options every subcommand shares, flattened into each subcommand's `Args`, so the
//! input, output and table-column options are the same everywhere.

use anyhow::{bail, Result};

use super::input::{InputFormat, TableCols, Walk};

#[derive(clap::Args, Debug)]
pub struct InputArgs {
    /// Input alignment info: in PAF or a custom table (TSV or PARQUET).
    #[arg(short = 'a', long, help_heading = "Inputs")]
    pub alignments: String,

    /// Input format: auto, paf, tsv, parquet
    #[arg(long = "in-format", value_enum, default_value_t = InputFormat::Auto, hide_possible_values = true, help_heading = "Inputs")]
    pub format: InputFormat,
}

impl InputArgs {
    /// The input format: as given, or from the file name.
    pub fn format(&self) -> Result<InputFormat> {
        self.format.resolve(&self.alignments)
    }
}

/// Each subcommand defines its own `-o/--output` and `--summary`, whose help says what they hold.
/// The summary must not overwrite the table.
pub fn check_outputs(output: &str, summary: Option<&str>) -> Result<()> {
    if summary.is_some_and(|s| s == output) {
        bail!("--summary and -o are the same file");
    }
    Ok(())
}

const TABLE: &str = "Table input columns (TSV, Parquet)";

#[derive(clap::Args, Debug)]
pub struct TableColumnArgs {
    /// Read ID column, written as query_name
    #[arg(long, default_value = "query_name", help_heading = TABLE)]
    pub id_col: String,

    /// Reference name column
    #[arg(long, default_value = "target_name_A", help_heading = TABLE)]
    pub chrom_col: String,

    /// Strand column (+/-; anything else is unmapped)
    #[arg(long, default_value = "strand_A", help_heading = TABLE)]
    pub strand_col: String,

    /// 0-based start column
    #[arg(long, default_value = "target_start_A", help_heading = TABLE)]
    pub start_col: String,

    /// End (exclusive) column
    #[arg(long, default_value = "target_end_A", help_heading = TABLE)]
    pub end_col: String,

    /// Read length column (optional; query_len and read_3p_clip are empty without it)
    #[arg(long, default_value = "query_len", help_heading = TABLE)]
    pub query_len_col: String,

    /// Query start column, 0-based on the read as sequenced, as in PAF (optional)
    #[arg(long, default_value = "query_start_A", help_heading = TABLE)]
    pub query_start_col: String,

    /// Query end (exclusive) column, on the read as sequenced (optional)
    #[arg(long, default_value = "query_end_A", help_heading = TABLE)]
    pub query_end_col: String,

    /// Mapping quality column (optional)
    #[arg(long, default_value = "mapq_A", help_heading = TABLE)]
    pub mapq_col: String,

    /// tp column (P, S, I or i), written as tp_tag (optional)
    #[arg(long, default_value = "tp_tag_A", help_heading = TABLE)]
    pub tp_col: String,
}

impl TableColumnArgs {
    /// The table columns, with `introns` (a cs or CIGAR column) for a subcommand that needs them.
    pub fn cols(&self, introns: Option<(String, Walk)>) -> TableCols {
        TableCols {
            id: self.id_col.clone(),
            chrom: self.chrom_col.clone(),
            strand: self.strand_col.clone(),
            start: self.start_col.clone(),
            end: self.end_col.clone(),
            tp: self.tp_col.clone(),
            query_len: self.query_len_col.clone(),
            query_start: self.query_start_col.clone(),
            query_end: self.query_end_col.clone(),
            mapq: self.mapq_col.clone(),
            introns,
        }
    }
}
