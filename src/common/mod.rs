//! Code shared by every subcommand.

pub mod io;
pub mod strand;

use parquet::file::metadata::KeyValue;

/// Key-value metadata for every Parquet table tarotx writes: its version, and the subcommand
/// that wrote the table (e.g. "label-by-gene").
pub fn parquet_metadata(command: &str) -> Vec<KeyValue> {
    vec![
        KeyValue::new("tarotx_version".to_string(), env!("CARGO_PKG_VERSION").to_string()),
        KeyValue::new("tarotx_command".to_string(), command.to_string()),
    ]
}
