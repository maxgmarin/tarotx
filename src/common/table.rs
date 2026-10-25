//! The per-alignment output table, written as TSV or Parquet: the leading columns (which alignment,
//! `input::LEAD_COLUMNS`) followed by a subcommand's label columns, one row at a time.

use std::fs::File;
use std::io::Write;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow_array::builder::{BooleanBuilder, Float64Builder, Int64Builder, StringBuilder, UInt64Builder};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;

use super::io::{open_output, Output};

/// Rows per Parquet row group. Set on the writer as well, because ArrowWriter otherwise keeps
/// appending writes to one row group of up to 1M rows and holds it in memory until it closes.
const ROW_GROUP_ROWS: usize = 65_536;

/// A column's type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Str,
    U64,
    I64,
    F64,
    Bool,
}

impl Kind {
    fn arrow(self) -> DataType {
        match self {
            Kind::Str => DataType::Utf8,
            Kind::U64 => DataType::UInt64,
            Kind::I64 => DataType::Int64,
            Kind::F64 => DataType::Float64,
            Kind::Bool => DataType::Boolean,
        }
    }
}

/// One cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Val<'a> {
    Str(&'a str),
    U64(u64),
    I64(i64),
    F64(f64),
    Bool(bool),
    Null,
}

impl Val<'_> {
    /// The TSV cell, as pandas writes it: floats in Python's repr form, booleans as True/False, null
    /// as empty.
    fn write_tsv(self, w: &mut dyn Write) -> std::io::Result<()> {
        match self {
            Val::Str(s) => w.write_all(s.as_bytes()),
            Val::U64(x) => write!(w, "{x}"),
            Val::I64(x) => write!(w, "{x}"),
            Val::F64(x) => w.write_all(float(x).as_bytes()),
            Val::Bool(b) => w.write_all(if b { b"True" } else { b"False" }),
            Val::Null => Ok(()),
        }
    }
}

/// Floats as pandas writes them (Python's repr): NaN empty, shortest round-trip digits, fixed
/// notation for exponents -4 to 15 (with ".0" when integral), else e.g. "4.4e-07".
pub fn float(x: f64) -> String {
    if x.is_nan() {
        return String::new();
    }
    if x.is_infinite() {
        return (if x > 0.0 { "inf" } else { "-inf" }).to_string();
    }
    let sci = format!("{:e}", x.abs()); // shortest digits, e.g. "4.4478841421113913e-7"
    let (mant, exp) = sci.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    let digits: String = mant.chars().filter(|c| *c != '.').collect();
    let sign = if x.is_sign_negative() { "-" } else { "" };
    if (-4..16).contains(&exp) {
        let body = if exp >= 0 {
            let int_len = exp as usize + 1;
            let int = format!("{:0<int_len$}", &digits[..digits.len().min(int_len)]);
            let frac = digits.get(int_len..).unwrap_or("");
            format!("{int}.{}", if frac.is_empty() { "0" } else { frac })
        } else {
            format!("0.{}{digits}", "0".repeat((-exp - 1) as usize))
        };
        format!("{sign}{body}")
    } else {
        let (d0, rest) = digits.split_at(1);
        let rest = if rest.is_empty() { String::new() } else { format!(".{rest}") };
        format!("{sign}{d0}{rest}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
    }
}

pub fn is_parquet_path(path: &str) -> bool {
    path.to_ascii_lowercase().ends_with(".parquet")
}

// ---------------------------------------------------------------- Parquet writer

enum Builder {
    Str(StringBuilder),
    U64(UInt64Builder),
    I64(Int64Builder),
    F64(Float64Builder),
    Bool(BooleanBuilder),
}

impl Builder {
    fn new(k: Kind) -> Builder {
        match k {
            Kind::Str => Builder::Str(StringBuilder::new()),
            Kind::U64 => Builder::U64(UInt64Builder::new()),
            Kind::I64 => Builder::I64(Int64Builder::new()),
            Kind::F64 => Builder::F64(Float64Builder::new()),
            Kind::Bool => Builder::Bool(BooleanBuilder::new()),
        }
    }

    fn append(&mut self, v: Val) {
        match (self, v) {
            (Builder::Str(b), Val::Str(s)) => b.append_value(s),
            (Builder::U64(b), Val::U64(x)) => b.append_value(x),
            (Builder::I64(b), Val::I64(x)) => b.append_value(x),
            (Builder::F64(b), Val::F64(x)) => b.append_option((!x.is_nan()).then_some(x)),
            (Builder::Bool(b), Val::Bool(x)) => b.append_value(x),
            (Builder::Str(b), Val::Null) => b.append_null(),
            (Builder::U64(b), Val::Null) => b.append_null(),
            (Builder::I64(b), Val::Null) => b.append_null(),
            (Builder::F64(b), Val::Null) => b.append_null(),
            (Builder::Bool(b), Val::Null) => b.append_null(),
            (_, v) => unreachable!("value {v:?} does not match its column type"),
        }
    }

    fn finish(&mut self) -> ArrayRef {
        match self {
            Builder::Str(b) => Arc::new(b.finish()),
            Builder::U64(b) => Arc::new(b.finish()),
            Builder::I64(b) => Arc::new(b.finish()),
            Builder::F64(b) => Arc::new(b.finish()),
            Builder::Bool(b) => Arc::new(b.finish()),
        }
    }
}

struct ParquetOut {
    writer: ArrowWriter<File>,
    schema: SchemaRef,
    columns: Vec<Builder>,
    buffered: usize,
}

impl ParquetOut {
    fn create(path: &str, columns: &[(String, Kind)], command: &str) -> Result<ParquetOut> {
        let fields: Vec<Field> = columns.iter().map(|(n, k)| Field::new(n, k.arrow(), true)).collect();
        let schema = Arc::new(Schema::new(fields));
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(ZstdLevel::try_new(3)?))
            .set_max_row_group_row_count(Some(ROW_GROUP_ROWS))
            .set_key_value_metadata(Some(super::parquet_metadata(command)))
            .build();
        let file = File::create(path).with_context(|| format!("cannot create '{path}'"))?;
        let writer = ArrowWriter::try_new(file, schema.clone(), Some(props)).context("creating the Parquet writer")?;
        let columns = columns.iter().map(|(_, k)| Builder::new(*k)).collect();
        Ok(ParquetOut { writer, schema, columns, buffered: 0 })
    }

    fn row(&mut self, cells: &[Val]) -> Result<()> {
        for (b, v) in self.columns.iter_mut().zip(cells) {
            b.append(*v);
        }
        self.buffered += 1;
        if self.buffered >= ROW_GROUP_ROWS {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        if self.buffered == 0 {
            return Ok(());
        }
        let arrays: Vec<ArrayRef> = self.columns.iter_mut().map(Builder::finish).collect();
        let batch = RecordBatch::try_new(self.schema.clone(), arrays).context("assembling a Parquet row group")?;
        self.writer.write(&batch).context("writing a Parquet row group")?;
        self.buffered = 0;
        Ok(())
    }

    fn close(mut self) -> Result<()> {
        self.flush()?;
        self.writer.close().context("closing the Parquet file")?;
        Ok(())
    }
}

// ---------------------------------------------------------------- the sink

enum Out {
    Tsv(Output),
    Parquet(Box<ParquetOut>),
}

/// The output table.
pub struct Sink {
    out: Out,
    n_columns: usize,
}

impl Sink {
    /// Start a table with these columns (leading, then labels). Parquet for names ending in
    /// ".parquet", recording `command` (the subcommand) in its metadata; otherwise TSV (stdout for
    /// "-", gzip for ".gz").
    pub fn open(path: &str, columns: &[(String, Kind)], command: &str) -> Result<Sink> {
        let out = if is_parquet_path(path) {
            Out::Parquet(Box::new(ParquetOut::create(path, columns, command)?))
        } else {
            let mut w = open_output(path)?;
            let names: Vec<&str> = columns.iter().map(|(n, _)| n.as_str()).collect();
            writeln!(w, "{}", names.join("\t"))?;
            Out::Tsv(w)
        };
        Ok(Sink { out, n_columns: columns.len() })
    }

    /// One row, with one cell per column.
    pub fn row(&mut self, cells: &[Val]) -> Result<()> {
        debug_assert_eq!(cells.len(), self.n_columns);
        match &mut self.out {
            Out::Tsv(w) => {
                for (i, v) in cells.iter().enumerate() {
                    if i > 0 {
                        w.write_all(b"\t")?;
                    }
                    v.write_tsv(w)?;
                }
                w.write_all(b"\n")?;
            }
            Out::Parquet(p) => p.row(cells)?,
        }
        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        match self.out {
            Out::Tsv(w) => w.finish()?,
            Out::Parquet(p) => p.close()?,
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::float;

    #[test]
    fn floats_like_pandas() {
        for (x, want) in [
            (1.0, "1.0"),
            (0.25, "0.25"),
            (-2.5, "-2.5"),
            (0.0, "0.0"),
            (100.0, "100.0"),
            (123456.789, "123456.789"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (4.4478841421113913e-07, "4.4478841421113913e-07"),
            (1e15, "1000000000000000.0"),
            (1e16, "1e+16"),
            (1.5e300, "1.5e+300"),
        ] {
            assert_eq!(float(x), want, "{x:e}");
        }
        assert_eq!(float(f64::NAN), "");
    }
}
