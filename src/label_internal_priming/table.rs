//! The per-alignment output table, written as TSV or Parquet.
//!
//! A table is the leading columns (which alignment) followed by the label columns, one row at a
//! time.

use std::fs::File;
use std::io::Write;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow_array::builder::{BooleanBuilder, Int64Builder, StringBuilder, UInt64Builder};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;

use crate::common::io::{open_output, Output};
use super::label::{Label, Refs};

/// Rows per Parquet row group. Set on the writer as well, because ArrowWriter otherwise keeps
/// appending writes to one row group of up to 1M rows and holds it in memory until it closes.
const ROW_GROUP_ROWS: usize = 65_536;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Str,
    U64,
    I64,
    Bool,
}

impl Kind {
    fn arrow(self) -> DataType {
        match self {
            Kind::Str => DataType::Utf8,
            Kind::U64 => DataType::UInt64,
            Kind::I64 => DataType::Int64,
            Kind::Bool => DataType::Boolean,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Val<'a> {
    Str(&'a str),
    U64(u64),
    I64(i64),
    Bool(bool),
    Null,
}

impl Val<'_> {
    /// The TSV cell: booleans as True/False (as pandas writes them), null as empty.
    fn write_tsv(self, w: &mut dyn Write) -> std::io::Result<()> {
        match self {
            Val::Str(s) => w.write_all(s.as_bytes()),
            Val::U64(x) => write!(w, "{x}"),
            Val::I64(x) => write!(w, "{x}"),
            Val::Bool(b) => w.write_all(if b { b"True" } else { b"False" }),
            Val::Null => Ok(()),
        }
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
    Bool(BooleanBuilder),
}

impl Builder {
    fn new(k: Kind) -> Builder {
        match k {
            Kind::Str => Builder::Str(StringBuilder::new()),
            Kind::U64 => Builder::U64(UInt64Builder::new()),
            Kind::I64 => Builder::I64(Int64Builder::new()),
            Kind::Bool => Builder::Bool(BooleanBuilder::new()),
        }
    }

    fn append(&mut self, v: Val) {
        match (self, v) {
            (Builder::Str(b), Val::Str(s)) => b.append_value(s),
            (Builder::U64(b), Val::U64(x)) => b.append_value(x),
            (Builder::I64(b), Val::I64(x)) => b.append_value(x),
            (Builder::Bool(b), Val::Bool(x)) => b.append_value(x),
            (Builder::Str(b), Val::Null) => b.append_null(),
            (Builder::U64(b), Val::Null) => b.append_null(),
            (Builder::I64(b), Val::Null) => b.append_null(),
            (Builder::Bool(b), Val::Null) => b.append_null(),
            (_, v) => unreachable!("value {v:?} does not match its column type"),
        }
    }

    fn finish(&mut self) -> ArrayRef {
        match self {
            Builder::Str(b) => Arc::new(b.finish()),
            Builder::U64(b) => Arc::new(b.finish()),
            Builder::I64(b) => Arc::new(b.finish()),
            Builder::Bool(b) => Arc::new(b.finish()),
        }
    }
}

struct ParquetOut {
    writer: ArrowWriter<File>,
    schema: SchemaRef,
    lead: Vec<Builder>,
    labels: Vec<Builder>,
    buffered: usize,
}

impl ParquetOut {
    fn create(path: &str, lead: &[(String, Kind)], refs: &Refs) -> Result<ParquetOut> {
        let label_fields = refs.fields();
        let fields: Vec<Field> =
            lead.iter().chain(label_fields.iter()).map(|(n, k)| Field::new(n, k.arrow(), true)).collect();
        let schema = Arc::new(Schema::new(fields));
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(ZstdLevel::try_new(3)?))
            .set_max_row_group_row_count(Some(ROW_GROUP_ROWS))
            .set_key_value_metadata(Some(crate::common::parquet_metadata("label-internal-priming")))
            .build();
        let file = File::create(path).with_context(|| format!("cannot create '{path}'"))?;
        let writer = ArrowWriter::try_new(file, schema.clone(), Some(props)).context("creating the Parquet writer")?;
        let lead = lead.iter().map(|(_, k)| Builder::new(*k)).collect();
        let labels = label_fields.iter().map(|(_, k)| Builder::new(*k)).collect();
        Ok(ParquetOut { writer, schema, lead, labels, buffered: 0 })
    }

    fn row(&mut self, lead: &[Val], label: &[Val]) -> Result<()> {
        for (b, v) in self.lead.iter_mut().zip(lead) {
            b.append(*v);
        }
        for (b, v) in self.labels.iter_mut().zip(label) {
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
        let arrays: Vec<ArrayRef> = self.lead.iter_mut().chain(self.labels.iter_mut()).map(Builder::finish).collect();
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

pub struct Sink {
    out: Out,
    exclude_flagged: bool,
}

/// Where the table goes. A runner opens it once it knows the leading columns, which for a table
/// input depend on its header.
pub struct OutSpec<'a> {
    pub path: &'a str,
    pub exclude_flagged: bool,
}

impl OutSpec<'_> {
    pub fn open(&self, lead: &[(String, Kind)], refs: &Refs) -> Result<Sink> {
        Sink::open(self.path, lead, refs, self.exclude_flagged)
    }
}

impl Sink {
    /// Start the table with these leading columns, then the label columns: Parquet for names
    /// ending in ".parquet"; otherwise TSV (stdout for "-", gzip for ".gz").
    pub fn open(path: &str, lead: &[(String, Kind)], refs: &Refs, exclude_flagged: bool) -> Result<Sink> {
        let out = if is_parquet_path(path) {
            Out::Parquet(Box::new(ParquetOut::create(path, lead, refs)?))
        } else {
            let mut w = open_output(path)?;
            let cols: Vec<String> = lead.iter().chain(refs.fields().iter()).map(|(n, _)| n.clone()).collect();
            writeln!(w, "{}", cols.join("\t"))?;
            Out::Tsv(w)
        };
        Ok(Sink { out, exclude_flagged })
    }

    /// One row: the leading cells, and the label (None for an unmapped alignment, whose label
    /// cells are empty).
    pub fn row(&mut self, lead: &[Val], l: Option<&Label>, refs: &Refs) -> Result<()> {
        if self.exclude_flagged && l.is_some_and(|l| l.flag) {
            return Ok(());
        }
        let labels = refs.values(l);
        match &mut self.out {
            Out::Tsv(w) => {
                for (i, v) in lead.iter().chain(&labels).enumerate() {
                    if i > 0 {
                        w.write_all(b"\t")?;
                    }
                    v.write_tsv(w)?;
                }
                w.write_all(b"\n")?;
            }
            Out::Parquet(p) => p.row(lead, &labels)?,
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
