//! Parquet table input. Only the columns used are read: ID, tp, chrom, strand, start and end, and
//! the record columns (read length, query span, mapq) when the table has them.

use std::fs::File;

use anyhow::{bail, Context, Result};
use arrow_array::cast::AsArray;
use arrow_array::types::*;
use arrow_array::{Array, RecordBatch};
use arrow_schema::DataType;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ProjectionMask;

use super::label::{Refs, Summary};
use super::run::{table_row, TableCols, TableRow};
use super::run::lead_columns;
use super::table::OutSpec;

const BATCH_ROWS: usize = 65_536;

/// String cells of one column, for the string layouts pandas, polars, DuckDB and maligno write.
enum Strs<'a> {
    Utf8(&'a arrow_array::StringArray),
    Large(&'a arrow_array::LargeStringArray),
    View(&'a arrow_array::StringViewArray),
    /// Dictionary-encoded (e.g. a pandas categorical): key per row, and the values.
    Dict(Vec<usize>, &'a dyn Array, Box<Strs<'a>>),
}

impl<'a> Strs<'a> {
    fn new(a: &'a dyn Array) -> Option<Strs<'a>> {
        Some(match a.data_type() {
            DataType::Utf8 => Strs::Utf8(a.as_string::<i32>()),
            DataType::LargeUtf8 => Strs::Large(a.as_string::<i64>()),
            DataType::Utf8View => Strs::View(a.as_string_view()),
            DataType::Dictionary(_, v) if matches!(**v, DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View) => {
                let d = a.as_any_dictionary();
                Strs::Dict(d.normalized_keys(), a, Box::new(Strs::new(d.values().as_ref())?))
            }
            _ => return None,
        })
    }

    fn get(&self, i: usize) -> Option<&'a str> {
        match self {
            Strs::Utf8(a) => a.is_valid(i).then(|| a.value(i)),
            Strs::Large(a) => a.is_valid(i).then(|| a.value(i)),
            Strs::View(a) => a.is_valid(i).then(|| a.value(i)),
            Strs::Dict(keys, a, values) => if a.is_valid(i) { values.get(keys[i]) } else { None },
        }
    }
}

/// A coordinate column: any integer type, or integers stored as strings.
enum Coords<'a> {
    Ints(Vec<Option<i64>>),
    Text(Strs<'a>),
}

impl<'a> Coords<'a> {
    fn new(a: &'a dyn Array, name: &str) -> Result<Coords<'a>> {
        macro_rules! prim {
            ($t:ty) => {{
                let p = a.as_primitive::<$t>();
                Coords::Ints((0..p.len()).map(|i| p.is_valid(i).then(|| p.value(i) as i64)).collect())
            }};
        }
        Ok(match a.data_type() {
            DataType::Int8 => prim!(Int8Type),
            DataType::Int16 => prim!(Int16Type),
            DataType::Int32 => prim!(Int32Type),
            DataType::Int64 => prim!(Int64Type),
            DataType::UInt8 => prim!(UInt8Type),
            DataType::UInt16 => prim!(UInt16Type),
            DataType::UInt32 => prim!(UInt32Type),
            DataType::UInt64 => prim!(UInt64Type),
            t => match Strs::new(a) {
                Some(s) => Coords::Text(s),
                None => bail!("column '{name}' has type {t}; coordinates must be integers"),
            },
        })
    }

    /// Row i: the value (None when null or not an integer), and the text it was read from.
    fn get(&self, i: usize) -> (Option<i64>, Option<&'a str>) {
        match self {
            Coords::Ints(v) => (v[i], None),
            Coords::Text(s) => {
                let t = s.get(i);
                (t.and_then(|t| t.parse().ok()), t)
            }
        }
    }
}

pub fn run_parquet(path: &str, cols: &TableCols, refs: &Refs, spec: &OutSpec) -> Result<Summary> {
    let file = File::open(path).with_context(|| format!("cannot open '{path}'"))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .with_context(|| format!("'{path}' is not a valid Parquet file"))?
        .with_batch_size(BATCH_ROWS);
    let schema = builder.schema().clone();
    let names: Vec<&str> = schema.fields().iter().map(|f| f.name().as_str()).collect();
    let ix = cols.find(&names, path)?;
    let text_cols: Vec<usize> = [ix.id, ix.chrom, ix.strand].into_iter().chain(ix.tp).collect();
    for &i in &text_cols {
        if Strs::new(&arrow_array::new_empty_array(schema.field(i).data_type())).is_none() {
            bail!("{path}: column '{}' has type {}; it must hold strings", names[i], schema.field(i).data_type());
        }
    }
    // Read only the columns used. The projected batches keep them in the file's order.
    let mut read: Vec<usize> = text_cols.iter().copied().chain([ix.start, ix.end]).chain(ix.record.into_iter().flatten()).collect();
    read.sort_unstable();
    read.dedup();
    let pos = |i: usize| read.iter().position(|&r| r == i).unwrap();
    let mask = ProjectionMask::roots(builder.parquet_schema(), read.iter().copied());

    let mut out = spec.open(&lead_columns(ix.record.is_some()), refs)?;
    let mut summary = Summary::default();
    let mut row = 0usize;
    for batch in builder.with_projection(mask).build().with_context(|| format!("reading '{path}'"))? {
        let batch: RecordBatch = batch.with_context(|| format!("reading a row group of '{path}'"))?;
        let strs = |i: usize| Strs::new(batch.column(pos(i)).as_ref()).unwrap();
        let (id, chrom, strand) = (strs(ix.id), strs(ix.chrom), strs(ix.strand));
        let tp = ix.tp.map(strs);
        let start = Coords::new(batch.column(pos(ix.start)).as_ref(), names[ix.start])?;
        let end = Coords::new(batch.column(pos(ix.end)).as_ref(), names[ix.end])?;
        let record: Option<Vec<Coords>> = ix
            .record
            .map(|r| r.iter().map(|&c| Coords::new(batch.column(pos(c)).as_ref(), names[c])).collect())
            .transpose()?;
        for i in 0..batch.num_rows() {
            let ((s, rs), (e, re)) = (start.get(i), end.get(i));
            let r = TableRow {
                id: id.get(i),
                tp: tp.as_ref().and_then(|t| t.get(i)),
                chrom: chrom.get(i),
                strand: strand.get(i),
                start: s,
                end: e,
                raw: [rs, re],
                record: record.as_ref().map(|r| [0, 1, 2, 3].map(|k| r[k].get(i).0)),
            };
            row += 1;
            table_row(&r, refs, &mut out, &mut summary).with_context(|| format!("{path}: row {row}"))?;
        }
    }
    out.finish()?;
    Ok(summary)
}
