//! Parquet table input.

use anyhow::{bail, Context, Result};
use arrow_array::cast::AsArray;
use arrow_array::types::*;
use arrow_array::{Array, RecordBatch};
use arrow_schema::DataType;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ProjectionMask;

use super::run::{begin_aln, lead, Ctx, Summary, TableCols};
use super::table::Sink;

/// String cells of one column, for the string layouts pandas, polars, DuckDB and maligno write.
pub enum Strs<'a> {
    Utf8(&'a arrow_array::StringArray),
    Large(&'a arrow_array::LargeStringArray),
    View(&'a arrow_array::StringViewArray),
    /// Dictionary-encoded (e.g. a pandas categorical): key per row, and the values.
    Dict(Vec<usize>, &'a dyn Array, Box<Strs<'a>>),
}

impl<'a> Strs<'a> {
    pub fn new(a: &'a dyn Array) -> Option<Strs<'a>> {
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

    pub fn get(&self, i: usize) -> Option<&'a str> {
        match self {
            Strs::Utf8(a) => a.is_valid(i).then(|| a.value(i)),
            Strs::Large(a) => a.is_valid(i).then(|| a.value(i)),
            Strs::View(a) => a.is_valid(i).then(|| a.value(i)),
            Strs::Dict(keys, a, values) => if a.is_valid(i) { values.get(keys[i]) } else { None },
        }
    }
}

/// Integer coordinates of one column: any integer type, or integers stored as strings. Nulls,
/// empty strings and non-integers (e.g. "*" on an unmapped row) are None.
pub fn ints(a: &dyn Array, name: &str) -> Result<Vec<Option<i64>>> {
    macro_rules! prim {
        ($t:ty) => {{
            let p = a.as_primitive::<$t>();
            (0..p.len()).map(|i| if p.is_valid(i) { Some(p.value(i) as i64) } else { None }).collect()
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
            Some(s) => (0..a.len()).map(|i| s.get(i).and_then(|x| x.parse::<i64>().ok())).collect(),
            None => bail!("column '{name}' has type {t}; coordinates must be integers"),
        },
    })
}

/// A cell as pandas' to_csv writes it.
pub fn cell(a: &dyn Array, i: usize) -> Result<String> {
    if a.is_null(i) {
        return Ok(String::new());
    }
    macro_rules! prim {
        ($t:ty) => {
            a.as_primitive::<$t>().value(i).to_string()
        };
    }
    Ok(match a.data_type() {
        DataType::Boolean => (if a.as_boolean().value(i) { "True" } else { "False" }).to_string(),
        DataType::Int8 => prim!(Int8Type),
        DataType::Int16 => prim!(Int16Type),
        DataType::Int32 => prim!(Int32Type),
        DataType::Int64 => prim!(Int64Type),
        DataType::UInt8 => prim!(UInt8Type),
        DataType::UInt16 => prim!(UInt16Type),
        DataType::UInt32 => prim!(UInt32Type),
        DataType::UInt64 => prim!(UInt64Type),
        DataType::Float32 => super::table::float(a.as_primitive::<Float32Type>().value(i) as f64),
        DataType::Float64 => super::table::float(a.as_primitive::<Float64Type>().value(i)),
        t => match Strs::new(a) {
            Some(s) => s.get(i).unwrap_or("").to_string(),
            None => bail!("a column of type {t} can't be written to TSV; write Parquet instead"),
        },
    })
}

const BATCH_ROWS: usize = 65_536;

/// A Parquet table. Only the ID, tp, coordinate and cs (or CIGAR) columns are read.
pub fn run_parquet(path: &str, cols: &TableCols, ctx: &Ctx, mut out: Sink) -> Result<Summary> {
    let file = std::fs::File::open(path).with_context(|| format!("cannot open '{path}'"))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .with_context(|| format!("'{path}' is not a valid Parquet file"))?
        .with_batch_size(BATCH_ROWS);
    let schema = builder.schema().clone();
    let names: Vec<&str> = schema.fields().iter().map(|f| f.name().as_str()).collect();
    let ix = cols.find(&names, path)?;
    for i in [ix.id, ix.chrom, ix.strand, ix.introns].into_iter().chain(ix.tp) {
        if Strs::new(&arrow_array::new_empty_array(schema.field(i).data_type())).is_none() {
            bail!("{path}: column '{}' has type {}; it must hold strings", names[i], schema.field(i).data_type());
        }
    }
    // A projected batch holds the columns read, in file order.
    let mut read: Vec<usize> = [ix.id, ix.chrom, ix.strand, ix.start, ix.end, ix.introns].into_iter().chain(ix.tp).collect();
    read.sort_unstable();
    read.dedup();
    let at = |i: usize| read.binary_search(&i).unwrap();
    let mask = ProjectionMask::roots(builder.parquet_schema(), read.iter().copied());
    begin_aln(&mut out)?;

    let mut summary = Summary::default();
    let mut row = 0u64; // for error messages
    for batch in builder.with_projection(mask).build().with_context(|| format!("reading '{path}'"))? {
        let batch: RecordBatch = batch.with_context(|| format!("reading a row group of '{path}'"))?;
        let col = |i: usize| batch.column(at(i)).as_ref();
        let strs = |i: usize| Strs::new(col(i)).unwrap();
        let (id, chrom, strand, introns) = (strs(ix.id), strs(ix.chrom), strs(ix.strand), strs(ix.introns));
        let tp = ix.tp.map(strs);
        let (start, end) = (ints(col(ix.start), names[ix.start])?, ints(col(ix.end), names[ix.end])?);
        for i in 0..batch.num_rows() {
            // the start and end cells as text, for the unmapped check and error messages
            let (rs, re) = (cell(col(ix.start), i).ok(), cell(col(ix.end), i).ok());
            let (l, placed) = cols
                .label_row(ctx, &mut summary, [chrom.get(i), strand.get(i), introns.get(i)], start[i], end[i], [rs.as_deref(), re.as_deref()])
                .with_context(|| format!("{path}: row {}", row + 1))?;
            out.row(&lead(id.get(i), tp.as_ref().and_then(|t| t.get(i)), placed.as_ref()), &l.values(ctx.ann))?;
            row += 1;
        }
    }
    out.finish()?;
    Ok(summary)
}
