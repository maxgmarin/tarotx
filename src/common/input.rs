//! The alignment input, the same for every subcommand: a minimap2 PAF, or a TSV or Parquet table.
//!
//! Every input alignment becomes one `Aln`, in input order: every PAF line and every table row. No
//! alignment is selected or dropped. The 12 leading columns of every output table
//! (`LEAD_COLUMNS`) describe which alignment a row is.

use std::fs::File;

use anyhow::{bail, Context, Result};
use arrow_array::cast::AsArray;
use arrow_array::types::*;
use arrow_array::{Array, RecordBatch};
use arrow_schema::DataType;
use clap::ValueEnum;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ProjectionMask;

use super::io::for_each_line;
use super::strand::Strand;
use super::table::{Kind, Val};

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum InputFormat {
    /// From the file name: .paf, .tsv/.txt (each may end in .gz), .parquet
    Auto,
    Paf,
    Tsv,
    Parquet,
}

impl InputFormat {
    /// The format of `path`: as given, or from the file name for `Auto`.
    pub fn resolve(self, path: &str) -> Result<InputFormat> {
        if self != InputFormat::Auto {
            if self == InputFormat::Parquet && path == "-" {
                bail!("Parquet input can't be read from stdin; give a file");
            }
            return Ok(self);
        }
        let p = path.to_ascii_lowercase();
        let p = p.strip_suffix(".gz").unwrap_or(&p);
        Ok(match p.rsplit('.').next() {
            Some("sam" | "bam" | "cram") => {
                bail!("SAM/BAM/CRAM is not read; convert it to PAF first with `maligno sam2paf -U` (keeps unmapped reads)")
            }
            Some("paf") => InputFormat::Paf,
            Some("tsv" | "txt" | "tab") => InputFormat::Tsv,
            Some("parquet") => InputFormat::Parquet,
            _ => bail!("cannot tell the format of '{path}' from its name; set --in-format paf, tsv or parquet"),
        })
    }
}

/// Which alignment string introns are walked from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Walk {
    Cs,
    Cigar,
}

/// Where a mapped alignment lies on the reference.
pub struct Placed<'a> {
    pub chrom: &'a str,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

/// One input alignment. Values the input doesn't have are None.
pub struct Aln<'a> {
    pub query_name: Option<&'a str>,
    pub tp: Option<&'a str>,
    /// None for an unmapped alignment.
    pub at: Option<Placed<'a>>,
    pub query_len: Option<u64>,
    /// The aligned part of the read, on the read as sequenced (as in PAF), and the mapq. None when
    /// unmapped.
    pub query_start: Option<u64>,
    pub query_end: Option<u64>,
    pub mapq: Option<u64>,
    /// The cs tag or CIGAR to walk introns from: from PAF's cs:Z or cg:Z tag, or a table's intron
    /// column when one was asked for.
    pub introns: Option<(&'a str, Walk)>,
}

/// The leading columns of every output table: which alignment a row is.
pub const LEAD_COLUMNS: [(&str, Kind); 12] = [
    ("query_name", Kind::Str),
    ("target_name", Kind::Str),
    ("target_start", Kind::U64),
    ("target_end", Kind::U64),
    ("strand", Kind::Str),
    ("query_len", Kind::U64),
    ("query_start", Kind::U64),
    ("query_end", Kind::U64),
    ("mapq", Kind::U64),
    ("tp_tag", Kind::Str),
    ("aln_type", Kind::Str),
    ("read_3p_clip", Kind::U64),
];

/// LEAD_COLUMNS, then a subcommand's label columns: all the columns of its output table.
pub fn table_columns(labels: Vec<(String, Kind)>) -> Vec<(String, Kind)> {
    LEAD_COLUMNS.iter().map(|(n, k)| (n.to_string(), *k)).chain(labels).collect()
}

/// The aln_type of an alignment with this tp tag.
fn aln_type_from_tp(tp: Option<&str>) -> &'static str {
    match tp {
        Some("P") => "primary",
        Some("S") => "secondary",
        Some("I" | "i") => "inversion",
        _ => "",
    }
}

impl<'a> Aln<'a> {
    /// The LEAD_COLUMNS cells.
    pub fn lead_cells(&self) -> Vec<Val<'a>> {
        let s = |x: Option<&'a str>| x.filter(|v| !v.is_empty()).map_or(Val::Null, Val::Str);
        let n = |x: Option<u64>| x.map_or(Val::Null, Val::U64);
        let at = self.at.as_ref();
        let aln_type = if at.is_some() { aln_type_from_tp(self.tp) } else { "unmapped" };
        // bases of the read past the alignment's 3' end
        let clip = self.query_len.zip(self.query_end).map(|(len, end)| len.saturating_sub(end));
        vec![
            s(self.query_name),
            s(at.map(|a| a.chrom)),
            n(at.map(|a| a.start)),
            n(at.map(|a| a.end)),
            s(at.map(|a| if a.strand == Strand::Plus { "+" } else { "-" })),
            n(self.query_len),
            n(self.query_start),
            n(self.query_end),
            n(self.mapq),
            s(self.tp),
            s(Some(aln_type)),
            n(clip),
        ]
    }
}

/// Call `f` with every alignment of `path`, in input order. Errors name the line (PAF, TSV) or the
/// row (Parquet).
pub fn for_each_alignment(path: &str, format: InputFormat, cols: &TableCols, mut f: impl FnMut(&Aln) -> Result<()>) -> Result<()> {
    match format {
        InputFormat::Paf | InputFormat::Auto => for_each_paf(path, &mut f),
        InputFormat::Tsv => for_each_tsv(path, cols, &mut f),
        InputFormat::Parquet => for_each_parquet(path, cols, &mut f),
    }
}

// ---------------------------------------------------------------- PAF

/// minimap2 PAF: every line is an alignment. minimap2 writes supplementary alignments as tp:A:P,
/// like primary ones. A line whose strand isn't '+'/'-' (e.g. '*') is unmapped.
fn for_each_paf(path: &str, f: &mut impl FnMut(&Aln) -> Result<()>) -> Result<()> {
    for_each_line(path, |ln, line| {
        if line.is_empty() {
            return Ok(());
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 12 {
            bail!("{path}:{ln}: a PAF line has at least 12 columns, found {}", c.len());
        }
        let num = |i: usize| -> Result<u64> {
            c[i].parse().with_context(|| format!("{path}:{ln}: column {} '{}' is not an integer", i + 1, c[i]))
        };
        let tag = |k: &str| c[12..].iter().find_map(|t| t.strip_prefix(k));
        let introns = match (tag("cs:Z:"), tag("cg:Z:")) {
            (Some(cs), _) => Some((cs, Walk::Cs)),
            (None, Some(cg)) => Some((cg, Walk::Cigar)),
            (None, None) => None,
        };
        let mut a = Aln {
            query_name: Some(c[0]),
            tp: tag("tp:A:"),
            at: None,
            query_len: Some(num(1)?),
            query_start: None,
            query_end: None,
            mapq: None,
            introns,
        };
        if let Some(strand) = Strand::parse(c[4]) {
            a.at = Some(Placed { chrom: c[5], start: num(7)?, end: num(8)?, strand });
            (a.query_start, a.query_end, a.mapq) = (Some(num(2)?), Some(num(3)?), Some(num(11)?));
        }
        f(&a).with_context(|| format!("{path}:{ln}"))
    })
}

// ---------------------------------------------------------------- TSV and Parquet tables

/// The column names of a TSV or Parquet table.
pub struct TableCols {
    pub id: String,
    pub chrom: String,
    pub strand: String,
    pub start: String,
    pub end: String,
    /// Optional: tp_tag is empty without it.
    pub tp: String,
    /// Optional: the read length, query span and mapq. The record columns they feed are empty
    /// without them.
    pub query_len: String,
    pub query_start: String,
    pub query_end: String,
    pub mapq: String,
    /// Required when set: the cs-tag or CIGAR column to walk introns from.
    pub introns: Option<(String, Walk)>,
}

/// Positions of a table's columns in its header.
struct TableIdx {
    id: usize,
    chrom: usize,
    strand: usize,
    start: usize,
    end: usize,
    tp: Option<usize>,
    /// query_len, query_start, query_end, mapq
    record: [Option<usize>; 4],
    introns: Option<(usize, Walk)>,
}

impl TableCols {
    /// Find the columns in `header`: the ID, chrom, strand, start and end (and the intron column,
    /// when set) must be there; the others are used when present.
    fn find(&self, header: &[&str], path: &str) -> Result<TableIdx> {
        let find = |name: &str| header.iter().position(|h| *h == name);
        let mut want = vec![self.id.as_str(), &self.chrom, &self.strand, &self.start, &self.end];
        if let Some((c, _)) = &self.introns {
            want.push(c);
        }
        let missing: Vec<&str> = want.iter().copied().filter(|c| find(c).is_none()).collect();
        if !missing.is_empty() {
            let intron_opts = if self.introns.is_some() { ", and --cs-col or --cigar-col" } else { "" };
            bail!(
                "{path}: no column {missing:?}; set the names with --id-col, --chrom-col, --strand-col, --start-col, \
                 --end-col{intron_opts}"
            );
        }
        let record = [&self.query_len, &self.query_start, &self.query_end, &self.mapq].map(|c| find(c));
        let absent: Vec<&str> = [&self.query_len, &self.query_start, &self.query_end, &self.mapq]
            .iter()
            .zip(record)
            .filter(|(_, i)| i.is_none())
            .map(|(c, _)| c.as_str())
            .collect();
        if !absent.is_empty() && absent.len() < 4 {
            eprintln!(
                "note: {path} has no column {absent:?}, so those leading columns are empty; set the names with \
                 --query-len-col, --query-start-col, --query-end-col, --mapq-col"
            );
        }
        Ok(TableIdx {
            id: find(&self.id).unwrap(),
            chrom: find(&self.chrom).unwrap(),
            strand: find(&self.strand).unwrap(),
            start: find(&self.start).unwrap(),
            end: find(&self.end).unwrap(),
            tp: find(&self.tp),
            record,
            introns: self.introns.as_ref().map(|(c, w)| (find(c).unwrap(), *w)),
        })
    }
}

/// One table row's cells, before the unmapped check: integers are None when empty or not an
/// integer, and `raw` keeps the start and end as written, when they were read as text.
struct Row<'a> {
    id: Option<&'a str>,
    tp: Option<&'a str>,
    chrom: Option<&'a str>,
    strand: Option<&'a str>,
    start: Option<i64>,
    end: Option<i64>,
    raw: [Option<&'a str>; 2],
    record: [Option<i64>; 4],
    introns: Option<(&'a str, Walk)>,
}

impl<'a> Row<'a> {
    /// The alignment. A row is unmapped when its strand isn't '+'/'-', its chrom is empty, or its
    /// start and end are both empty; otherwise its start and end must be a valid 0-based interval.
    fn aln(self) -> Result<Aln<'a>> {
        let empty = |v: Option<i64>, raw: Option<&str>| v.is_none() && raw.is_none_or(str::is_empty);
        let shown = |v: Option<i64>, raw: Option<&str>| raw.map(String::from).or(v.map(|v| v.to_string())).unwrap_or_default();
        let strand = self.strand.and_then(Strand::parse);
        let chrom = self.chrom.filter(|c| !c.is_empty());
        let at = match (strand, chrom, self.start, self.end) {
            (Some(strand), Some(chrom), Some(s), Some(e)) if s >= 0 && e > s => {
                Some(Placed { chrom, start: s as u64, end: e as u64, strand })
            }
            (Some(_), Some(_), _, _) if !(empty(self.start, self.raw[0]) && empty(self.end, self.raw[1])) => bail!(
                "start '{}' and end '{}' are not a valid 0-based interval",
                shown(self.start, self.raw[0]),
                shown(self.end, self.raw[1])
            ),
            _ => None,
        };
        let [query_len, query_start, query_end, mapq] = self.record.map(|x| x.and_then(|x| u64::try_from(x).ok()));
        let mapped = at.is_some();
        Ok(Aln {
            query_name: self.id,
            tp: self.tp,
            at,
            query_len,
            query_start: query_start.filter(|_| mapped),
            query_end: query_end.filter(|_| mapped),
            mapq: mapq.filter(|_| mapped),
            introns: self.introns,
        })
    }
}

/// A TSV table with a header line.
fn for_each_tsv(path: &str, cols: &TableCols, f: &mut impl FnMut(&Aln) -> Result<()>) -> Result<()> {
    let mut idx: Option<TableIdx> = None;
    for_each_line(path, |ln, line| {
        let c: Vec<&str> = line.split('\t').collect();
        let Some(ix) = &idx else {
            idx = Some(cols.find(&c, path)?);
            return Ok(());
        };
        if line.is_empty() {
            return Ok(());
        }
        let get = |i: usize| c.get(i).copied();
        let int = |i: usize| get(i).and_then(|s| s.parse::<i64>().ok());
        let row = Row {
            id: get(ix.id),
            tp: ix.tp.and_then(get),
            chrom: get(ix.chrom),
            strand: get(ix.strand),
            start: int(ix.start),
            end: int(ix.end),
            raw: [get(ix.start), get(ix.end)],
            record: ix.record.map(|i| i.and_then(int)),
            introns: ix.introns.and_then(|(i, w)| get(i).map(|t| (t, w))),
        };
        row.aln().and_then(|a| f(&a)).with_context(|| format!("{path}:{ln}"))
    })?;
    if idx.is_none() {
        bail!("{path}: empty table (no header line)");
    }
    Ok(())
}

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

/// An integer column: any integer type, or integers stored as strings.
enum Ints<'a> {
    Ints(Vec<Option<i64>>),
    Text(Strs<'a>),
}

impl<'a> Ints<'a> {
    fn new(a: &'a dyn Array, name: &str) -> Result<Ints<'a>> {
        macro_rules! prim {
            ($t:ty) => {{
                let p = a.as_primitive::<$t>();
                Ints::Ints((0..p.len()).map(|i| p.is_valid(i).then(|| p.value(i) as i64)).collect())
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
                Some(s) => Ints::Text(s),
                None => bail!("column '{name}' has type {t}; it must hold integers"),
            },
        })
    }

    /// Row i: the value (None when null or not an integer), and the text it was read from.
    fn get(&self, i: usize) -> (Option<i64>, Option<&'a str>) {
        match self {
            Ints::Ints(v) => (v[i], None),
            Ints::Text(s) => {
                let t = s.get(i);
                (t.and_then(|t| t.parse().ok()), t)
            }
        }
    }
}

const BATCH_ROWS: usize = 65_536;

/// A Parquet table. Only the columns used are read.
fn for_each_parquet(path: &str, cols: &TableCols, f: &mut impl FnMut(&Aln) -> Result<()>) -> Result<()> {
    let file = File::open(path).with_context(|| format!("cannot open '{path}'"))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .with_context(|| format!("'{path}' is not a valid Parquet file"))?
        .with_batch_size(BATCH_ROWS);
    let schema = builder.schema().clone();
    let names: Vec<&str> = schema.fields().iter().map(|f| f.name().as_str()).collect();
    let ix = cols.find(&names, path)?;
    let text_cols: Vec<usize> = [ix.id, ix.chrom, ix.strand].into_iter().chain(ix.tp).chain(ix.introns.map(|(i, _)| i)).collect();
    for &i in &text_cols {
        if Strs::new(&arrow_array::new_empty_array(schema.field(i).data_type())).is_none() {
            bail!("{path}: column '{}' has type {}; it must hold strings", names[i], schema.field(i).data_type());
        }
    }
    // Read only the columns used. The projected batches keep them in the file's order.
    let mut read: Vec<usize> = text_cols.iter().copied().chain([ix.start, ix.end]).chain(ix.record.into_iter().flatten()).collect();
    read.sort_unstable();
    read.dedup();
    let pos = |i: usize| read.binary_search(&i).unwrap();
    let mask = ProjectionMask::roots(builder.parquet_schema(), read.iter().copied());

    let mut row = 0usize;
    for batch in builder.with_projection(mask).build().with_context(|| format!("reading '{path}'"))? {
        let batch: RecordBatch = batch.with_context(|| format!("reading a row group of '{path}'"))?;
        let col = |i: usize| batch.column(pos(i)).as_ref();
        let strs = |i: usize| Strs::new(col(i)).unwrap();
        let (id, chrom, strand) = (strs(ix.id), strs(ix.chrom), strs(ix.strand));
        let tp = ix.tp.map(strs);
        let introns = ix.introns.map(|(i, w)| (strs(i), w));
        let start = Ints::new(col(ix.start), names[ix.start])?;
        let end = Ints::new(col(ix.end), names[ix.end])?;
        let mut record: [Option<Ints>; 4] = [None, None, None, None];
        for (k, c) in ix.record.iter().enumerate() {
            if let Some(c) = *c {
                record[k] = Some(Ints::new(col(c), names[c])?);
            }
        }
        for i in 0..batch.num_rows() {
            row += 1;
            let ((s, rs), (e, re)) = (start.get(i), end.get(i));
            let r = Row {
                id: id.get(i),
                tp: tp.as_ref().and_then(|t| t.get(i)),
                chrom: chrom.get(i),
                strand: strand.get(i),
                start: s,
                end: e,
                raw: [rs, re],
                record: [0, 1, 2, 3].map(|k| record[k].as_ref().and_then(|c| c.get(i).0)),
                introns: introns.as_ref().and_then(|(t, w)| t.get(i).map(|x| (x, *w))),
            };
            r.aln().and_then(|a| f(&a)).with_context(|| format!("{path}: row {row}"))?;
        }
    }
    Ok(())
}
