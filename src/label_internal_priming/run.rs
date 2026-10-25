//! The PAF and TSV runners. Each streams its alignments, labels them, writes the table rows (and,
//! for PAF, the tagged alignments) and returns the summary counts.
//!
//! One row per input alignment, in input order: every PAF line and every table row, with no
//! selection among a read's alignments. Every input gets the same six leading columns, as
//! label-by-gene writes them (which alignment a row is). PAF adds six columns from the alignment
//! record, and so do tables that have the read length, query span and mapq (e.g. a maligno compare
//! table). Unmapped alignments get a row with empty coordinates and empty labels.

use anyhow::{bail, Context, Result};
use std::io::Write;

use crate::common::io::{for_each_line, open_output};
use crate::common::strand::Strand;
use super::label::{Label, Refs, Summary};
use super::table::{Kind, OutSpec, Sink, Val};

/// Leading columns for every input: which alignment a row is.
const LEAD_COLUMNS: [(&str, Kind); 6] = [
    ("query_name", Kind::Str),
    ("tp_tag", Kind::Str),
    ("target_name", Kind::Str),
    ("target_start", Kind::U64),
    ("target_end", Kind::U64),
    ("strand", Kind::Str),
];

/// More leading columns from the alignment record: for PAF input, and for tables that have the read
/// length, query span and mapq.
const RECORD_COLUMNS: [(&str, Kind); 6] = [
    ("query_len", Kind::U64),
    ("query_start", Kind::U64),
    ("query_end", Kind::U64),
    ("mapq", Kind::U64),
    ("aln_type", Kind::Str),
    ("read_3p_clip", Kind::U64),
];

/// The leading columns: LEAD_COLUMNS, then RECORD_COLUMNS when the input has them.
pub fn lead_columns(record_columns: bool) -> Vec<(String, Kind)> {
    let extra: &[(&str, Kind)] = if record_columns { &RECORD_COLUMNS } else { &[] };
    LEAD_COLUMNS.iter().chain(extra).map(|(n, k)| (n.to_string(), *k)).collect()
}

/// Where a mapped alignment lies.
struct Placed<'a> {
    chrom: &'a str,
    start: u64,
    end: u64,
    strand: Strand,
}

impl Placed<'_> {
    fn label<'r>(&self, refs: &'r Refs) -> Label<'r> {
        let strand = if self.strand == Strand::Plus { "+" } else { "-" };
        refs.label(self.chrom, strand, self.start as i64, self.end as i64)
    }
}

/// The LEAD_COLUMNS cells. The coordinates are empty for an unmapped alignment.
fn lead<'a>(query_name: Option<&'a str>, tp: Option<&'a str>, at: Option<&Placed<'a>>) -> Vec<Val<'a>> {
    let str_or_null = |s: Option<&'a str>| s.filter(|s| !s.is_empty()).map_or(Val::Null, Val::Str);
    let mut v = vec![str_or_null(query_name), str_or_null(tp)];
    match at {
        Some(p) => v.extend([
            Val::Str(p.chrom),
            Val::U64(p.start),
            Val::U64(p.end),
            Val::Str(if p.strand == Strand::Plus { "+" } else { "-" }),
        ]),
        None => v.extend([Val::Null; 4]),
    }
    v
}

/// One PAF line. Query coordinates are on the read as sequenced, so for an oriented read
/// `query_len - query_end` bases lie past its 3' end.
struct Aln<'a> {
    query_name: &'a str,
    tp: Option<&'a str>,
    query_len: u64,
    /// None for an unmapped record.
    mapped: Option<Mapped<'a>>,
}

struct Mapped<'a> {
    at: Placed<'a>,
    query_start: u64,
    query_end: u64,
    mapq: u32,
    aln_type: &'a str,
}

/// The aln_type of a PAF line or table row, from its tp tag.
fn aln_type_from_tp(tp: Option<&str>) -> &'static str {
    match tp {
        Some("P") => "primary",
        Some("S") => "secondary",
        Some("I" | "i") => "inversion",
        _ => "",
    }
}

/// The RECORD_COLUMNS cells. `span` is a mapped alignment's query start, query end and mapq (None
/// when unmapped, which leaves them empty and makes aln_type "unmapped"). read_3p_clip needs the
/// read length and the query end.
fn record_cells<'a>(query_len: Option<u64>, span: Option<[Option<u64>; 3]>, aln_type: &'a str) -> [Val<'a>; 6] {
    let opt = |x: Option<u64>| x.map_or(Val::Null, Val::U64);
    match span {
        Some([query_start, query_end, mapq]) => [
            opt(query_len),
            opt(query_start),
            opt(query_end),
            opt(mapq),
            Val::Str(aln_type),
            opt(query_len.zip(query_end).map(|(len, end)| len.saturating_sub(end))),
        ],
        None => [opt(query_len), Val::Null, Val::Null, Val::Null, Val::Str("unmapped"), Val::Null],
    }
}

/// Label one alignment record (None when unmapped) and write its row.
fn emit<'r>(a: &Aln, refs: &'r Refs, out: &mut Sink, summary: &mut Summary) -> Result<Option<Label<'r>>> {
    let l = a.mapped.as_ref().map(|m| m.at.label(refs));
    summary.add(refs, a.mapped.as_ref().map(|m| m.at.chrom).zip(l.as_ref()));
    let mut cells = lead(Some(a.query_name), a.tp, a.mapped.as_ref().map(|m| &m.at));
    let span = a.mapped.as_ref().map(|m| [Some(m.query_start), Some(m.query_end), Some(m.mapq as u64)]);
    cells.extend(record_cells(Some(a.query_len), span, a.mapped.as_ref().map_or("", |m| m.aln_type)));
    out.row(&cells, l.as_ref(), refs)?;
    Ok(l)
}

// ---------------------------------------------------------------- TSV and Parquet tables

/// Columns of a TSV or Parquet table.
pub struct TableCols {
    pub id: String,
    /// Optional: when the table has no such column, tp_tag is empty.
    pub tp: String,
    pub chrom: String,
    pub strand: String,
    pub start: String,
    pub end: String,
    /// Optional, all four together: the read length, the query span (on the read as sequenced,
    /// as in PAF) and the mapq. With them, a table gets RECORD_COLUMNS too.
    pub query_len: String,
    pub query_start: String,
    pub query_end: String,
    pub mapq: String,
}

/// Positions of a table's columns in its header.
pub struct TableIdx {
    pub id: usize,
    pub tp: Option<usize>,
    pub chrom: usize,
    pub strand: usize,
    pub start: usize,
    pub end: usize,
    /// The query_len, query_start, query_end and mapq columns, when the table has all four.
    pub record: Option<[usize; 4]>,
}

impl TableCols {
    /// Find the ID, chrom, strand, start and end columns in `header`, and the tp and record columns
    /// if present.
    pub fn find(&self, header: &[&str], path: &str) -> Result<TableIdx> {
        let want = [self.id.as_str(), &self.chrom, &self.strand, &self.start, &self.end];
        let find = |name: &str| header.iter().position(|h| *h == name);
        let missing: Vec<&str> = want.iter().copied().filter(|c| find(c).is_none()).collect();
        if !missing.is_empty() {
            bail!("{path}: no column {missing:?}; set the names with --id-col, --chrom-col, --strand-col, --start-col, --end-col");
        }
        let [id, chrom, strand, start, end] = want.map(|c| find(c).unwrap());
        let rec = [&self.query_len, &self.query_start, &self.query_end, &self.mapq].map(|c| (c.as_str(), find(c)));
        let record = match rec.iter().filter(|(_, i)| i.is_none()).map(|(c, _)| *c).collect::<Vec<_>>() {
            missing if missing.is_empty() => Some(rec.map(|(_, i)| i.unwrap())),
            missing if missing.len() < 4 => {
                eprintln!(
                    "note: {path} has no column {missing:?}, so the record columns (query_len ... read_3p_clip) are \
                     left out; set the names with --query-len-col, --query-start-col, --query-end-col, --mapq-col"
                );
                None
            }
            _ => None,
        };
        Ok(TableIdx { id, tp: find(&self.tp), chrom, strand, start, end, record })
    }
}

/// One table row's cells: ID, tp, chrom, strand, and the start and end as integers (None when
/// empty or not an integer) and, when they were read as text, as written.
pub struct TableRow<'a> {
    pub id: Option<&'a str>,
    pub tp: Option<&'a str>,
    pub chrom: Option<&'a str>,
    pub strand: Option<&'a str>,
    pub start: Option<i64>,
    pub end: Option<i64>,
    pub raw: [Option<&'a str>; 2],
    /// query_len, query_start, query_end and mapq, when the table has those columns.
    pub record: Option<[Option<i64>; 4]>,
}

/// Label one table row and write it. A row is unmapped when its strand isn't '+'/'-', its chrom is
/// empty, or its start and end are both empty; otherwise its start and end must be a valid 0-based
/// interval.
pub fn table_row(r: &TableRow, refs: &Refs, out: &mut Sink, summary: &mut Summary) -> Result<()> {
    let strand = r.strand.and_then(Strand::parse);
    let chrom = r.chrom.filter(|c| !c.is_empty());
    let empty = |v: Option<i64>, raw: Option<&str>| v.is_none() && raw.is_none_or(str::is_empty);
    let shown = |v: Option<i64>, raw: Option<&str>| raw.map(String::from).or(v.map(|v| v.to_string())).unwrap_or_default();
    let at = match (strand, chrom, r.start, r.end) {
        (Some(strand), Some(chrom), Some(s), Some(e)) if s >= 0 && e > s => {
            Some(Placed { chrom, start: s as u64, end: e as u64, strand })
        }
        (Some(_), Some(_), _, _) if !(empty(r.start, r.raw[0]) && empty(r.end, r.raw[1])) => bail!(
            "start '{}' and end '{}' are not a valid 0-based interval",
            shown(r.start, r.raw[0]),
            shown(r.end, r.raw[1])
        ),
        _ => None,
    };
    let l = at.as_ref().map(|p| p.label(refs));
    summary.add(refs, at.as_ref().map(|p| p.chrom).zip(l.as_ref()));
    let mut cells = lead(r.id, r.tp, at.as_ref());
    if let Some(rec) = r.record {
        let [query_len, query_start, query_end, mapq] = rec.map(|x| x.and_then(|x| u64::try_from(x).ok()));
        let span = at.as_ref().map(|_| [query_start, query_end, mapq]);
        cells.extend(record_cells(query_len, span, aln_type_from_tp(r.tp)));
    }
    out.row(&cells, l.as_ref(), refs)
}

/// A TSV table with a header line. The output starts once the header shows which leading columns
/// the table can fill.
pub fn run_tsv(path: &str, cols: &TableCols, refs: &Refs, spec: &OutSpec) -> Result<Summary> {
    let mut summary = Summary::default();
    let mut started: Option<(TableIdx, Sink)> = None;
    for_each_line(path, |ln, line| {
        let f: Vec<&str> = line.split('\t').collect();
        let Some((ix, out)) = &mut started else {
            let ix = cols.find(&f, path)?;
            let out = spec.open(&lead_columns(ix.record.is_some()), refs)?;
            started = Some((ix, out));
            return Ok(());
        };
        if line.is_empty() {
            return Ok(());
        }
        let get = |i: usize| f.get(i).copied();
        let int = |i: usize| get(i).and_then(|s| s.parse::<i64>().ok());
        let r = TableRow {
            id: get(ix.id),
            tp: ix.tp.and_then(get),
            chrom: get(ix.chrom),
            strand: get(ix.strand),
            start: int(ix.start),
            end: int(ix.end),
            raw: [get(ix.start), get(ix.end)],
            record: ix.record.map(|r| r.map(int)),
        };
        table_row(&r, refs, out, &mut summary).with_context(|| format!("{path}:{ln}"))
    })?;
    let Some((_, out)) = started else { bail!("{path}: empty table (no header line)") };
    out.finish()?;
    Ok(summary)
}

// ---------------------------------------------------------------- PAF

/// minimap2 PAF. Every line gets a row, as given: minimap2 marks supplementary alignments `tp:A:P`
/// like primary ones, so the pieces of a split read are separate rows. A line whose strand isn't
/// '+'/'-' is unmapped. With `tagged`, every line is also written there with the tags appended.
pub fn run_paf(path: &str, tagged: Option<&str>, refs: &Refs, spec: &OutSpec) -> Result<Summary> {
    let mut out = spec.open(&lead_columns(true), refs)?;
    let mut summary = Summary::default();
    let mut tagged = tagged.map(open_output).transpose()?;
    for_each_line(path, |ln, line| {
        if line.is_empty() {
            return Ok(());
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            bail!("{path}:{ln}: a PAF line has at least 12 columns, found {}", f.len());
        }
        let num = |i: usize| -> Result<u64> {
            f[i].parse().with_context(|| format!("{path}:{ln}: column {} '{}' is not an integer", i + 1, f[i]))
        };
        let tp = f[12..].iter().find_map(|t| t.strip_prefix("tp:A:"));
        let aln_type = aln_type_from_tp(tp);
        let mapped = match Strand::parse(f[4]) {
            None => None,
            Some(strand) => Some(Mapped {
                at: Placed { chrom: f[5], start: num(7)?, end: num(8)?, strand },
                query_start: num(2)?,
                query_end: num(3)?,
                mapq: num(11)? as u32,
                aln_type,
            }),
        };
        let a = Aln { query_name: f[0], tp, query_len: num(1)?, mapped };
        let l = emit(&a, refs, &mut out, &mut summary)?;
        if let Some(w) = &mut tagged {
            // an unmapped line is copied as it is
            match &l {
                Some(l) => w.write_all(tagged_paf_line(&f, l).as_bytes())?,
                None => writeln!(w, "{line}")?,
            }
        }
        Ok(())
    })?;
    if let Some(w) = tagged {
        w.finish()?;
    }
    out.finish()?;
    Ok(summary)
}

/// The PAF line with this subcommand's tags appended, replacing any pd/td/ip tags it already had.
fn tagged_paf_line(f: &[&str], l: &Label) -> String {
    let mut line: Vec<String> = f
        .iter()
        .enumerate()
        .filter(|(i, t)| *i < 12 || !(t.starts_with("pd:") || t.starts_with("td:") || t.starts_with("ip:")))
        .map(|(_, t)| t.to_string())
        .collect();
    if let Some((d, _)) = l.polya {
        line.push(format!("pd:i:{d}"));
    }
    if let Some((d, _)) = l.tx_end {
        line.push(format!("td:i:{d}"));
    }
    line.push(format!("ip:A:{}", if l.flag { 'Y' } else { 'N' }));
    let mut s = line.join("\t");
    s.push('\n');
    s
}
