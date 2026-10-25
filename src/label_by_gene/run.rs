//! One runner per input type. Each streams its alignments, labels them, writes the table rows and
//! returns the summary counts.

use anyhow::{bail, Context, Result};

use super::annot::Annotation;
use crate::common::strand::Strand;
use super::blocks::{Structure, Walk};
use std::io::Write;

use crate::common::io::{for_each_line, open_output};
use super::label::{Support, GeneLabel, Params, Status};
use super::table::{Kind, Sink, Val};

/// Everything a runner labels against.
pub struct Ctx<'a> {
    pub ann: &'a Annotation,
    pub params: &'a Params,
}

/// Counts for the summary printed at the end of a run.
#[derive(Default)]
pub struct Summary {
    n: u64,
    /// [status][support level], with index 3 for no support level, indexed like `Status::ALL`, `Support::ALL`.
    counts: [[u64; 4]; 7],
    /// Mapped alignments on chromosomes the annotation doesn't have.
    off_annotation: u64,
    mapped: u64,
}

impl Summary {
    fn add(&mut self, ann: &Annotation, chrom: Option<&str>, l: &GeneLabel) {
        self.n += 1;
        let s = Status::ALL.iter().position(|&s| s == l.status).unwrap();
        let e = l.support.map_or(3, |e| Support::ALL.iter().position(|&k| k == e).unwrap());
        self.counts[s][e] += 1;
        if let Some(c) = chrom {
            self.mapped += 1;
            self.off_annotation += !ann.has_chrom(c) as u64;
        }
    }

    pub fn print(&self, ann: &Annotation) {
        let pct = |x: u64| if self.n == 0 { 0.0 } else { x as f64 * 100.0 / self.n as f64 };
        eprintln!("{} alignments labelled, {} genes, {} annotated introns", self.n, ann.genes.len(), ann.n_introns);
        if self.mapped > 0 && self.off_annotation == self.mapped {
            eprintln!("WARNING: no alignment is on a chromosome named in the annotation (chr1 vs 1?)");
        } else if self.off_annotation as f64 > self.mapped as f64 * 0.01 {
            eprintln!(
                "note: {:.1}% of mapped alignments are on chromosomes absent from the annotation",
                self.off_annotation as f64 * 100.0 / self.mapped as f64
            );
        }
        // status by support level, as count (percent of all alignments)
        let cell = |x: u64| format!("{x} ({:.1}%)", pct(x));
        let header: Vec<&str> = Support::ALL.iter().map(|e| e.name()).chain(["no gene", "total"]).collect();
        let w = Status::ALL.iter().map(|s| s.name().len()).max().unwrap_or(0) + 2;
        eprintln!("{:<w$}{}", "", header.iter().map(|h| format!("{h:>18}")).collect::<String>());
        for (s, row) in Status::ALL.iter().zip(self.counts) {
            let cells: String = row.iter().chain([&row.iter().sum()]).map(|&x| format!("{:>18}", cell(x))).collect();
            eprintln!("{:<w$}{cells}", s.name());
        }
    }
}

impl Summary {
    /// The status × support-level counts as a TSV: every combination that can occur, zeros included.
    pub fn write_tsv(&self, path: &str) -> Result<()> {
        let mut w = open_output(path)?;
        writeln!(w, "GeneAssignment_Status\tGeneAssignment_SupportLevel\tcount\tpercent")?;
        for (si, s) in Status::ALL.iter().enumerate() {
            let levels: &[usize] = match s {
                Status::UniqueGene | Status::Tied => &[0, 1, 2],
                Status::BestOfMultiple => &[0, 1], // intronic-only genes always tie
                _ => &[3],
            };
            for &e in levels {
                let n = self.counts[si][e];
                let pct = if self.n == 0 { 0.0 } else { n as f64 * 100.0 / self.n as f64 };
                let level = Support::ALL.get(e).map_or("", |l| l.name());
                writeln!(w, "{}\t{level}\t{n}\t{pct:.2}", s.name())?;
            }
        }
        w.finish().with_context(|| format!("cannot write '{path}'"))
    }
}

fn label_mapped(ctx: &Ctx, summary: &mut Summary, chrom: &str, strand: Strand, s: &Structure) -> GeneLabel {
    let l = GeneLabel::compute(ctx.ann, ctx.params, chrom, strand, s);
    summary.add(ctx.ann, Some(chrom), &l);
    l
}

fn label_unmapped(ctx: &Ctx, summary: &mut Summary) -> GeneLabel {
    let l = GeneLabel::unmapped();
    summary.add(ctx.ann, None, &l);
    l
}

/// Leading columns of the table, the same for every input: which alignment the labels are for.
const ALN_COLUMNS: [(&str, Kind); 6] = [
    ("query_name", Kind::Str),
    ("tp_tag", Kind::Str),
    ("target_name", Kind::Str),
    ("target_start", Kind::U64),
    ("target_end", Kind::U64),
    ("strand", Kind::Str),
];

pub fn begin_aln(out: &mut Sink) -> Result<()> {
    let lead: Vec<(String, Kind)> = ALN_COLUMNS.iter().map(|(n, k)| (n.to_string(), *k)).collect();
    out.begin(&lead)
}

/// Where an alignment lies, for the leading columns.
pub struct Placed<'a> {
    pub chrom: &'a str,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

/// The leading cells. The coordinates are empty for an unmapped alignment, and `tp_tag` when the
/// input has none.
pub fn lead<'a>(query_name: Option<&'a str>, tp_tag: Option<&'a str>, at: Option<&Placed<'a>>) -> [Val<'a>; 6] {
    let s = |x: Option<&'a str>| x.filter(|v| !v.is_empty()).map_or(Val::Null, Val::Str);
    [
        s(query_name),
        s(tp_tag),
        s(at.map(|a| a.chrom)),
        at.map_or(Val::Null, |a| Val::U64(a.start)),
        at.map_or(Val::Null, |a| Val::U64(a.end)),
        at.map_or(Val::Null, |a| Val::Str(if a.strand == Strand::Plus { "+" } else { "-" })),
    ]
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
    /// The cs-tag or CIGAR column the introns are walked from.
    pub introns: String,
    pub walk: Walk,
}

/// Positions of a table's columns in its header.
pub struct TableIdx {
    pub id: usize,
    pub tp: Option<usize>,
    pub chrom: usize,
    pub strand: usize,
    pub start: usize,
    pub end: usize,
    /// The cs or CIGAR column.
    pub introns: usize,
}

impl TableCols {
    /// Find the ID, coordinate and cs (or CIGAR) columns in `header`, and the tp column if present.
    pub fn find(&self, header: &[&str], path: &str) -> Result<TableIdx> {
        let want = [self.id.as_str(), &self.chrom, &self.strand, &self.start, &self.end, &self.introns];
        let find = |name: &str| header.iter().position(|h| *h == name);
        let missing: Vec<&str> = want.iter().copied().filter(|c| find(c).is_none()).collect();
        if !missing.is_empty() {
            bail!(
                "{path}: no column {missing:?}; set the names with --id-col, --chrom-col, --strand-col, --start-col, \
                 --end-col, and --cs-col or --cigar-col"
            );
        }
        let [id, chrom, strand, start, end, introns] = want.map(|c| find(c).unwrap());
        Ok(TableIdx { id, tp: find(&self.tp), chrom, strand, start, end, introns })
    }

    /// Label one table row. A row is unmapped when its strand isn't '+'/'-' or its start or end is
    /// empty; otherwise its coordinates must be integers and its cs or CIGAR must cover the span.
    /// Returns the label and, for a mapped row, where it lies.
    pub fn label_row<'a>(
        &self,
        ctx: &Ctx,
        summary: &mut Summary,
        cells: [Option<&'a str>; 3],
        start: Option<i64>,
        end: Option<i64>,
        raw: [Option<&str>; 2],
    ) -> Result<(GeneLabel, Option<Placed<'a>>)> {
        let [chrom, strand, introns] = cells;
        let strand = strand.and_then(Strand::parse);
        let (Some(strand), Some(chrom)) = (strand, chrom) else { return Ok((label_unmapped(ctx, summary), None)) };
        let (start, end) = match (start, end) {
            (Some(s), Some(e)) if s >= 0 && e > s => (s as u64, e as u64),
            _ if raw.iter().all(|r| r.is_none_or(str::is_empty)) => return Ok((label_unmapped(ctx, summary), None)),
            _ => bail!("start '{}' and end '{}' are not a valid 0-based interval", raw[0].unwrap_or(""), raw[1].unwrap_or("")),
        };
        let s = Structure::walk(start, end, introns.unwrap_or(""), self.walk)?;
        Ok((label_mapped(ctx, summary, chrom, strand, &s), Some(Placed { chrom, start, end, strand })))
    }
}

/// A TSV table with a header line.
pub fn run_tsv(path: &str, cols: &TableCols, ctx: &Ctx, mut out: Sink) -> Result<Summary> {
    let mut summary = Summary::default();
    let mut idx: Option<TableIdx> = None;
    for_each_line(path, |ln, line| {
        let f: Vec<&str> = line.split('\t').collect();
        let Some(ix) = &idx else {
            idx = Some(cols.find(&f, path)?);
            begin_aln(&mut out)?;
            return Ok(());
        };
        if line.is_empty() {
            return Ok(());
        }
        let get = |i: usize| f.get(i).copied();
        let int = |i: usize| get(i).and_then(|s| s.parse::<i64>().ok());
        let (l, at) = cols
            .label_row(
                ctx,
                &mut summary,
                [get(ix.chrom), get(ix.strand), get(ix.introns)],
                int(ix.start),
                int(ix.end),
                [get(ix.start), get(ix.end)],
            )
            .with_context(|| format!("{path}:{ln}"))?;
        out.row(&lead(get(ix.id), ix.tp.and_then(get), at.as_ref()), &l.values(ctx.ann))?;
        Ok(())
    })?;
    if idx.is_none() {
        bail!("{path}: empty table (no header line)");
    }
    out.finish()?;
    Ok(summary)
}

// ---------------------------------------------------------------- PAF

/// minimap2 PAF. Every line is labelled: PAF marks supplementary alignments `tp:A:P` too. Introns
/// come from the cs tag, or from the cg tag (CIGAR) when there is no cs tag.
pub fn run_paf(path: &str, ctx: &Ctx, mut out: Sink) -> Result<Summary> {
    let mut summary = Summary::default();
    begin_aln(&mut out)?;
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
        let tag = |k: &str| f[12..].iter().find_map(|t| t.strip_prefix(k));
        let (tstart, tend) = (num(7)?, num(8)?);
        let (l, at) = match Strand::parse(f[4]) {
            None => (label_unmapped(ctx, &mut summary), None),
            Some(strand) => {
                let s = match (tag("cs:Z:"), tag("cg:Z:")) {
                    (Some(cs), _) => Structure::walk(tstart, tend, cs, Walk::Cs),
                    (None, Some(cg)) => Structure::walk(tstart, tend, cg, Walk::Cigar),
                    (None, None) => bail!("{path}:{ln}: no cs:Z or cg:Z tag to find introns in; run minimap2 with -c or --cs"),
                }
                .with_context(|| format!("{path}:{ln}"))?;
                (label_mapped(ctx, &mut summary, f[5], strand, &s), Some(Placed { chrom: f[5], start: tstart, end: tend, strand }))
            }
        };
        out.row(&lead(Some(f[0]), tag("tp:A:"), at.as_ref()), &l.values(ctx.ann))?;
        Ok(())
    })?;
    out.finish()?;
    Ok(summary)
}
