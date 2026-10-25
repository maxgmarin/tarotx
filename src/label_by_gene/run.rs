//! Labelling every input alignment, and the summary counts.

use std::io::Write;

use anyhow::{bail, Context, Result};

use super::annot::Annotation;
use super::blocks::Structure;
use super::label::{GeneLabel, Params, Status, Support};
use crate::common::input::{for_each_alignment, table_columns, Aln, InputFormat, TableCols};
use crate::common::io::open_output;
use crate::common::table::Sink;

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
        eprintln!(
            "{} alignments ({} unmapped), {} genes, {} annotated introns",
            self.n,
            self.n - self.mapped,
            ann.genes.len(),
            ann.n_introns
        );
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

/// Label every alignment of `path` and write the table to `output`. Returns the summary counts.
pub fn run(path: &str, format: InputFormat, cols: &TableCols, ann: &Annotation, params: &Params, output: &str) -> Result<Summary> {
    let mut out = Sink::open(output, &table_columns(GeneLabel::fields()), "label-by-gene")?;
    let mut summary = Summary::default();
    for_each_alignment(path, format, cols, |a| {
        let l = label(a, ann, params)?;
        summary.add(ann, a.at.as_ref().map(|p| p.chrom), &l);
        let mut cells = a.lead_cells();
        cells.extend(l.values(ann));
        out.row(&cells)
    })?;
    out.finish()?;
    Ok(summary)
}

/// The gene label of one alignment. A mapped alignment needs its cs tag or CIGAR, whose walk must
/// cover its span.
fn label(a: &Aln, ann: &Annotation, params: &Params) -> Result<GeneLabel> {
    let Some(at) = &a.at else { return Ok(GeneLabel::unmapped()) };
    let Some((text, walk)) = a.introns else {
        bail!("no cs:Z or cg:Z tag to find introns in; run minimap2 with -c or --cs")
    };
    let s = Structure::walk(at.start, at.end, text, walk)?;
    Ok(GeneLabel::compute(ann, params, at.chrom, at.strand, &s))
}
