//! Per-alignment labels: what lies at the alignment's 3' end, and the internal-priming call.

use std::io::Write;

use anyhow::Result;

use super::points::{Tract, TxEnd, TxEnds, Tracts};
use crate::common::io::open_output;
use crate::common::strand::Strand;
use super::table::{Kind, Val};

pub struct Params {
    pub polya_max_dist: i64,
    pub tes_max_dist: i64,
}

/// Everything an alignment is labelled against.
pub struct Refs {
    pub tracts: Tracts,
    /// None when no annotation was given: the flag then rests on poly(A) adjacency alone.
    pub tx_ends: Option<TxEnds>,
    pub params: Params,
}

pub struct Label<'a> {
    /// 0-based 3' boundary: target end for '+', target start for '-'.
    pub end3: i64,
    /// Signed distance to the nearest same-strand annotated end, in transcript orientation
    /// (< 0: the read stops short of it), and that end.
    pub tx_end: Option<(i64, &'a TxEnd)>,
    /// Distance to the nearest same-strand tract edge, and that tract.
    pub polya: Option<(i64, &'a Tract)>,
    pub near_polya: bool,
    pub near_tx_end: bool,
    pub flag: bool,
}

impl Refs {
    pub fn label(&self, chrom: &str, strand: &str, start: i64, end: i64) -> Label<'_> {
        let strand = Strand::parse(strand);
        let end3 = if strand == Some(Strand::Plus) { end } else { start };
        let tx_end = strand.and_then(|s| {
            let e = self.tx_ends.as_ref()?.nearest(chrom, s, end3)?;
            let raw = end3 - e.pos;
            Some((if s == Strand::Plus { raw } else { -raw }, e))
        });
        let polya = strand.and_then(|s| {
            let hit = self.tracts.nearest_edge(chrom, s, end3)?;
            Some(((end3 - hit.pos).abs(), hit.tract))
        });
        let near_polya = polya.is_some_and(|(d, _)| d <= self.params.polya_max_dist);
        let near_tx_end = tx_end.is_some_and(|(d, _)| d.abs() <= self.params.tes_max_dist);
        Label { end3, tx_end, polya, near_polya, near_tx_end, flag: near_polya && !near_tx_end }
    }

    /// Names and types of the label columns, in output order.
    pub fn fields(&self) -> Vec<(String, Kind)> {
        let p = &self.params;
        let anno = self.tx_ends.is_some();
        let mut c = vec![("TxAln_3p_end".to_string(), Kind::I64)];
        if anno {
            c.push(("AnnoTranscriptEnd_signed_dist".into(), Kind::I64));
            c.push(("AnnoTranscriptEnd_transcript_id".into(), Kind::Str));
            c.push(("AnnoTranscriptEnd_gene_name".into(), Kind::Str));
        }
        c.push(("GenomicPolyA_edge_dist".into(), Kind::I64));
        c.push(("GenomicPolyA_tract_len".into(), Kind::I64));
        c.push(("GenomicPolyA_score".into(), Kind::I64));
        c.push((format!("GenomicPolyA_within_{}bp", p.polya_max_dist), Kind::Bool));
        if anno {
            c.push((format!("AnnoTranscriptEnd_within_{}bp", p.tes_max_dist), Kind::Bool));
        }
        c.push(("internal_priming_flag".into(), Kind::Bool));
        c
    }

    /// The label values, matching `fields()`; all empty for an unmapped alignment (None).
    pub fn values<'l>(&self, l: Option<&'l Label>) -> Vec<Val<'l>> {
        let Some(l) = l else { return vec![Val::Null; self.fields().len()] };
        let opt_str = |s: &'l Option<String>| s.as_deref().map_or(Val::Null, Val::Str);
        let mut v = vec![Val::I64(l.end3)];
        if self.tx_ends.is_some() {
            match l.tx_end {
                Some((d, e)) => v.extend([Val::I64(d), opt_str(&e.transcript_id), opt_str(&e.gene_name)]),
                None => v.extend([Val::Null, Val::Null, Val::Null]),
            }
        }
        match l.polya {
            Some((d, t)) => v.extend([Val::I64(d), Val::I64((t.end - t.start) as i64), Val::I64(t.score)]),
            None => v.extend([Val::Null, Val::Null, Val::Null]),
        }
        v.push(Val::Bool(l.near_polya));
        if self.tx_ends.is_some() {
            v.push(Val::Bool(l.near_tx_end));
        }
        v.push(Val::Bool(l.flag));
        v
    }
}

/// Counts for the summary printed at the end of a run, and written with --summary.
#[derive(Default)]
pub struct Summary {
    /// Every row: the labelled alignments and the unmapped ones.
    pub n: u64,
    pub unmapped: u64,
    /// Labelled alignments by [near_polya][near_tx_end].
    pub table: [[u64; 2]; 2],
    pub flagged: u64,
    pub off_tracts: u64,
    pub off_annotation: u64,
}

impl Summary {
    /// Count one row: an alignment on `chrom` with its label, or an unmapped one (None).
    pub fn add(&mut self, refs: &Refs, placed: Option<(&str, &Label)>) {
        self.n += 1;
        let Some((chrom, l)) = placed else {
            self.unmapped += 1;
            return;
        };
        self.table[l.near_polya as usize][l.near_tx_end as usize] += 1;
        self.flagged += l.flag as u64;
        if !refs.tracts.has_chrom(chrom) {
            self.off_tracts += 1;
        }
        if refs.tx_ends.as_ref().is_some_and(|t| !t.has_chrom(chrom)) {
            self.off_annotation += 1;
        }
    }

    fn pct(&self, x: u64) -> f64 {
        if self.n == 0 { 0.0 } else { x as f64 * 100.0 / self.n as f64 }
    }

    /// The counts on stderr, as percent of all rows.
    pub fn print(&self, refs: &Refs) {
        let p = &refs.params;
        let cell = |x: u64| format!("{x} ({:.1}%)", self.pct(x));
        let mapped = self.n - self.unmapped;
        let ends = refs.tx_ends.as_ref().map_or(String::new(), |t| format!("{} annotated transcript ends, ", t.n));
        eprintln!("{} alignments ({} unmapped), {ends}{} poly(A) tracts", self.n, self.unmapped, refs.tracts.n);
        for (name, off) in [("poly(A) tracts", Some(self.off_tracts)), ("annotation", refs.tx_ends.as_ref().map(|_| self.off_annotation))] {
            match off {
                Some(off) if mapped > 0 && off == mapped => {
                    eprintln!("WARNING: no alignment is on a chromosome named in the {name} (chr1 vs 1?)")
                }
                Some(off) if off as f64 > 0.01 * mapped as f64 => {
                    eprintln!("note: {:.1}% of alignments are on chromosomes absent from the {name}", off as f64 * 100.0 / mapped as f64)
                }
                _ => {}
            }
        }
        if refs.tx_ends.is_some() {
            let t = &self.table;
            let w = 18;
            eprintln!("{:>32}  {:>w$}  {:>w$}  {:>w$}", format!("AnnoTranscriptEnd_within_{}bp:", p.tes_max_dist), "False", "True", "total");
            for (k, name) in [(0, "False"), (1, "True")] {
                eprintln!(
                    "{:>32}  {:>w$}  {:>w$}  {:>w$}",
                    format!("GenomicPolyA_within_{}bp {name}", p.polya_max_dist),
                    cell(t[k][0]),
                    cell(t[k][1]),
                    cell(t[k][0] + t[k][1])
                );
            }
            eprintln!("{:>32}  {:>w$}  {:>w$}  {:>w$}", "total", cell(t[0][0] + t[1][0]), cell(t[0][1] + t[1][1]), cell(mapped));
        } else {
            eprintln!(
                "GenomicPolyA_within_{}bp: {} (no annotation given, so the flag is poly(A) adjacency alone)",
                p.polya_max_dist,
                cell(self.table[1][0])
            );
        }
        eprintln!("internal_priming_flag: {}", cell(self.flagged));
        eprintln!("unmapped (no labels): {}", cell(self.unmapped));
    }

    /// The counts as a TSV: one row per combination of the label columns that can occur, zeros
    /// included, then the unmapped alignments (empty label cells). Percent is of all rows.
    pub fn write_tsv(&self, path: &str, refs: &Refs) -> Result<()> {
        let p = &refs.params;
        let tf = |b: bool| if b { "True" } else { "False" };
        let mut w = open_output(path)?;
        let anno = refs.tx_ends.is_some();
        let mut head = vec![format!("GenomicPolyA_within_{}bp", p.polya_max_dist)];
        if anno {
            head.push(format!("AnnoTranscriptEnd_within_{}bp", p.tes_max_dist));
        }
        head.extend(["internal_priming_flag", "count", "percent"].map(String::from));
        writeln!(w, "{}", head.join("\t"))?;
        let mut row = |cells: Vec<&str>, n: u64| writeln!(w, "{}\t{n}\t{:.2}", cells.join("\t"), self.pct(n));
        for near_polya in [false, true] {
            if anno {
                for near_tx_end in [false, true] {
                    let flag = near_polya && !near_tx_end;
                    row(vec![tf(near_polya), tf(near_tx_end), tf(flag)], self.table[near_polya as usize][near_tx_end as usize])?;
                }
            } else {
                row(vec![tf(near_polya), tf(near_polya)], self.table[near_polya as usize][0])?;
            }
        }
        row(vec![""; if anno { 3 } else { 2 }], self.unmapped)?;
        w.finish()?;
        Ok(())
    }
}
