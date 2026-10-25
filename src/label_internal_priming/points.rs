//! Reference points that alignment 3' ends are compared against: poly(A) tract edges and
//! annotated transcript 3' ends, kept per chromosome and strand, sorted by position.
//! All coordinates are 0-based half-open.

use std::collections::HashMap;

use anyhow::{bail, Context, Result};

use crate::common::io::for_each_line;
use crate::common::strand::Strand;

/// Index of a strand in the per-chromosome [plus, minus] pairs below.
fn idx(strand: Strand) -> usize {
    strand as usize
}

/// Index of the point nearest to `x` among `n` points sorted by `pos`. When two points are
/// equally close, the one beyond the 3' end wins (the higher one on '+', the lower one on '-').
pub fn nearest(n: usize, pos: impl Fn(usize) -> i64, x: i64, strand: Strand) -> Option<usize> {
    if n == 0 {
        return None;
    }
    // first point at or after x
    let (mut lo, mut hi) = (0, n);
    while lo < hi {
        let mid = (lo + hi) / 2;
        if pos(mid) < x {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    let right = lo.min(n - 1);
    let left = lo.saturating_sub(1);
    let (d_left, d_right) = ((x - pos(left)).abs(), (pos(right) - x).abs());
    Some(if d_right < d_left || (d_right == d_left && strand == Strand::Plus) { right } else { left })
}

// ---------------------------------------------------------------- poly(A) tracts

#[derive(Clone, Copy, Debug)]
pub struct Tract {
    pub start: u32,
    pub end: u32,
    pub score: i64,
}

/// Tracts of one chromosome and strand, sorted by start.
#[derive(Default)]
struct StrandTracts {
    tracts: Vec<Tract>,
    /// Edge positions and their tract, sorted. Only built when tracts overlap; otherwise the edges
    /// are already in order as start0, end0, start1, end1, ... and are read from `tracts`.
    edges: Option<Vec<(u32, u32)>>,
}

impl StrandTracts {
    fn n_edges(&self) -> usize {
        self.edges.as_ref().map_or(2 * self.tracts.len(), Vec::len)
    }
    fn edge(&self, k: usize) -> (i64, &Tract) {
        match &self.edges {
            Some(e) => (e[k].0 as i64, &self.tracts[e[k].1 as usize]),
            None => {
                let t = &self.tracts[k / 2];
                (if k.is_multiple_of(2) { t.start } else { t.end } as i64, t)
            }
        }
    }
}

/// The nearest tract edge to a 3' end.
pub struct EdgeHit<'a> {
    pub pos: i64,
    pub tract: &'a Tract,
}

#[derive(Default)]
pub struct Tracts {
    by_chrom: HashMap<String, [StrandTracts; 2]>,
    pub n: usize,
}

impl Tracts {
    /// polyscan poly(A) tracts: BED6 (`chrom start end name score strand`) or 5-column
    /// `bedtools merge -s -c 5,6 -o max,distinct` output (`chrom start end score strand`).
    pub fn load(path: &str) -> Result<Tracts> {
        let mut t = Tracts::default();
        let mut ncols = None;
        for_each_line(path, |ln, line| {
            if line.is_empty() || line.starts_with('#') {
                return Ok(());
            }
            let f: Vec<&str> = line.split('\t').collect();
            if *ncols.get_or_insert(f.len()) != f.len() {
                bail!("{path}:{ln}: expected {} columns like the first line, found {}", ncols.unwrap(), f.len());
            }
            let (score, strand) = match f.len() {
                6 => (f[4], f[5]),
                5 => (f[3], f[4]),
                n => bail!("{path}: expected 5 or 6 columns, found {n}"),
            };
            let Some(strand) = Strand::parse(strand) else {
                bail!("{path}:{ln}: strand column contains '{strand}', expected '+' or '-'");
            };
            let num = |s: &str, what: &str| {
                s.parse::<u32>().with_context(|| format!("{path}:{ln}: bad {what} '{s}'"))
            };
            let tract = Tract {
                start: num(f[1], "start")?,
                end: num(f[2], "end")?,
                score: score.parse().with_context(|| format!("{path}:{ln}: bad score '{score}'"))?,
            };
            if !t.by_chrom.contains_key(f[0]) {
                t.by_chrom.insert(f[0].to_string(), Default::default());
            }
            t.by_chrom.get_mut(f[0]).unwrap()[idx(strand)].tracts.push(tract);
            t.n += 1;
            Ok(())
        })?;
        for st in t.by_chrom.values_mut().flatten() {
            st.tracts.sort_by_key(|t| (t.start, t.end));
            let interleaved = st.tracts.iter().all(|t| t.start <= t.end)
                && st.tracts.windows(2).all(|w| w[0].end <= w[1].start);
            if !interleaved {
                let mut e: Vec<(u32, u32)> = (0..st.tracts.len() as u32)
                    .flat_map(|i| [(st.tracts[i as usize].start, i), (st.tracts[i as usize].end, i)])
                    .collect();
                e.sort_by_key(|&(p, _)| p);
                st.edges = Some(e);
            }
        }
        Ok(t)
    }

    pub fn has_chrom(&self, chrom: &str) -> bool {
        self.by_chrom.contains_key(chrom)
    }

    /// The same-strand tract edge (start or end) nearest to `x`.
    pub fn nearest_edge(&self, chrom: &str, strand: Strand, x: i64) -> Option<EdgeHit<'_>> {
        let st = &self.by_chrom.get(chrom)?[idx(strand)];
        let k = nearest(st.n_edges(), |k| st.edge(k).0, x, strand)?;
        let (pos, tract) = st.edge(k);
        Some(EdgeHit { pos, tract })
    }
}

// ---------------------------------------------------------------- annotated transcript ends

pub struct TxEnd {
    pub pos: i64,
    pub transcript_id: Option<String>,
    pub gene_name: Option<String>,
}

#[derive(Default)]
pub struct TxEnds {
    by_chrom: HashMap<String, [Vec<TxEnd>; 2]>,
    pub n: usize,
}

/// Value of `key` in a GFF3 attribute column (`k1=v1;k2=v2`).
fn attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attrs.split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='))
}

/// True when the comma-separated `tag` attribute contains exactly `tag`.
fn has_tag(attrs: &str, tag: &str) -> bool {
    attr(attrs, "tag").is_some_and(|tags| tags.split(',').any(|t| t == tag))
}

impl TxEnds {
    /// One point per GFF3 `transcript` 3' end (end for '+', start for '-').
    /// `basic_only` keeps transcripts tagged `basic`; `drop_end_nf` drops those tagged `mRNA_end_NF`.
    /// A GFF3 without transcript lines (e.g. a gene/exon-only one), or with none left after these
    /// filters, is an error.
    pub fn load(path: &str, basic_only: bool, drop_end_nf: bool) -> Result<TxEnds> {
        let mut t = TxEnds::default();
        let mut in_fasta = false;
        let mut seen = 0u64;
        for_each_line(path, |ln, line| {
            if in_fasta || line.starts_with('#') || line.is_empty() {
                in_fasta |= line.starts_with("##FASTA");
                return Ok(());
            }
            let f: Vec<&str> = line.splitn(9, '\t').collect();
            if f.len() < 9 {
                bail!("{path}:{ln}: expected 9 GFF3 columns, found {}", f.len());
            }
            if f[2] != "transcript" {
                return Ok(());
            }
            seen += 1;
            let Some(strand) = Strand::parse(f[6]) else { return Ok(()) };
            let attrs = f[8];
            if (basic_only && !has_tag(attrs, "basic")) || (drop_end_nf && has_tag(attrs, "mRNA_end_NF")) {
                return Ok(());
            }
            let num = |s: &str, what: &str| {
                s.parse::<i64>().with_context(|| format!("{path}:{ln}: bad {what} '{s}'"))
            };
            let pos = match strand {
                Strand::Plus => num(f[4], "end")?,
                Strand::Minus => num(f[3], "start")? - 1,
            };
            let end = TxEnd {
                pos,
                transcript_id: attr(attrs, "transcript_id").map(str::to_string),
                gene_name: attr(attrs, "gene_name").map(str::to_string),
            };
            if !t.by_chrom.contains_key(f[0]) {
                t.by_chrom.insert(f[0].to_string(), Default::default());
            }
            t.by_chrom.get_mut(f[0]).unwrap()[idx(strand)].push(end);
            t.n += 1;
            Ok(())
        })?;
        if seen == 0 {
            bail!("{path} has no transcript lines; label-internal-priming needs the full GENCODE GFF3 (a gene/exon-only GFF3 has no transcripts)");
        }
        if t.n == 0 {
            bail!("{path}: none of its {seen} transcript lines is left after --gencode-basic-only and dropping mRNA_end_NF transcripts");
        }
        for v in t.by_chrom.values_mut().flatten() {
            v.sort_by_key(|e| e.pos); // stable: ties keep file order
        }
        Ok(t)
    }

    pub fn has_chrom(&self, chrom: &str) -> bool {
        self.by_chrom.contains_key(chrom)
    }

    /// The same-strand annotated transcript end nearest to `x`.
    pub fn nearest(&self, chrom: &str, strand: Strand, x: i64) -> Option<&TxEnd> {
        let v = &self.by_chrom.get(chrom)?[idx(strand)];
        nearest(v.len(), |i| v[i].pos, x, strand).map(|i| &v[i])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nearest_ties_go_beyond_the_3p_end() {
        let p = [10i64, 20, 30];
        let at = |i: usize| p[i];
        assert_eq!(nearest(3, at, 15, Strand::Plus), Some(1));
        assert_eq!(nearest(3, at, 15, Strand::Minus), Some(0));
        assert_eq!(nearest(3, at, 20, Strand::Minus), Some(1));
        assert_eq!(nearest(3, at, 0, Strand::Minus), Some(0));
        assert_eq!(nearest(3, at, 99, Strand::Plus), Some(2));
        assert_eq!(nearest(3, at, 24, Strand::Plus), Some(1));
        assert_eq!(nearest(3, at, 26, Strand::Minus), Some(2));
        assert_eq!(nearest(0, at, 5, Strand::Plus), None);
    }

    #[test]
    fn overlapping_tracts_use_sorted_edges() {
        let p = std::env::temp_dir().join(format!("tarotx-overlap-{}.bed", std::process::id()));
        // unmerged: [100,150) overlaps [120,130), so edges are 100,120,130,150
        std::fs::write(&p, "c\t100\t150\tA\t90\t+\nc\t120\t130\tA\t80\t+\n").unwrap();
        let t = Tracts::load(p.to_str().unwrap()).unwrap();
        let hit = |x| t.nearest_edge("c", Strand::Plus, x).map(|h| (h.pos, h.tract.score));
        assert_eq!(hit(131), Some((130, 80)));
        assert_eq!(hit(148), Some((150, 90)));
        assert_eq!(hit(111), Some((120, 80)));
        assert_eq!(t.nearest_edge("c", Strand::Minus, 100).map(|h| h.pos), None);
        std::fs::remove_file(p).unwrap();
    }

    #[test]
    fn gff_attributes() {
        let a = "ID=x;transcript_id=ENST1;gene_name=ABC;tag=basic,mRNA_end_NF";
        assert_eq!(attr(a, "transcript_id"), Some("ENST1"));
        assert_eq!(attr(a, "ID"), Some("x"));
        assert_eq!(attr(a, "gene_id"), None);
        assert!(has_tag(a, "basic") && has_tag(a, "mRNA_end_NF") && !has_tag(a, "mRNA"));
    }
}
