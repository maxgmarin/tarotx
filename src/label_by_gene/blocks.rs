//! An alignment's structure on the reference: its introns, and its aligned blocks (the span minus
//! the introns). Introns are always walked from the alignment itself: a minimap2 `cs` tag (`~`
//! operations) or a CIGAR (`N` operations, e.g. PAF's cg tag). All coordinates are 0-based half-open.
//!
//! The CIGAR and cs walks are ported from maligno (`cigar_junctions.rs`, `cs_parser.rs`), with
//! their tests. Both also return the reference length walked, which must equal end - start.

use anyhow::{bail, Result};

/// Introns as (start, end) relative to the alignment start, and the reference length walked, from
/// a CIGAR string. `M = X D N` consume the reference; `I S H P` don't.
pub fn introns_from_cigar(cigar: &str) -> Result<(Vec<(u64, u64)>, u64)> {
    let bytes = cigar.as_bytes();
    let (mut pos, mut r_pos) = (0, 0u64);
    let mut out = Vec::new();
    while pos < bytes.len() {
        let start = pos;
        while pos < bytes.len() && bytes[pos].is_ascii_digit() {
            pos += 1;
        }
        if pos == start || pos >= bytes.len() {
            bail!("malformed CIGAR '{cigar}'");
        }
        let n: u64 = cigar[start..pos].parse()?;
        match bytes[pos] {
            b'N' => {
                out.push((r_pos, r_pos + n));
                r_pos += n;
            }
            b'M' | b'=' | b'X' | b'D' => r_pos += n,
            b'I' | b'S' | b'H' | b'P' => {}
            op => bail!("malformed CIGAR '{cigar}': unknown operation '{}'", op as char),
        }
        pos += 1;
    }
    Ok((out, r_pos))
}

/// The same from a minimap2 cs tag: `:N` and `=ACGT` matches, `*xy` substitutions and `-seq`
/// deletions consume the reference, `+seq` insertions don't, and `~xxNNNyy` is an intron of NNN bp.
pub fn introns_from_cs(cs: &str) -> (Vec<(u64, u64)>, u64) {
    let bytes = cs.as_bytes();
    let (mut pos, mut r_pos) = (0, 0u64);
    let mut out = Vec::new();
    let run = |pos: &mut usize, keep: fn(u8) -> bool| {
        let s = *pos;
        while *pos < bytes.len() && keep(bytes[*pos]) {
            *pos += 1;
        }
        *pos - s
    };
    let uint = |pos: &mut usize| {
        let mut n = 0u64;
        while *pos < bytes.len() && bytes[*pos].is_ascii_digit() {
            n = n * 10 + (bytes[*pos] - b'0') as u64;
            *pos += 1;
        }
        n
    };
    while pos < bytes.len() {
        let op = bytes[pos];
        pos += 1;
        match op {
            b':' => r_pos += uint(&mut pos),
            b'=' | b'-' => r_pos += run(&mut pos, |b| b.is_ascii_alphabetic()) as u64,
            b'+' => {
                run(&mut pos, |b| b.is_ascii_alphabetic());
            }
            b'*' => {
                pos += 2;
                r_pos += 1;
            }
            b'~' => {
                pos = (pos + 2).min(bytes.len()); // splice-motif prefix, e.g. "gt"
                let n = uint(&mut pos);
                pos = (pos + 2).min(bytes.len()); // suffix, e.g. "ag"
                out.push((r_pos, r_pos + n));
                r_pos += n;
            }
            _ => {}
        }
    }
    (out, r_pos)
}

/// Which alignment string the introns are walked from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Walk {
    Cs,
    Cigar,
}

/// An alignment's introns and aligned blocks, in absolute coordinates.
#[derive(Debug, PartialEq, Eq)]
pub struct Structure {
    pub blocks: Vec<(u32, u32)>,
    pub introns: Vec<(u32, u32)>,
}

impl Structure {
    /// From the span [start, end) and absolute introns, which must be ascending, non-overlapping
    /// and inside the span. Inputs go through `walk`; this is for the walk and for tests.
    pub fn new(start: u64, end: u64, introns: &[(u64, u64)]) -> Result<Structure> {
        let mut blocks = Vec::with_capacity(introns.len() + 1);
        let mut at = start;
        for &(s, e) in introns {
            if s < at || e <= s || e > end {
                bail!("introns {introns:?} are not ascending, non-overlapping and inside the alignment {start}-{end}");
            }
            if s > at {
                blocks.push((at as u32, s as u32));
            }
            at = e;
        }
        if end > at {
            blocks.push((at as u32, end as u32));
        }
        let introns = introns.iter().map(|&(s, e)| (s as u32, e as u32)).collect();
        Ok(Structure { blocks, introns })
    }

    /// From the span [start, end) and the alignment's cs tag or CIGAR, whose walk must cover the
    /// span exactly.
    pub fn walk(start: u64, end: u64, text: &str, kind: Walk) -> Result<Structure> {
        match kind {
            Walk::Cs => Structure::from_walk(start, end, introns_from_cs(text), "the cs tag"),
            Walk::Cigar => Structure::from_walk(start, end, introns_from_cigar(text)?, "the CIGAR"),
        }
    }

    /// From introns relative to the start (a CIGAR or cs walk), checking the walked length.
    fn from_walk(start: u64, end: u64, (rel, walked): (Vec<(u64, u64)>, u64), what: &str) -> Result<Structure> {
        if start + walked != end {
            bail!("{what} covers {walked} bp of the reference, but the alignment is {start}-{end} ({} bp)", end - start);
        }
        let abs: Vec<(u64, u64)> = rel.iter().map(|&(s, e)| (start + s, start + e)).collect();
        Structure::new(start, end, &abs)
    }

    pub fn aligned_bp(&self) -> u64 {
        self.blocks.iter().map(|&(s, e)| (e - s) as u64).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cigar(c: &str) -> Vec<(u64, u64)> {
        introns_from_cigar(c).unwrap().0
    }

    fn cs(c: &str) -> Vec<(u64, u64)> {
        introns_from_cs(c).0
    }

    // ---- ported from maligno cigar_junctions.rs
    #[test]
    fn cigar_introns() {
        assert!(cigar("").is_empty());
        assert!(cigar("100M").is_empty());
        assert!(cigar("50M5I45M").is_empty());
        assert_eq!(cigar("100M50N100M"), [(100, 150)]);
        assert_eq!(cigar("100M50N100M200N50M"), [(100, 150), (250, 450)]);
        assert_eq!(cigar("100M5D45N50M"), [(105, 150)]); // deletion advances the reference
        assert_eq!(cigar("50M10I50N50M"), [(50, 100)]); // insertion doesn't
        assert_eq!(cigar("10S100M50N100M10S"), [(100, 150)]);
        assert_eq!(cigar("10H100M50N100M10H"), [(100, 150)]);
        assert_eq!(cigar("50=50N50X"), [(50, 100)]);
        assert_eq!(cigar("50M10P50N50M"), [(50, 100)]);
        assert_eq!(introns_from_cigar("100M50N100M").unwrap().1, 250);
        // maligno drops trailing digits silently; here a malformed CIGAR is an error
        assert!(introns_from_cigar("50M50N50M99").is_err());
        assert!(introns_from_cigar("50Q").is_err());
    }

    // ---- ported from maligno cs_parser.rs
    #[test]
    fn cs_introns() {
        assert!(cs("").is_empty());
        assert!(cs(":50*ag:30").is_empty());
        assert_eq!(cs(":100~ct50ag:30"), [(100, 150)]);
        assert_eq!(cs(":50-acg~ct100ag:20"), [(53, 153)]); // deletion advances the reference
        assert_eq!(cs(":50+ac~ct100ag:20"), [(50, 150)]); // insertion doesn't
        assert_eq!(cs(":50*ag~ct100ag:20"), [(51, 151)]); // substitution: one reference base
        assert_eq!(cs(":50~ct100ag:30~gt75ac:20"), [(50, 150), (180, 255)]);
        assert_eq!(cs("=ACGT~ct50ag:20"), [(4, 54)]); // long-form match
        assert_eq!(cs(":50+aaaaa~ct100ag:20"), [(50, 150)]);
        assert_eq!(introns_from_cs(":50~ct100ag:30~gt75ac:20").1, 275);
    }

    #[test]
    fn cs_and_cigar_agree() {
        // the same alignment: 50 match, 2 bp insertion, 30 match, 3 bp deletion, 100 bp intron, 20 match
        let (c, cl) = introns_from_cigar("5S50M2I30M3D100N20M").unwrap();
        let (s, sl) = introns_from_cs(":50+ac:30-tgc~gt100ag:20");
        assert_eq!((c, cl), (s, sl));
    }

    #[test]
    fn blocks_from_introns() {
        let s = Structure::new(100, 1000, &[(200, 300), (400, 900)]).unwrap();
        assert_eq!(s.blocks, [(100, 200), (300, 400), (900, 1000)]);
        assert_eq!(s.aligned_bp(), 300);
        assert_eq!(Structure::new(100, 200, &[]).unwrap().blocks, [(100, 200)]);
        assert!(Structure::new(100, 1000, &[(400, 900), (200, 300)]).is_err()); // not ascending
        assert!(Structure::new(100, 1000, &[(50, 300)]).is_err()); // outside the span
        let w = Structure::walk(1000, 1270, ":50~ct100ag:30~gt75ac:20", Walk::Cs);
        assert!(w.unwrap_err().to_string().contains("covers 275 bp"));
        let w = Structure::walk(1000, 1275, ":50~ct100ag:30~gt75ac:20", Walk::Cs).unwrap();
        assert_eq!(w.introns, [(1050, 1150), (1180, 1255)]);
        assert_eq!(Structure::walk(1000, 1275, "50M100N30M75N20M", Walk::Cigar).unwrap(), w);
        assert!(Structure::walk(1000, 1100, "", Walk::Cs).is_err()); // an empty cs covers nothing
    }
}
