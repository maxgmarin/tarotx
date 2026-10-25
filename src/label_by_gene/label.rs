//! The gene assignment of one alignment: which gene it belongs to and how clear that choice is
//! (`GeneAssignment_Status`), and how the gene is supported (`GeneAssignment_SupportLevel`). The rules
//! are in the README ("How classification works").

use super::annot::Annotation;
use crate::common::strand::Strand;
use super::blocks::Structure;
use crate::common::table::{float, Kind, Val};

pub struct Params {
    /// Share of the read's aligned bp on a gene's exons that makes its support `exonic`.
    pub min_exonic_frac: f64,
    /// How far ahead in exonic fraction the top gene must be, at equal matched introns.
    pub tie_margin: f64,
    pub near_max_dist: u64,
}

/// How a gene is supported, weakest first (the order is used for ranking).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Support {
    IntronicOnly,
    Exonic,
    KnownJunctions,
}

impl Support {
    pub const ALL: [Support; 3] = [Support::KnownJunctions, Support::Exonic, Support::IntronicOnly];

    pub fn name(self) -> &'static str {
        match self {
            Support::KnownJunctions => "known_junctions",
            Support::Exonic => "exonic",
            Support::IntronicOnly => "intronic_only",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    UniqueGene,
    BestOfMultiple,
    Tied,
    NearbyGene,
    AntisenseOnly,
    Intergenic,
    Unmapped,
}

impl Status {
    pub const ALL: [Status; 7] =
        [Status::UniqueGene, Status::BestOfMultiple, Status::Tied, Status::NearbyGene, Status::AntisenseOnly, Status::Intergenic, Status::Unmapped];

    pub fn name(self) -> &'static str {
        match self {
            Status::UniqueGene => "unique_gene",
            Status::BestOfMultiple => "best_of_multiple",
            Status::Tied => "tied",
            Status::NearbyGene => "nearby_gene",
            Status::AntisenseOnly => "antisense_only",
            Status::Intergenic => "intergenic",
            Status::Unmapped => "unmapped",
        }
    }
}

/// A same-strand gene the read touches, with its support.
#[derive(Clone, Copy, Debug)]
struct Candidate {
    gene: usize,
    matched: u64,
    exonic_bp: u64,
    level: Support,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    /// The read lies past the gene's 3' end.
    DownstreamOf3p,
    /// The read lies before the gene's 5' end.
    UpstreamOf5p,
}

impl Side {
    fn name(self) -> &'static str {
        match self {
            Side::DownstreamOf3p => "downstream_of_3p",
            Side::UpstreamOf5p => "upstream_of_5p",
        }
    }
}

#[derive(Debug)]
pub struct GeneLabel {
    pub status: Status,
    pub support: Option<Support>,
    /// The chosen gene: the top-ranked same-strand gene, or the nearest one for nearby_gene.
    pub gene: Option<usize>,
    read_introns: u64,
    /// None for unmapped alignments.
    aligned_bp: Option<u64>,
    /// Same-strand touched genes, in rank order. The first `n_best` are the best gene(s).
    candidates: Vec<Candidate>,
    /// Genes tied for the best call: 1 for unique_gene and best_of_multiple, 2+ for tied, 0 without
    /// a gene.
    n_best: usize,
    near: Option<(u64, Side)>,
    /// GeneAssignment_Candidates, GeneAssignment_AntisenseGenes
    text: [String; 2],
}

/// Bases of the sorted blocks on the sorted exon intervals.
fn overlap_bp(blocks: &[(u32, u32)], exons: &[(u32, u32)]) -> u64 {
    let (mut i, mut j, mut bp) = (0, 0, 0u64);
    while i < blocks.len() && j < exons.len() {
        let (s, e) = (blocks[i].0.max(exons[j].0), blocks[i].1.min(exons[j].1));
        if e > s {
            bp += (e - s) as u64;
        }
        if blocks[i].1 < exons[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    bp
}

/// True when the gene body [start, end) overlaps one of the sorted blocks.
fn touches_blocks(blocks: &[(u32, u32)], start: u32, end: u32) -> bool {
    let k = blocks.partition_point(|b| b.1 <= start);
    blocks.get(k).is_some_and(|b| b.0 < end)
}

/// A gene's name, or its ID when it has none.
fn name(ann: &Annotation, g: usize) -> &str {
    let gene = &ann.genes[g];
    gene.name.as_deref().unwrap_or(&gene.id)
}

/// Rounded to 3 decimals as numpy does (ties to even).
fn round3(x: f64) -> f64 {
    (x * 1000.0).round_ties_even() / 1000.0
}

/// The nearest same-strand gene to the blocks, among `genes`: (gene, distance, side). The distance
/// is the gap to the nearest aligned block; equal distances go to the smaller gene_id.
fn nearest(ann: &Annotation, blocks: &[(u32, u32)], genes: impl Iterator<Item = usize>) -> Option<(usize, u64, Side)> {
    let mut best: Option<(u64, &str, usize, Side)> = None;
    for g in genes {
        let gene = &ann.genes[g];
        // the closest block before the gene, and after it
        let k = blocks.partition_point(|b| b.1 <= gene.start);
        let before = k.checked_sub(1).map(|i| (gene.start - blocks[i].1) as u64);
        let after = blocks.get(k).map(|b| (b.0 - gene.end) as u64);
        let (d, read_after) = match (before, after) {
            (Some(b), Some(a)) if a < b => (a, true),
            (Some(b), _) => (b, false),
            (None, Some(a)) => (a, true),
            (None, None) => continue,
        };
        // a block after the gene: past its 3' end on '+', before its 5' end on '-'
        let side = match (gene.strand == Strand::Plus, read_after) {
            (true, true) | (false, false) => Side::DownstreamOf3p,
            _ => Side::UpstreamOf5p,
        };
        if best.is_none_or(|b| (d, gene.id.as_str()) < (b.0, b.1)) {
            best = Some((d, gene.id.as_str(), g, side));
        }
    }
    best.map(|(d, _, g, side)| (g, d, side))
}

impl GeneLabel {
    pub fn unmapped() -> GeneLabel {
        GeneLabel {
            status: Status::Unmapped,
            support: None,
            gene: None,
            read_introns: 0,
            aligned_bp: None,
            candidates: Vec::new(),
            n_best: 0,
            near: None,
            text: Default::default(),
        }
    }

    pub fn compute(ann: &Annotation, p: &Params, chrom: &str, strand: Strand, s: &Structure) -> GeneLabel {
        let blocks = &s.blocks;
        let aligned_bp = s.aligned_bp();
        let mut label = GeneLabel { read_introns: s.introns.len() as u64, aligned_bp: Some(aligned_bp), ..GeneLabel::unmapped() };
        let (Some(first), Some(last)) = (blocks.first(), blocks.last()) else {
            label.status = Status::Intergenic;
            return label;
        };
        let (span_start, span_end) = (first.0, last.1);

        // 1. Touched genes: their body overlaps an aligned block.
        let (mut same, mut antisense, mut not_touched) = (Vec::new(), Vec::new(), Vec::new());
        for g in ann.overlapping(chrom, span_start, span_end) {
            let gene = &ann.genes[g];
            match (touches_blocks(blocks, gene.start, gene.end), gene.strand == strand) {
                (true, true) => same.push(g),
                (true, false) => antisense.push(g),
                (false, true) => not_touched.push(g), // inside one of the read's introns
                (false, false) => {}
            }
        }

        // 2-3. Support for each same-strand gene.
        let frac = |bp: u64| if aligned_bp > 0 { bp as f64 / aligned_bp as f64 } else { 0.0 };
        let mut cands: Vec<Candidate> = same
            .iter()
            .map(|&g| {
                let exonic_bp = overlap_bp(blocks, ann.exons(g));
                let annotated = ann.introns(g);
                let matched = s.introns.iter().filter(|i| annotated.binary_search(i).is_ok()).count() as u64;
                let level = if matched > 0 {
                    Support::KnownJunctions
                } else if frac(exonic_bp) >= p.min_exonic_frac {
                    Support::Exonic
                } else {
                    Support::IntronicOnly
                };
                Candidate { gene: g, matched, exonic_bp, level }
            })
            .collect();

        // 4. Rank: support level, matched introns, exonic bp, then gene_id to fix the order.
        cands.sort_by(|a, b| {
            (b.level, b.matched, b.exonic_bp)
                .cmp(&(a.level, a.matched, a.exonic_bp))
                .then_with(|| ann.genes[a.gene].id.cmp(&ann.genes[b.gene].id))
        });

        if let Some(top) = cands.first().copied() {
            // 5. The best gene(s): the top gene, and every gene tied with it, i.e. at the same level
            // with the same matched introns and less than the margin behind in exonic fraction.
            // Intronic-only genes all tie with each other.
            // The competitors: every gene with exonic support (either exonic level), or, when there
            // is none, every touched gene. A gene only grazed intronically doesn't compete with an
            // exonic one.
            let competing = if top.level == Support::IntronicOnly {
                cands.len()
            } else {
                cands.iter().filter(|c| c.level != Support::IntronicOnly).count()
            };
            label.n_best = cands
                .iter()
                .take_while(|c| {
                    // an exact tie always counts, so --tie-margin 0 means exact ties only
                    let close = c.exonic_bp == top.exonic_bp || frac(top.exonic_bp) - frac(c.exonic_bp) < p.tie_margin;
                    c.level == top.level && (top.level == Support::IntronicOnly || (c.matched == top.matched && close))
                })
                .count();
            // 6. Status: one competitor, a clear winner, or a tie.
            label.status = match (competing, label.n_best) {
                (1, _) => Status::UniqueGene,
                (_, 1) => Status::BestOfMultiple,
                _ => Status::Tied,
            };
            label.support = Some(top.level);
            label.gene = Some(top.gene);
        } else {
            // 7. No same-strand gene touched: nearby_gene, then antisense_only, then intergenic.
            let flanking = ann.flanking(chrom, strand, span_start, span_end);
            match nearest(ann, blocks, flanking.into_iter().chain(not_touched)).filter(|n| n.1 <= p.near_max_dist) {
                Some((g, d, side)) => {
                    label.status = Status::NearbyGene;
                    label.gene = Some(g);
                    label.near = Some((d, side));
                }
                None if !antisense.is_empty() => label.status = Status::AntisenseOnly,
                None => label.status = Status::Intergenic,
            }
        }

        let n_introns = label.read_introns;
        label.text = [
            cands
                .iter()
                .map(|c| {
                    let f = float(round3(frac(c.exonic_bp)));
                    format!("{} ({}) {}/{n_introns} {f}", name(ann, c.gene), ann.genes[c.gene].id, c.matched)
                })
                .collect::<Vec<_>>()
                .join("; "),
            antisense.iter().map(|&g| format!("{} ({})", name(ann, g), ann.genes[g].id)).collect::<Vec<_>>().join("; "),
        ];
        label.candidates = cands;
        label
    }

    /// Names and types of the label columns, in output order.
    pub fn fields() -> Vec<(String, Kind)> {
        use Kind::*;
        [
            ("GeneAssignment_Status", Str),
            ("GeneAssignment_SupportLevel", Str),
            ("GeneAssignment_GeneID", Str),
            ("GeneAssignment_GeneName", Str),
            ("GeneAssignment_GeneType", Str),
            ("GeneAssignment_IntronsMatched", I64),
            ("GeneAssignment_ReadIntrons", I64),
            ("GeneAssignment_ExonicFrac", F64),
            ("GeneAssignment_AlignedBp", I64),
            ("GeneAssignment_NumSameStrandGenes", I64),
            ("GeneAssignment_NumBestGenes", I64),
            ("GeneAssignment_Candidates", Str),
            ("GeneAssignment_AntisenseGenes", Str),
            ("GeneAssignment_NearbyGeneDist", I64),
            ("GeneAssignment_NearbyGeneSide", Str),
        ]
        .into_iter()
        .map(|(n, k)| (n.to_string(), k))
        .collect()
    }

    /// The label values, matching `fields()`. Empty means not applicable.
    pub fn values<'a>(&'a self, ann: &'a Annotation) -> Vec<Val<'a>> {
        let opt = |x: Option<&'a str>| x.map_or(Val::Null, Val::Str);
        let text = |x: &'a String| if x.is_empty() { Val::Null } else { Val::Str(x) };
        let int = |x: Option<u64>| x.map_or(Val::Null, |v| Val::I64(v as i64));
        let top = self.candidates.first().filter(|_| self.support.is_some());
        let g = self.gene.map(|g| &ann.genes[g]);
        let frac = top.zip(self.aligned_bp.filter(|&a| a > 0)).map(|(t, a)| round3(t.exonic_bp as f64 / a as f64));
        vec![
            Val::Str(self.status.name()),
            opt(self.support.map(Support::name)),
            opt(g.map(|g| g.id.as_str())),
            opt(self.gene.map(|i| name(ann, i))),
            opt(g.and_then(|g| g.gene_type.as_deref())),
            int(top.map(|t| t.matched)),
            Val::I64(self.read_introns as i64),
            frac.map_or(Val::Null, Val::F64),
            int(self.aligned_bp),
            Val::I64(self.candidates.len() as i64),
            Val::I64(self.n_best as i64),
            text(&self.text[0]),
            text(&self.text[1]),
            int(self.near.map(|n| n.0)),
            opt(self.near.map(|n| n.1.name())),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::label_by_gene::annot::tests::{load, load_opts, GFF};
    use crate::label_by_gene::annot::Options;
    use Support::*;

    const P: Params = Params { min_exonic_frac: 0.2, tie_margin: 0.1, near_max_dist: 5000 };

    /// Extra genes: G8 and G9 share part of G1's last exon (for the margin), G10 has a long intron
    /// and no gene_name.
    const EXTRA: &str = "chrT\tt\tgene\t101\t1000\t.\t+\t.\tgene_id=G8;gene_name=EIGHT;gene_type=lncRNA
chrT\tt\texon\t901\t950\t.\t+\t.\tParent=T8;gene_id=G8
chrT\tt\tgene\t101\t1000\t.\t+\t.\tgene_id=G9;gene_name=NINE;gene_type=lncRNA
chrT\tt\texon\t901\t985\t.\t+\t.\tParent=T9;gene_id=G9
chrT\tt\tgene\t5001\t8000\t.\t+\t.\tgene_id=G10;gene_type=lncRNA
chrT\tt\texon\t5001\t5100\t.\t+\t.\tParent=T10;gene_id=G10
chrT\tt\texon\t7901\t8000\t.\t+\t.\tParent=T10;gene_id=G10
";

    /// Label a read given as its blocks.
    fn label(ann: &Annotation, p: &Params, strand: Strand, blocks: &[(u64, u64)]) -> GeneLabel {
        let introns: Vec<(u64, u64)> = blocks.windows(2).map(|w| (w[0].1, w[1].0)).collect();
        let s = Structure::new(blocks[0].0, blocks.last().unwrap().1, &introns).unwrap();
        GeneLabel::compute(ann, p, "chrT", strand, &s)
    }

    /// Status, support and chosen gene ID, for compact asserts.
    fn sum(ann: &Annotation, l: &GeneLabel) -> (Status, Option<Support>, Option<String>) {
        (l.status, l.support, l.gene.map(|g| ann.genes[g].id.clone()))
    }

    fn col<'a>(ann: &'a Annotation, l: &'a GeneLabel, c: &str) -> Val<'a> {
        let i = GeneLabel::fields().iter().position(|(n, _)| n == c).unwrap();
        l.values(ann)[i]
    }

    /// The tied genes: the first NumBestGenes candidates of a tied call (empty otherwise).
    fn tied(ann: &Annotation, l: &GeneLabel) -> String {
        let n = if l.status == Status::Tied { l.n_best } else { 0 };
        l.candidates[..n].iter().map(|c| ann.genes[c.gene].id.as_str()).collect::<Vec<_>>().join(",")
    }

    fn no_readthrough(text: &str) -> Annotation {
        load_opts(text, &Options { gene_types: &[], skip_readthrough: true })
    }

    #[test]
    fn known_junction_tie_and_best() {
        let a = load(GFF, &[]);
        // intron (200,300) is annotated in G1, G3 (pseudogene) and G4 (readthrough): a three-way tie
        let l = label(&a, &P, Strand::Plus, &[(150, 200), (300, 400)]);
        assert_eq!(sum(&a, &l), (Status::Tied, Some(KnownJunctions), Some("G1".into())));
        assert_eq!((col(&a, &l, "GeneAssignment_NumSameStrandGenes"), col(&a, &l, "GeneAssignment_NumBestGenes")), (Val::I64(3), Val::I64(3)));
        assert_eq!(tied(&a, &l), "G1,G3,G4");
        assert!(matches!(col(&a, &l, "GeneAssignment_Candidates"), Val::Str(c) if c.starts_with("ONE (G1) 1/1 1.0; THREE (G3) 1/1 1.0; ONE-RT (G4)")));
        // both introns are in G1 (and readthrough G4), only the first in G3: without G4, G1 is
        // best_of_multiple
        let b = no_readthrough(GFF);
        let l = label(&b, &P, Strand::Plus, &[(150, 200), (300, 400), (900, 990)]);
        assert_eq!(sum(&b, &l), (Status::BestOfMultiple, Some(KnownJunctions), Some("G1".into())));
        assert_eq!((col(&b, &l, "GeneAssignment_NumSameStrandGenes"), col(&b, &l, "GeneAssignment_NumBestGenes")), (Val::I64(2), Val::I64(1)));
        assert_eq!(col(&b, &l, "GeneAssignment_IntronsMatched"), Val::I64(2));
        assert_eq!(col(&b, &l, "GeneAssignment_Candidates"), Val::Str("ONE (G1) 2/2 1.0; THREE (G3) 1/2 0.625"));
        assert_eq!(tied(&b, &l), "");
    }

    #[test]
    fn unique_gene_among_competitors() {
        let a = load(GFF, &[]);
        // intron (250,400) is only in G1; G3 and G4 are touched but with no exonic bp (intronic_only)
        let l = label(&a, &P, Strand::Plus, &[(200, 250), (400, 440)]);
        assert_eq!(sum(&a, &l), (Status::UniqueGene, Some(KnownJunctions), Some("G1".into())));
        assert_eq!((col(&a, &l, "GeneAssignment_NumSameStrandGenes"), col(&a, &l, "GeneAssignment_NumBestGenes")), (Val::I64(3), Val::I64(1)));
        assert_eq!(col(&a, &l, "GeneAssignment_ExonicFrac"), Val::F64(1.0));
        assert_eq!(col(&a, &l, "GeneAssignment_AlignedBp"), Val::I64(90));
        // G1 matches both introns; G4 matches none but has 50 of 150 bp on its exons (exonic):
        // a competitor, so G1 is best_of_multiple, not unique_gene
        let l = label(&a, &P, Strand::Plus, &[(200, 250), (400, 450), (900, 950)]);
        assert_eq!(sum(&a, &l), (Status::BestOfMultiple, Some(KnownJunctions), Some("G1".into())));
        assert_eq!(col(&a, &l, "GeneAssignment_Candidates"), Val::Str("ONE (G1) 2/2 1.0; ONE-RT (G4) 0/2 0.333; THREE (G3) 0/2 0.0"));
    }

    #[test]
    fn exonic_margin() {
        let a = no_readthrough(&format!("{GFF}{EXTRA}"));
        // on [910,990): G1 1.0, G9 0.9375 (within 0.10: tied), G8 0.5; G3 has no exon there
        let l = label(&a, &P, Strand::Plus, &[(910, 990)]);
        assert_eq!(sum(&a, &l), (Status::Tied, Some(Exonic), Some("G1".into())));
        assert_eq!(tied(&a, &l), "G1,G9");
        assert_eq!(col(&a, &l, "GeneAssignment_NumBestGenes"), Val::I64(2));
        // without G9, G1 is 0.5 ahead of G8: best_of_multiple
        let extra: Vec<&str> = EXTRA.lines().filter(|l| !l.contains("G9")).collect();
        let b = no_readthrough(&format!("{GFF}{}\n", extra.join("\n")));
        let l = label(&b, &P, Strand::Plus, &[(910, 990)]);
        assert_eq!(sum(&b, &l), (Status::BestOfMultiple, Some(Exonic), Some("G1".into())));
        // with no margin, G1's 0.0625 lead over G9 is enough; an exact tie (G1 and G4) stays tied
        let zero = Params { tie_margin: 0.0, ..P };
        assert_eq!(label(&a, &zero, Strand::Plus, &[(910, 990)]).status, Status::BestOfMultiple);
        let with_g4 = load(GFF, &[]);
        let l = label(&with_g4, &zero, Strand::Plus, &[(910, 990)]);
        assert_eq!((l.status, tied(&with_g4, &l)), (Status::Tied, "G1,G4".to_string()));
        // a gap of exactly the margin is not a tie
        let wide = Params { tie_margin: 0.5, ..P };
        assert_eq!(label(&b, &wide, Strand::Plus, &[(910, 990)]).status, Status::BestOfMultiple);
        // under 20% of the read on exons is intronic_only: 10 of 60 bp on G1's exon
        assert_eq!(label(&a, &P, Strand::Plus, &[(440, 500)]).support, Some(IntronicOnly));
    }

    #[test]
    fn intronic() {
        let a = load(&format!("{GFF}{EXTRA}"), &[]);
        // inside G10's long intron: one gene, intronic_only; G10 has no name, so its ID stands in
        let l = label(&a, &P, Strand::Plus, &[(6000, 6100)]);
        assert_eq!(sum(&a, &l), (Status::UniqueGene, Some(IntronicOnly), Some("G10".into())));
        assert_eq!(col(&a, &l, "GeneAssignment_GeneName"), Val::Str("G10"));
        assert_eq!(col(&a, &l, "GeneAssignment_IntronsMatched"), Val::I64(0));
        // inside the bodies of G1, G3, G4, G8, G9 with no exon: all tied
        let l = label(&a, &P, Strand::Plus, &[(460, 480)]);
        assert_eq!(sum(&a, &l), (Status::Tied, Some(IntronicOnly), Some("G1".into())));
        assert_eq!(tied(&a, &l), "G1,G3,G4,G8,G9");
        assert_eq!(col(&a, &l, "GeneAssignment_NumBestGenes"), Val::I64(5));
    }

    #[test]
    fn gene_inside_a_read_intron_is_not_touched() {
        let a = load(GFF, &[]);
        // the read's intron (540,620) holds all of G5 (550,600); G2 is on the other strand
        let l = label(&a, &P, Strand::Plus, &[(520, 540), (620, 640)]);
        let cands = col(&a, &l, "GeneAssignment_Candidates");
        assert!(matches!(cands, Val::Str(s) if !s.contains("G5")), "{cands:?}");
        assert_eq!(col(&a, &l, "GeneAssignment_AntisenseGenes"), Val::Str("TWO (G2)"));
        assert_eq!(l.support, Some(IntronicOnly));
    }

    #[test]
    fn near_antisense_none_unmapped() {
        let a = load(GFF, &[]);
        // 200 bp past the 3' ends of G1, G3, G4 (a tie: the smaller gene_id); G6 is 1,700 bp away
        let l = label(&a, &P, Strand::Plus, &[(1200, 1300)]);
        assert_eq!(sum(&a, &l), (Status::NearbyGene, None, Some("G1".into())));
        assert_eq!((col(&a, &l, "GeneAssignment_NearbyGeneDist"), col(&a, &l, "GeneAssignment_NearbyGeneSide")), (Val::I64(200), Val::Str("downstream_of_3p")));
        assert_eq!((col(&a, &l, "GeneAssignment_IntronsMatched"), col(&a, &l, "GeneAssignment_NumBestGenes")), (Val::Null, Val::I64(0)));
        // 100 bp before G6's 5' end
        let l = label(&a, &P, Strand::Plus, &[(2800, 2900)]);
        assert_eq!((l.gene.map(|g| a.genes[g].id.as_str()), l.near), (Some("G6"), Some((100, Side::UpstreamOf5p))));
        // a '-' read on G1's exon: G2 is 300 bp away, so nearby_gene within 5 kb, antisense_only within
        // 100 bp
        let near100 = Params { near_max_dist: 100, ..P };
        assert_eq!(label(&a, &P, Strand::Minus, &[(150, 200)]).status, Status::NearbyGene);
        let l = label(&a, &near100, Strand::Minus, &[(150, 200)]);
        assert_eq!(sum(&a, &l), (Status::AntisenseOnly, None, None));
        assert_eq!(col(&a, &l, "GeneAssignment_AntisenseGenes"), Val::Str("ONE (G1); THREE (G3); ONE-RT (G4)"));
        assert_eq!(label(&a, &near100, Strand::Plus, &[(9100, 9200)]).status, Status::Intergenic);
        let u = GeneLabel::unmapped();
        let v = u.values(&a);
        assert_eq!(v.len(), GeneLabel::fields().len());
        assert_eq!((v[0], v[1], v[6], v[8], v[9], v[10]), (Val::Str("unmapped"), Val::Null, Val::I64(0), Val::Null, Val::I64(0), Val::I64(0)));
    }

    #[test]
    fn overlap_helpers() {
        assert_eq!(overlap_bp(&[(0, 10), (20, 30)], &[(5, 25)]), 10);
        assert_eq!(overlap_bp(&[(0, 100)], &[(10, 20), (30, 40), (90, 200)]), 30);
        assert!(touches_blocks(&[(0, 10), (50, 60)], 55, 70));
        assert!(!touches_blocks(&[(0, 10), (50, 60)], 10, 50));
    }
}
