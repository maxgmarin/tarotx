//! The gene annotation: genes from GFF3 `gene` lines, and for each gene its merged exons and its
//! annotated introns, built from the `exon` lines. All coordinates are 0-based half-open.

use std::collections::HashMap;
use std::ops::Range;

use anyhow::{bail, Context, Result};

use crate::common::io::for_each_line;
use crate::common::strand::Strand;

#[derive(Debug)]
pub struct Gene {
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    pub id: String,
    pub name: Option<String>,
    pub gene_type: Option<String>,
    exons: Range<usize>,
    introns: Range<usize>,
}

/// Genes of one chromosome: a range of `Annotation::genes` (sorted by start), the longest gene
/// length, and per strand the gene indices sorted by start and by end, for the nearest-gene search.
#[derive(Default)]
struct Chrom {
    genes: Range<usize>,
    max_len: u32,
    by_start: [Vec<usize>; 2],
    by_end: [Vec<usize>; 2],
}

pub struct Options<'a> {
    /// Keep only these gene types (all when empty).
    pub gene_types: &'a [String],
    /// Leave out genes tagged `readthrough_gene`.
    pub skip_readthrough: bool,
}

pub struct Annotation {
    pub genes: Vec<Gene>,
    chroms: HashMap<String, Chrom>,
    exon_union: Vec<(u32, u32)>,
    introns: Vec<(u32, u32)>,
    pub n_introns: usize,
}

/// Value of `key` in a GFF3 attribute column (`k1=v1;k2=v2`).
fn attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attrs.split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='))
}

struct RawGene {
    chrom: String,
    gene: Gene,
}

/// An exon line, with its gene and transcript as interned keys.
struct RawExon {
    gene: u32,
    tx: u32,
    start: u32,
    end: u32,
}

type Intervals = Vec<(u32, u32)>;

fn intern(map: &mut HashMap<String, u32>, s: &str) -> u32 {
    if let Some(&k) = map.get(s) {
        return k;
    }
    let k = map.len() as u32;
    map.insert(s.to_string(), k);
    k
}

impl Annotation {
    /// Read the `gene` and `exon` lines of a GFF3 (gzipped or not). Exons name their gene by the
    /// `gene_id` attribute and their transcript by `Parent`, as in GENCODE.
    pub fn load(path: &str, opts: &Options) -> Result<Annotation> {
        let mut gene_keys: HashMap<String, u32> = HashMap::new();
        let mut tx_keys: HashMap<String, u32> = HashMap::new();
        let mut genes: Vec<Option<RawGene>> = Vec::new(); // by gene key
        let mut exons: Vec<RawExon> = Vec::new();
        let mut in_fasta = false;
        for_each_line(path, |ln, line| {
            if in_fasta || line.starts_with('#') || line.is_empty() {
                in_fasta |= line.starts_with("##FASTA");
                return Ok(());
            }
            let f: Vec<&str> = line.splitn(9, '\t').collect();
            if f.len() < 9 {
                bail!("{path}:{ln}: expected 9 GFF3 columns, found {}", f.len());
            }
            if f[2] != "gene" && f[2] != "exon" {
                return Ok(());
            }
            let num = |s: &str, what: &str| s.parse::<u32>().with_context(|| format!("{path}:{ln}: bad {what} '{s}'"));
            let (start, end) = (num(f[3], "start")? - 1, num(f[4], "end")?);
            let attrs = f[8];
            let Some(gene_id) = attr(attrs, "gene_id") else {
                bail!("{path}:{ln}: {} line has no gene_id attribute (a GENCODE-style GFF3 is needed)", f[2]);
            };
            let key = intern(&mut gene_keys, gene_id) as usize;
            if genes.len() <= key {
                genes.resize_with(key + 1, || None);
            }
            if f[2] == "exon" {
                let Some(parent) = attr(attrs, "Parent") else { bail!("{path}:{ln}: exon line has no Parent attribute") };
                let tx = intern(&mut tx_keys, parent);
                exons.push(RawExon { gene: key as u32, tx, start, end });
                return Ok(());
            }
            let Some(strand) = Strand::parse(f[6]) else { return Ok(()) }; // unstranded genes can't be assigned
            if genes[key].is_some() {
                bail!("{path}:{ln}: gene_id '{gene_id}' appears on more than one gene line");
            }
            let gene_type = attr(attrs, "gene_type");
            if !opts.gene_types.is_empty() && !gene_type.is_some_and(|t| opts.gene_types.iter().any(|k| k == t)) {
                return Ok(());
            }
            if opts.skip_readthrough && attr(attrs, "tag").is_some_and(|t| t.split(',').any(|t| t == "readthrough_gene")) {
                return Ok(());
            }
            genes[key] = Some(RawGene {
                chrom: f[0].to_string(),
                gene: Gene {
                    start,
                    end,
                    strand,
                    id: gene_id.to_string(),
                    name: attr(attrs, "gene_name").map(str::to_string),
                    gene_type: gene_type.map(str::to_string),
                    exons: 0..0,
                    introns: 0..0,
                },
            });
            Ok(())
        })?;

        // Exon union per gene: sort by gene and start; a new interval starts past the running end.
        exons.retain(|e| genes.get(e.gene as usize).is_some_and(Option::is_some));
        exons.sort_unstable_by_key(|e| (e.gene, e.start, e.end));
        let mut union: Vec<Vec<(u32, u32)>> = (0..genes.len()).map(|_| Vec::new()).collect();
        for e in &exons {
            let u = &mut union[e.gene as usize];
            match u.last_mut() {
                Some(last) if e.start <= last.1 => last.1 = last.1.max(e.end),
                _ => u.push((e.start, e.end)),
            }
        }
        // Annotated introns: gaps between consecutive exons of each transcript, unique per gene.
        exons.sort_unstable_by_key(|e| (e.tx, e.start, e.end));
        let mut introns: Vec<Vec<(u32, u32)>> = (0..genes.len()).map(|_| Vec::new()).collect();
        for w in exons.windows(2) {
            if w[0].tx == w[1].tx && w[1].start > w[0].end {
                introns[w[1].gene as usize].push((w[0].end, w[1].start));
            }
        }
        drop(exons);

        // Genes sorted by chromosome (in first-seen order) and start; per-gene lists in flat arrays.
        let mut chrom_order: HashMap<String, usize> = HashMap::new();
        for g in genes.iter().flatten() {
            let n = chrom_order.len();
            chrom_order.entry(g.chrom.clone()).or_insert(n);
        }
        let mut raw: Vec<(RawGene, Intervals, Intervals)> = genes
            .into_iter()
            .zip(union)
            .zip(introns)
            .filter_map(|((g, u), i)| g.map(|g| (g, u, i)))
            .collect();
        raw.sort_by(|a, b| {
            (chrom_order[&a.0.chrom], a.0.gene.start, a.0.gene.end, &a.0.gene.id)
                .cmp(&(chrom_order[&b.0.chrom], b.0.gene.start, b.0.gene.end, &b.0.gene.id))
        });
        let mut ann =
            Annotation { genes: Vec::new(), chroms: HashMap::new(), exon_union: Vec::new(), introns: Vec::new(), n_introns: 0 };
        for (i, (rg, u, mut intr)) in raw.into_iter().enumerate() {
            let mut g = rg.gene;
            intr.sort_unstable();
            intr.dedup();
            g.exons = ann.exon_union.len()..ann.exon_union.len() + u.len();
            ann.exon_union.extend(u);
            g.introns = ann.introns.len()..ann.introns.len() + intr.len();
            ann.introns.extend(intr);
            let c = ann.chroms.entry(rg.chrom).or_default();
            if c.genes.is_empty() {
                c.genes = i..i;
            }
            c.genes.end = i + 1;
            c.max_len = c.max_len.max(g.end - g.start);
            ann.genes.push(g);
        }
        ann.n_introns = ann.introns.len();
        for c in ann.chroms.values_mut() {
            for s in 0..2 {
                let mut idx: Vec<usize> = c.genes.clone().filter(|&i| ann.genes[i].strand as usize == s).collect();
                c.by_start[s] = idx.clone(); // genes are already sorted by start
                idx.sort_by_key(|&i| ann.genes[i].end);
                c.by_end[s] = idx;
            }
        }
        Ok(ann)
    }

    pub fn has_chrom(&self, chrom: &str) -> bool {
        self.chroms.contains_key(chrom)
    }

    /// Indices of the genes (either strand) whose body overlaps [start, end), in start order.
    pub fn overlapping(&self, chrom: &str, start: u32, end: u32) -> Vec<usize> {
        let Some(c) = self.chroms.get(chrom) else { return Vec::new() };
        let genes = &self.genes[c.genes.clone()];
        // a gene can only overlap if it starts in [start - max_len, end)
        let lo = genes.partition_point(|g| g.start < start.saturating_sub(c.max_len));
        let hi = genes.partition_point(|g| g.start < end);
        (lo..hi).filter(|&i| genes[i].end > start).map(|i| c.genes.start + i).collect()
    }

    /// Same-strand genes bordering [start, end) on a chromosome: the gene(s) with the largest end
    /// at or before `start`, and the gene(s) with the smallest start at or after `end`.
    pub fn flanking(&self, chrom: &str, strand: Strand, start: u32, end: u32) -> Vec<usize> {
        let Some(c) = self.chroms.get(chrom) else { return Vec::new() };
        let (by_end, by_start) = (&c.by_end[strand as usize], &c.by_start[strand as usize]);
        let mut out = Vec::new();
        let left = by_end.partition_point(|&i| self.genes[i].end <= start);
        if left > 0 {
            let best = self.genes[by_end[left - 1]].end;
            out.extend(by_end[..left].iter().rev().take_while(|&&i| self.genes[i].end == best));
        }
        let right = by_start.partition_point(|&i| self.genes[i].start < end);
        if let Some(&first) = by_start.get(right) {
            let best = self.genes[first].start;
            out.extend(by_start[right..].iter().take_while(|&&i| self.genes[i].start == best));
        }
        out
    }

    /// The gene's exons, merged across its transcripts, sorted.
    pub fn exons(&self, gene: usize) -> &[(u32, u32)] {
        &self.exon_union[self.genes[gene].exons.clone()]
    }

    /// The gene's annotated introns, unique and sorted.
    pub fn introns(&self, gene: usize) -> &[(u32, u32)] {
        &self.introns[self.genes[gene].introns.clone()]
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    /// A small annotation on chrT, used here and by the labelling tests.
    pub const GFF: &str = "##gff-version 3
chrT\tt\tgene\t101\t1000\t.\t+\t.\tID=G1;gene_id=G1;gene_name=ONE;gene_type=protein_coding
chrT\tt\texon\t101\t200\t.\t+\t.\tParent=T1a;gene_id=G1
chrT\tt\texon\t301\t400\t.\t+\t.\tParent=T1a;gene_id=G1
chrT\tt\texon\t901\t1000\t.\t+\t.\tParent=T1a;gene_id=G1
chrT\tt\texon\t151\t250\t.\t+\t.\tParent=T1b;gene_id=G1
chrT\tt\texon\t401\t450\t.\t+\t.\tParent=T1b;gene_id=G1
chrT\tt\texon\t901\t1000\t.\t+\t.\tParent=T1b;gene_id=G1
chrT\tt\tgene\t501\t700\t.\t-\t.\tID=G2;gene_id=G2;gene_name=TWO;gene_type=lncRNA
chrT\tt\texon\t501\t700\t.\t-\t.\tParent=T2;gene_id=G2
chrT\tt\tgene\t101\t1000\t.\t+\t.\tID=G3;gene_id=G3;gene_name=THREE;gene_type=processed_pseudogene
chrT\tt\texon\t101\t200\t.\t+\t.\tParent=T3;gene_id=G3
chrT\tt\texon\t301\t400\t.\t+\t.\tParent=T3;gene_id=G3
chrT\tt\tgene\t101\t1000\t.\t+\t.\tID=G4;gene_id=G4;gene_name=ONE-RT;gene_type=protein_coding;tag=readthrough_gene
chrT\tt\texon\t101\t200\t.\t+\t.\tParent=T4;gene_id=G4
chrT\tt\texon\t301\t400\t.\t+\t.\tParent=T4;gene_id=G4
chrT\tt\texon\t901\t1000\t.\t+\t.\tParent=T4;gene_id=G4
chrT\tt\tgene\t551\t600\t.\t+\t.\tID=G5;gene_id=G5;gene_name=SMALL;gene_type=snRNA
chrT\tt\texon\t551\t600\t.\t+\t.\tParent=T5;gene_id=G5
chrT\tt\tgene\t3001\t4000\t.\t+\t.\tID=G6;gene_id=G6;gene_name=FAR;gene_type=protein_coding
chrT\tt\texon\t3001\t4000\t.\t+\t.\tParent=T6;gene_id=G6
";

    pub fn load(text: &str, gene_types: &[String]) -> Annotation {
        load_opts(text, &Options { gene_types, skip_readthrough: false })
    }

    pub fn load_opts(text: &str, opts: &Options) -> Annotation {
        let p = std::env::temp_dir().join(format!("tarotx-annot-{}-{:?}.gff3", std::process::id(), std::thread::current().id()));
        std::fs::write(&p, text).unwrap();
        let a = Annotation::load(p.to_str().unwrap(), opts);
        std::fs::remove_file(p).unwrap();
        a.unwrap()
    }

    fn gene(a: &Annotation, id: &str) -> usize {
        a.genes.iter().position(|g| g.id == id).unwrap()
    }

    #[test]
    fn exon_union_and_introns() {
        let a = load(GFF, &[]);
        let g1 = gene(&a, "G1");
        // T1a [100,200) [300,400) [900,1000) and T1b [150,250) [400,450) [900,1000)
        assert_eq!(a.exons(g1), [(100, 250), (300, 450), (900, 1000)]);
        assert_eq!(a.introns(g1), [(200, 300), (250, 400), (400, 900), (450, 900)]);
        assert_eq!(a.introns(gene(&a, "G2")), []);
        assert_eq!(a.n_introns, 4 + 1 + 2);
    }

    #[test]
    fn gene_flags() {
        let a = load(GFF, &[]);
        let no_rt = load_opts(GFF, &Options { gene_types: &[], skip_readthrough: true });
        assert!(no_rt.genes.iter().all(|g| g.id != "G4") && no_rt.genes.len() == a.genes.len() - 1);
        assert_eq!(a.genes[gene(&a, "G2")].strand, Strand::Minus);
        assert_eq!(a.genes[gene(&a, "G1")].name.as_deref(), Some("ONE"));
        let only = load(GFF, &["lncRNA".to_string()]);
        assert_eq!(only.genes.len(), 1);
        assert_eq!(only.exons(0), [(500, 700)]);
    }

    #[test]
    fn overlap_and_flanking_queries() {
        let a = load(GFF, &[]);
        let ids = |v: Vec<usize>| v.into_iter().map(|i| a.genes[i].id.as_str()).collect::<Vec<_>>();
        assert_eq!(ids(a.overlapping("chrT", 560, 570)), ["G1", "G3", "G4", "G2", "G5"]);
        assert_eq!(ids(a.overlapping("chrT", 1000, 1100)), Vec::<&str>::new());
        assert_eq!(ids(a.overlapping("chrX", 0, 100)), Vec::<&str>::new());
        // '+' genes ending at or before 2000, and starting at or after 2100
        let mut plus = ids(a.flanking("chrT", Strand::Plus, 2000, 2100));
        plus.sort();
        assert_eq!(plus, ["G1", "G3", "G4", "G6"]);
        assert_eq!(ids(a.flanking("chrT", Strand::Minus, 2000, 2100)), ["G2"]);
    }

    #[test]
    fn gff_errors() {
        let p = std::env::temp_dir().join(format!("tarotx-bad-{}.gff3", std::process::id()));
        let opts = Options { gene_types: &[], skip_readthrough: false };
        std::fs::write(&p, "chrT\tt\texon\t1\t10\t.\t+\t.\tParent=T1\n").unwrap();
        let e = Annotation::load(p.to_str().unwrap(), &opts).err().unwrap();
        assert!(e.to_string().contains("no gene_id"), "{e}");
        std::fs::write(&p, "chrT\tt\tgene\t1\t10\t.\t+\t.\tgene_id=G\nchrT\tt\tgene\t1\t10\t.\t+\t.\tgene_id=G\n").unwrap();
        let e = Annotation::load(p.to_str().unwrap(), &opts).err().unwrap();
        assert!(e.to_string().contains("more than one gene line"), "{e}");
        std::fs::remove_file(p).unwrap();
    }
}
