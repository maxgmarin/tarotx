//! End-to-end tests of `tarotx label-by-gene` on a hand-built chromosome (`tests/data/label_by_gene`). `anno.gff3` is the annotation of the
//! unit tests (genes G1-G6 on chrT). The same eight reads are in `aln.paf` (cs tags), `aln_cg.paf` (cg
//! tags), `aln.tsv` and `aln.parquet` (snappy; both tables have `cs_A` and `cigar_A` columns); the
//! PAFs also have a supplementary alignment of r1 on their last line.
//!
//! G1, G3 (pseudogene) and G4 (readthrough) share most exons and introns, so many reads tie.
//!
//! r1 [150,200)+[300,400)  intron (200,300) in G1, G3, G4       tied         known_junctions
//! r2 [910,990)            on the last exon of G1 and G4        tied         exonic  (G3: no exon here)
//! r3 [460,480)            in the bodies of G1, G3, G4, no exon tied         intronic_only
//! r4 [1200,1300)          200 bp past G1's 3' end              nearby_gene  G1
//! r5 [150,200) '-'        on G1 (antisense); G2 300 bp away    nearby_gene  G2
//! r6 [9100,9200)          5,100 bp from G6                     intergenic
//! r7 unmapped                                                  unmapped
//! r8 [440,500)            10 of 60 bp on G1's exons (< 20%)    tied         intronic_only

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

fn data(f: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/label_by_gene").join(f).to_string_lossy().into_owned()
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("tarotx-label-by-gene-test-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d.join(name)
}

/// Run `tarotx label-by-gene` with `args`, returning (stdout, stderr). Panics if it fails.
fn label_by_gene(args: &[&str]) -> (String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_tarotx")).arg("label-by-gene").args(args).output().unwrap();
    let (o, e) = (String::from_utf8(out.stdout).unwrap(), String::from_utf8(out.stderr).unwrap());
    assert!(out.status.success(), "tarotx label-by-gene {args:?} failed:\n{e}");
    (o, e)
}

/// Run `tarotx label-by-gene` expecting it to fail; returns stderr.
fn label_by_gene_err(args: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_tarotx")).arg("label-by-gene").args(args).output().unwrap();
    assert!(!out.status.success(), "tarotx label-by-gene {args:?} should have failed");
    String::from_utf8(out.stderr).unwrap()
}

type Rows = Vec<HashMap<String, String>>;

/// TSV text as rows of column -> value.
fn rows(tsv: &str) -> Rows {
    let mut lines = tsv.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines.map(|l| header.iter().map(|h| h.to_string()).zip(l.split('\t').map(String::from)).collect()).collect()
}

fn col(rows: &Rows, c: &str) -> Vec<String> {
    rows.iter().map(|r| r[c].clone()).collect()
}

/// Label `aln` to TSV on stdout, with extra options.
fn run(aln: &str, extra: &[&str]) -> Rows {
    let (a, g) = (data(aln), data("anno.gff3"));
    let mut args = vec!["-a", a.as_str(), "-g", g.as_str(), "-o", "-"];
    args.extend(extra);
    rows(&label_by_gene(&args).0)
}

const STATUS: [&str; 8] = ["tied", "tied", "tied", "nearby_gene", "nearby_gene", "intergenic", "unmapped", "tied"];

#[test]
fn paf_labels() {
    let r = run("aln.paf", &[]);
    // every line, unmapped and supplementary included
    assert_eq!(col(&r, "query_name"), ["r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r1"]);
    let r: Rows = r[..8].to_vec();
    assert_eq!(col(&r, "GeneAssignment_Status"), STATUS);
    assert_eq!(col(&r, "GeneAssignment_SupportLevel"), ["known_junctions", "exonic", "intronic_only", "", "", "", "", "intronic_only"]);
    assert_eq!(col(&r, "GeneAssignment_GeneID"), ["G1", "G1", "G1", "G1", "G2", "", "", "G1"]);
    assert_eq!(col(&r, "GeneAssignment_GeneName"), ["ONE", "ONE", "ONE", "ONE", "TWO", "", "", "ONE"]);
    assert_eq!(col(&r, "GeneAssignment_IntronsMatched"), ["1", "0", "0", "", "", "", "", "0"]);
    assert_eq!(col(&r, "GeneAssignment_ReadIntrons"), ["1", "0", "0", "0", "0", "0", "0", "0"]);
    assert_eq!(col(&r, "GeneAssignment_ExonicFrac"), ["1.0", "1.0", "0.0", "", "", "", "", "0.167"]);
    assert_eq!(col(&r, "GeneAssignment_AlignedBp"), ["150", "80", "20", "100", "50", "100", "", "60"]);
    assert_eq!(col(&r, "GeneAssignment_NumSameStrandGenes"), ["3", "3", "3", "0", "0", "0", "0", "3"]);
    assert_eq!(col(&r, "GeneAssignment_NumBestGenes"), ["3", "2", "3", "0", "0", "0", "0", "3"]);
    assert_eq!(col(&r, "tp_tag"), ["P"; 8]);
    assert_eq!(r[1]["GeneAssignment_Candidates"], "ONE (G1) 0/0 1.0; ONE-RT (G4) 0/0 1.0; THREE (G3) 0/0 0.0");
    assert_eq!(r[4]["GeneAssignment_AntisenseGenes"], "ONE (G1); THREE (G3); ONE-RT (G4)");
    assert_eq!(col(&r, "GeneAssignment_NearbyGeneDist"), ["", "", "", "200", "300", "", "", ""]);
    assert_eq!(col(&r, "GeneAssignment_NearbyGeneSide")[3..5], ["downstream_of_3p", "downstream_of_3p"]);
    assert_eq!(r[6]["target_name"], "");
}

#[test]
fn every_input_format_agrees() {
    let gene_cols = |r: &Rows| -> Vec<Vec<(String, String)>> {
        r.iter()
            .map(|row| {
                let mut v: Vec<(String, String)> =
                    row.iter().filter(|(k, _)| k.starts_with("GeneAssignment_")).map(|(k, v)| (k.clone(), v.clone())).collect();
                v.sort();
                v
            })
            .collect()
    };
    let paf = run("aln.paf", &[]);
    assert_eq!(paf.len(), 9);
    assert_eq!(gene_cols(&run("aln_cg.paf", &[])), gene_cols(&paf));
    // the tables hold the first eight reads
    assert_eq!(gene_cols(&run("aln.tsv", &[])), gene_cols(&paf)[..8]);
    assert_eq!(gene_cols(&run("aln.parquet", &[])), gene_cols(&paf)[..8]);
    assert_eq!(gene_cols(&run("aln.tsv", &["--cigar-col", "cigar_A"])), gene_cols(&paf)[..8]);
    assert_eq!(gene_cols(&run("aln.parquet", &["--cigar-col", "cigar_A"])), gene_cols(&paf)[..8]);
}

#[test]
fn one_standard_table_for_every_input() {
    // the same leading columns whatever the input, and none of a table's other columns
    let lead = ["query_name", "tp_tag", "target_name", "target_start", "target_end", "strand", "GeneAssignment_Status"];
    let g = data("anno.gff3");
    for out in [
        label_by_gene(&["-a", &data("aln.paf"), "-g", &g, "-o", "-"]).0,
        label_by_gene(&["-a", &data("aln.tsv"), "-g", &g, "-o", "-"]).0,
        label_by_gene(&["-a", &data("aln.parquet"), "-g", &g, "-o", "-"]).0,
    ] {
        let header: Vec<&str> = out.lines().next().unwrap().split('\t').collect();
        assert_eq!(header.len(), 6 + 15);
        assert_eq!(header[..7], lead);
    }
    // a table keeps every row in order, unmapped included; it has no tp column here, so tp_tag is empty.
    // Parquet in, Parquet out gives the same rows.
    let pq = tmp("labels.parquet");
    label_by_gene(&["-a", &data("aln.parquet"), "-g", &g, "-o", pq.to_str().unwrap()]);
    assert_eq!(std::fs::read(&pq).unwrap()[..4], *b"PAR1");
    let r = run("aln.parquet", &[]);
    assert_eq!(r, run("aln.tsv", &[]));
    assert_eq!(col(&r, "query_name"), ["r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8"]);
    assert_eq!(col(&r, "GeneAssignment_Status"), STATUS);
    assert_eq!((r[0]["tp_tag"].as_str(), r[0]["target_start"].as_str(), r[0]["strand"].as_str()), ("", "150", "+"));
    assert_eq!((r[6]["target_name"].as_str(), r[6]["target_start"].as_str()), ("", ""));
}

#[test]
fn table_intron_columns() {
    // a table with its own column names; the introns are walked from its cs or CIGAR column
    let t = tmp("cs.tsv");
    let head = "id\tchr\tstr\tbeg\tfin\tcs\tcig\n";
    std::fs::write(&t, format!("{head}r1\tchrT\t+\t150\t400\t:50~gt100ag:100\t50M100N100M\nrx\tchrT\t+\t150\t401\t:50~gt100ag:100\t50M100N100M\n"))
        .unwrap();
    let (tp, g) = (t.to_str().unwrap().to_string(), data("anno.gff3"));
    let cols = ["--id-col", "id", "--chrom-col", "chr", "--strand-col", "str", "--start-col", "beg", "--end-col", "fin"];
    let args = |extra: &[&'static str]| -> Vec<String> {
        ["-a", &tp, "-g", &g, "-o", "-"].iter().chain(&cols).chain(extra).map(|s| s.to_string()).collect()
    };
    let fails = |extra: &[&'static str]| label_by_gene_err(&args(extra).iter().map(String::as_str).collect::<Vec<_>>());
    // the second row's cs and CIGAR cover 250 bp but its span is 251: an error naming the row
    let err = fails(&["--cs-col", "cs"]);
    assert!(err.contains("cs.tsv:3") && err.contains("the cs tag covers 250 bp"), "{err}");
    let err = fails(&["--cigar-col", "cig"]);
    assert!(err.contains("cs.tsv:3") && err.contains("the CIGAR covers 250 bp"), "{err}");
    std::fs::write(&t, format!("{head}r1\tchrT\t+\t150\t400\t:50~gt100ag:100\t50M100N100M\n")).unwrap();
    for extra in [["--cs-col", "cs"], ["--cigar-col", "cig"]] {
        let r = rows(&label_by_gene(&args(&extra).iter().map(String::as_str).collect::<Vec<_>>()).0);
        assert_eq!((r[0]["GeneAssignment_Status"].as_str(), r[0]["GeneAssignment_SupportLevel"].as_str()), ("tied", "known_junctions"));
    }
    // the cs column must be there (here the default cs_A isn't); the two columns can't both be set
    assert!(fails(&[]).contains("no column [\"cs_A\"]"));
    assert!(fails(&["--cs-col", "cs", "--cigar-col", "cig"]).contains("cannot be used with"));
    // precomputed junction strings are not read
    std::fs::write(&t, "query_name\ttarget_name_A\tstrand_A\ttarget_start_A\ttarget_end_A\tgenomic_junctions_A\nr1\tchrT\t+\t150\t400\t((200, 300),)\n")
        .unwrap();
    let err = label_by_gene_err(&["-a", &tp, "-g", &g, "-o", "-"]);
    assert!(err.contains("no column [\"cs_A\"]") && err.contains("--cigar-col"), "{err}");
}

#[test]
fn options_and_errors() {
    let (paf, g) = (data("aln.paf"), data("anno.gff3"));
    // a smaller near distance: r4 is intergenic, and r5 (on G1, opposite strand) antisense_only
    let r = run("aln.paf", &["--near-max-dist", "100"]);
    assert_eq!(col(&r, "GeneAssignment_Status")[3..5], ["intergenic", "antisense_only"]);
    // without the readthrough gene G4, r2 has one exonic gene; r1 still ties G1 with G3
    let r = run("aln.paf", &["--skip-readthrough"]);
    assert_eq!((r[1]["GeneAssignment_Status"].as_str(), r[1]["GeneAssignment_SupportLevel"].as_str()), ("unique_gene", "exonic"));
    assert_eq!(r[0]["GeneAssignment_NumBestGenes"], "2"); // G1 and G3
    assert!(r[0]["GeneAssignment_Candidates"].starts_with("ONE (G1) 1/1 1.0; THREE (G3) 1/1 1.0"));
    // protein-coding genes only, and no readthrough: r1 has one gene left
    let r = run("aln.paf", &["--skip-readthrough", "--gene-types", "protein_coding"]);
    assert_eq!(r[0]["GeneAssignment_Status"], "unique_gene");
    // with no margin only exact ties are tied, and r2 fits G1 and G4 exactly equally
    let r = run("aln.paf", &["--tie-margin", "0"]);
    assert_eq!(r[1]["GeneAssignment_Status"], "tied");
    // a lower exonic minimum makes r8 exonic (10 of 60 bp on G1, more than G3 and G4)
    let r = run("aln.paf", &["--min-exonic-frac", "0.1"]);
    assert_eq!((r[7]["GeneAssignment_Status"].as_str(), r[7]["GeneAssignment_SupportLevel"].as_str()), ("unique_gene", "exonic"));
    assert!(label_by_gene_err(&["-a", &paf, "-g", &g, "-o", "-", "--labels-only"]).contains("unexpected argument"));
    assert!(label_by_gene_err(&["-a", &paf, "-g", &g, "-o", "-", "--tagged-out", "x.paf"]).contains("unexpected argument"));
    assert!(label_by_gene_err(&["-a", "x.bam", "-g", &g, "-o", "-"]).contains("convert it to PAF first"));
    let bad = tmp("nocs.paf");
    std::fs::write(&bad, "r\t10\t0\t10\t+\tchrT\t10000\t0\t10\t10\t10\t60\ttp:A:P\n").unwrap();
    assert!(label_by_gene_err(&["-a", bad.to_str().unwrap(), "-g", &g, "-o", "-"]).contains("no cs:Z or cg:Z tag"));
}

#[test]
fn summary_file() {
    let s = tmp("summary.tsv");
    let (paf, g) = (data("aln.paf"), data("anno.gff3"));
    label_by_gene(&["-a", &paf, "-g", &g, "-o", "-", "--summary", s.to_str().unwrap()]);
    let text = std::fs::read_to_string(&s).unwrap();
    let rows = rows(&text);
    assert_eq!(rows.len(), 12); // every combination that can occur, zeros included
    let get = |st: &str, lvl: &str| {
        let r = rows.iter().find(|r| r["GeneAssignment_Status"] == st && r["GeneAssignment_SupportLevel"] == lvl).unwrap();
        (r["count"].clone(), r["percent"].clone())
    };
    // 9 PAF lines: r1, r2 and the r1 supplementary (on r2's exon) tie with known junctions or exonic
    // support, r3 and r8 tie intronic-only, r4 and r5 are nearby_gene, r6 intergenic, r7 unmapped
    assert_eq!(get("tied", "known_junctions"), ("1".into(), "11.11".into()));
    assert_eq!(get("tied", "exonic"), ("2".into(), "22.22".into()));
    assert_eq!(get("tied", "intronic_only"), ("2".into(), "22.22".into()));
    assert_eq!(get("nearby_gene", "").0, "2");
    assert_eq!(get("unique_gene", "exonic").0, "0");
    let total: u64 = rows.iter().map(|r| r["count"].parse::<u64>().unwrap()).sum();
    assert_eq!(total, 9);
    assert!(label_by_gene_err(&["-a", &paf, "-g", &g, "-o", "x.tsv", "--summary", "x.tsv"]).contains("same file"));
}
