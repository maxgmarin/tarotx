//! End-to-end tests of `tarotx label-internal-priming` on a hand-built chromosome
//! (`tests/data/label_internal_priming`). Expected values are worked out by hand:
//!
//! tracts (0-based): + [100,110) score 90, + [200,230) score 100, - [300,310) score 80
//! annotated 3' ends: TX1 '+' at 150 (basic); TX2 '-' at 400 (basic, mRNA_end_NF); TX3 '-' at 280
//!
//! `aln.paf` holds nine alignments: r1-r7, r8 (unmapped), and a supplementary alignment of r2
//! (also tp:A:P, as minimap2 writes it). `aln.tsv` is a table with its own column names.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

fn data(f: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/label_internal_priming").join(f).to_string_lossy().into_owned()
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("tarotx-label-internal-priming-test-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d.join(name)
}

/// `tarotx label-internal-priming`, ready for its arguments.
fn subcommand() -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_tarotx"));
    c.arg("label-internal-priming");
    c
}

/// Run the subcommand, returning (stdout, stderr). Panics if it fails.
fn label_internal_priming(args: &[&str]) -> (String, String) {
    let out = subcommand().args(args).output().unwrap();
    let (o, e) = (String::from_utf8(out.stdout).unwrap(), String::from_utf8(out.stderr).unwrap());
    assert!(out.status.success(), "tarotx label-internal-priming {args:?} failed:\n{e}");
    (o, e)
}

/// TSV text as rows of column -> value.
fn rows(tsv: &str) -> Vec<HashMap<String, String>> {
    let mut lines = tsv.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| header.iter().map(|h| h.to_string()).zip(l.split('\t').map(String::from)).collect())
        .collect()
}

fn col(rows: &[HashMap<String, String>], c: &str) -> Vec<String> {
    rows.iter().map(|r| r[c].clone()).collect()
}

fn base(aln: &str) -> Vec<String> {
    ["-a", &data(aln), "-p", &data("tracts.bed"), "-g", &data("anno.gff3"), "-o", "-"].map(String::from).to_vec()
}

fn run(args: &[String]) -> Vec<HashMap<String, String>> {
    rows(&label_internal_priming(&args.iter().map(String::as_str).collect::<Vec<_>>()).0)
}

#[test]
fn paf_rows_as_given() {
    let r = run(&base("aln.paf"));
    // one row per line, in input order: r8 is unmapped (empty coordinates and labels), and the last
    // line is a supplementary piece of r2, written tp:A:P by minimap2 like a primary
    assert_eq!(col(&r, "query_name"), ["r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r2"]);
    assert_eq!(col(&r, "tp_tag"), ["P", "P", "P", "P", "P", "P", "P", "", "P"]);
    assert_eq!(col(&r, "aln_type"), ["primary", "primary", "primary", "primary", "primary", "primary", "primary", "unmapped", "primary"]);
    assert_eq!(col(&r, "InternalPriming_ReadEnd3p"), ["100", "98", "215", "233", "310", "600", "10", "", "198"]);
    assert_eq!(col(&r, "InternalPriming_PolyAEdgeDist"), ["0", "2", "15", "3", "0", "290", "", "", "2"]);
    // r3 ends 15 bp from both edges of [200,230): the edge beyond its 3' end wins, so the tract is that one
    assert_eq!(col(&r, "InternalPriming_PolyATractLen"), ["10", "10", "30", "30", "10", "10", "", "", "30"]);
    assert_eq!(col(&r, "InternalPriming_PolyAScore"), ["90", "90", "100", "100", "80", "80", "", "", "100"]);
    assert_eq!(col(&r, "InternalPriming_TxEndSignedDist"), ["-50", "-52", "65", "83", "-30", "-320", "", "", "48"]);
    // TX2 is mRNA_end_NF, so dropped: the '-' reads see TX3
    assert_eq!(col(&r, "InternalPriming_TxEndTranscriptID"), ["TX1", "TX1", "TX1", "TX1", "TX3", "TX3", "", "", "TX1"]);
    assert_eq!(col(&r, "InternalPriming_PolyAWithin3bp"), ["True", "True", "False", "True", "True", "False", "False", "", "True"]);
    assert_eq!(col(&r, "InternalPriming_TxEndWithin50bp"), ["True", "False", "False", "False", "True", "False", "False", "", "True"]);
    assert_eq!(col(&r, "InternalPriming_Flag"), ["False", "True", "False", "True", "False", "False", "False", "", "False"]);
    // query coordinates on the read as sequenced; r5 is '-' with 3 bp past its 3' end
    assert_eq!(col(&r, "query_len"), ["55", "88", "65", "83", "100", "100", "10", "10", "88"]);
    assert_eq!(col(&r, "query_start"), ["5", "0", "0", "0", "7", "0", "0", "", "10"]);
    assert_eq!(col(&r, "query_end"), ["55", "78", "65", "83", "97", "100", "10", "", "88"]);
    // r2's first piece stops 10 bases short of the read's end; its supplementary piece holds the 3' end
    assert_eq!(col(&r, "read_3p_clip"), ["0", "10", "0", "0", "3", "0", "0", "", "0"]);
    assert_eq!(col(&r, "target_name"), ["chrT", "chrT", "chrT", "chrT", "chrT", "chrT", "chrX", "", "chrT"]);
    assert_eq!(col(&r, "target_end"), ["100", "98", "215", "233", "400", "700", "10", "", "198"]);
}

#[test]
fn same_leading_columns_for_every_input() {
    // the same 12 leading columns as label-by-gene, from a PAF and from a table (whose other columns
    // are not copied)
    let lead = ["query_name", "target_name", "target_start", "target_end", "strand", "query_len", "query_start", "query_end", "mapq", "tp_tag", "aln_type", "read_3p_clip"];
    let header = |args: &[String]| -> Vec<String> {
        let (out, _) = label_internal_priming(&args.iter().map(String::as_str).collect::<Vec<_>>());
        out.lines().next().unwrap().split('\t').map(String::from).collect()
    };
    let mut a = base("aln.tsv");
    a.extend(TSV_COLS.map(String::from));
    for h in [header(&base("aln.paf")), header(&a)] {
        assert_eq!(h[..12], lead);
        assert_eq!(h[12], "InternalPriming_ReadEnd3p");
    }
    // the fixture table has no read length, query span or mapq: those columns are empty
    let r = run(&a);
    assert_eq!((r[0]["query_len"].as_str(), r[0]["read_3p_clip"].as_str(), r[0]["aln_type"].as_str()), ("", "", "primary"));
}

#[test]
fn table_record_columns() {
    // a table with the read length, query span and mapq (maligno's names by default) gets the
    // record columns too, with the same values as the PAF
    let t = tmp("record.tsv");
    std::fs::write(
        &t,
        "query_name\tquery_len\ttarget_name_A\tstrand_A\ttarget_start_A\ttarget_end_A\tquery_start_A\tquery_end_A\tmapq_A\ttp_tag_A\n\
         r2\t88\tchrT\t+\t20\t98\t0\t78\t60\tP\n\
         r5\t100\tchrT\t-\t310\t400\t7\t97\t60\tP\n\
         r8\t10\t*\t*\t\t\t\t\t\t\n",
    )
    .unwrap();
    let mut a = base("aln.tsv");
    a[1] = t.to_string_lossy().into_owned();
    let tab = run(&a);
    let paf = run(&base("aln.paf"));
    for (row, i) in tab.iter().zip([1, 4, 7]) {
        for c in ["query_len", "query_start", "query_end", "mapq", "aln_type", "read_3p_clip", "InternalPriming_ReadEnd3p", "InternalPriming_Flag"] {
            assert_eq!(row[c], paf[i][c], "column {c}, read {}", row["query_name"]);
        }
    }
    assert_eq!(col(&tab, "read_3p_clip"), ["10", "3", ""]);
    assert_eq!(col(&tab, "aln_type"), ["primary", "primary", "unmapped"]);

    // with only some of the four columns, the others are empty, and a note says which are missing
    let mut b = a.clone();
    b.extend(["--mapq-col", "no_such_column"].map(String::from));
    let (out, err) = label_internal_priming(&b.iter().map(String::as_str).collect::<Vec<_>>());
    assert!(err.contains("no column [\"no_such_column\"]") && err.contains("--mapq-col"), "{err}");
    let r = rows(&out);
    assert_eq!((r[0]["mapq"].as_str(), r[0]["query_len"].as_str(), r[0]["read_3p_clip"].as_str()), ("", "88", "10"));
}

/// The column options for `aln.tsv`.
const TSV_COLS: [&str; 12] =
    ["--id-col", "id", "--tp-col", "tp", "--chrom-col", "chr", "--strand-col", "str", "--start-col", "beg", "--end-col", "fin"];

#[test]
fn tsv_input() {
    let mut a = base("aln.tsv");
    a.extend(TSV_COLS.map(String::from));
    let r = run(&a);
    assert_eq!(col(&r, "query_name"), ["r1", "r2", "r5", "r6", "r9", "r10"]);
    assert_eq!(col(&r, "tp_tag"), ["P", "P", "P", "S", "P", ""]);
    assert_eq!(col(&r, "InternalPriming_Flag"), ["False", "True", "False", "False", "", ""]);
    // strand '.' (r9) and '*' with no coordinates (r10) are unmapped: empty coordinates and labels
    assert_eq!(col(&r, "target_name"), ["chrT", "chrT", "chrT", "chrT", "", ""]);
    assert_eq!(col(&r, "InternalPriming_ReadEnd3p"), ["100", "98", "310", "600", "", ""]);

    // a stranded row whose coordinates aren't a 0-based interval is an error that names the line
    let bad = tmp("bad.tsv");
    std::fs::write(&bad, "id\tchr\tstr\tbeg\tfin\nr1\tchrT\t+\t50\t100\nr2\tchrT\t+\tx\t98\n").unwrap();
    let mut b = base("aln.tsv");
    b[1] = bad.to_string_lossy().into_owned();
    b.extend(TSV_COLS.map(String::from));
    let err = label_internal_priming_err(&b.iter().map(String::as_str).collect::<Vec<_>>());
    assert!(err.contains("bad.tsv:3") && err.contains("not a valid 0-based interval"), "{err}");
}

#[test]
fn tsv_missing_column_is_an_error() {
    let out = subcommand()
        .args(["-a", &data("aln.tsv"), "-p", &data("tracts.bed"), "-o", "-"])
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--chrom-col"));
}

#[test]
fn annotation_options() {
    let mut a = base("aln.paf");
    a.extend(["--keep-end-nf", "--polya-max-dist", "0", "--tes-max-dist", "60"].map(String::from));
    let r = run(&a);
    // with TX2 kept, r6 (3' end 600) is nearest to TX2's end at 400
    assert_eq!(r[5]["InternalPriming_TxEndSignedDist"], "-200");
    assert_eq!(r[5]["InternalPriming_TxEndTranscriptID"], "TX2");
    assert_eq!(col(&r, "InternalPriming_PolyAWithin0bp"), ["True", "False", "False", "False", "True", "False", "False", "", "False"]);
    assert_eq!(col(&r, "InternalPriming_TxEndWithin60bp"), ["True", "True", "False", "False", "True", "False", "False", "", "True"]);

    // basic only: TX3 is not basic and TX2 is mRNA_end_NF, so '-' reads have no annotated end
    let mut a = base("aln.paf");
    a.push("--gencode-basic-only".into());
    let r = run(&a);
    assert_eq!(r[4]["InternalPriming_TxEndSignedDist"], "");
    assert_eq!(r[4]["InternalPriming_Flag"], "True");
}

#[test]
fn without_annotation_flag_is_polya_adjacency() {
    let (out, err) = label_internal_priming(&["-a", &data("aln.paf"), "-p", &data("tracts.bed"), "-o", "-"]);
    let r = rows(&out);
    assert!(!out.lines().next().unwrap().contains("InternalPriming_TxEnd"));
    assert_eq!(col(&r, "InternalPriming_Flag"), col(&r, "InternalPriming_PolyAWithin3bp"));
    assert!(err.contains("poly(A) adjacency alone"));
}

#[test]
fn chromosome_name_mismatch_warns() {
    let (_, err) = label_internal_priming(&["-a", &data("aln.paf"), "-p", &data("tracts.bed"), "-o", "-"]);
    // r7 is on chrX, absent from the tracts: 1 of the 8 mapped alignments
    assert!(err.contains("12.5% of alignments are on chromosomes absent from the poly(A) tracts"), "{err}");
}

#[test]
fn gzipped_input_and_output() {
    let gz = tmp("out.tsv.gz");
    let mut a = base("aln.paf");
    a[7] = gz.to_string_lossy().into_owned();
    label_internal_priming(&a.iter().map(String::as_str).collect::<Vec<_>>());
    let bytes = std::fs::read(&gz).unwrap();
    assert_eq!(&bytes[..2], &[0x1f, 0x8b]);
    // and read it back as gzipped TSV input
    let (out, _) = label_internal_priming(&[
        "-a", gz.to_str().unwrap(), "--in-format", "tsv", "-p", &data("tracts.bed"), "-o", "-", "--tp-col", "tp_tag",
        "--chrom-col", "target_name", "--strand-col", "strand", "--start-col", "target_start", "--end-col", "target_end",
    ]);
    let r = rows(&out);
    assert_eq!(r.len(), 9);
    assert_eq!(r[7]["target_name"], ""); // the unmapped row stays unmapped
}

#[test]
fn sam_bam_cram_are_rejected() {
    for f in ["sample.sam", "sample.bam", "sample.cram"] {
        let err = label_internal_priming_err(&["-a", f, "-p", &data("tracts.bed"), "-o", "-"]);
        assert!(err.contains("convert it to PAF first"), "{err}");
    }
}

/// Read a Parquet file back as (column names, types, rows of cells rendered like the TSV).
fn read_parquet(path: &Path) -> (Vec<String>, Vec<String>, Vec<Vec<String>>) {
    use arrow_array::cast::AsArray;
    use arrow_array::types::{Int64Type, UInt64Type};
    use arrow_array::Array;
    use arrow_schema::DataType;
    let r = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(std::fs::File::open(path).unwrap())
        .unwrap();
    let schema = r.schema().clone();
    let names = schema.fields().iter().map(|f| f.name().clone()).collect();
    let types = schema.fields().iter().map(|f| f.data_type().to_string()).collect();
    let mut rows = Vec::new();
    for b in r.build().unwrap() {
        let b = b.unwrap();
        for i in 0..b.num_rows() {
            rows.push(
                b.columns()
                    .iter()
                    .map(|a| {
                        if a.is_null(i) {
                            return String::new();
                        }
                        match a.data_type() {
                            DataType::Utf8 => a.as_string::<i32>().value(i).to_string(),
                            DataType::Int64 => a.as_primitive::<Int64Type>().value(i).to_string(),
                            DataType::UInt64 => a.as_primitive::<UInt64Type>().value(i).to_string(),
                            DataType::Boolean => (if a.as_boolean().value(i) { "True" } else { "False" }).to_string(),
                            t => panic!("unexpected type {t}"),
                        }
                    })
                    .collect(),
            );
        }
    }
    (names, types, rows)
}

/// The same rows as `rows()` gives for a TSV.
fn parquet_rows(path: &Path) -> Vec<HashMap<String, String>> {
    let (names, _, rows) = read_parquet(path);
    rows.into_iter().map(|r| names.iter().cloned().zip(r).collect()).collect()
}

#[test]
fn parquet_output_matches_tsv() {
    let pq = tmp("paf.parquet");
    let mut a = base("aln.paf");
    let tsv = run(&a);
    a[7] = pq.to_string_lossy().into_owned();
    label_internal_priming(&a.iter().map(String::as_str).collect::<Vec<_>>());
    let (names, types, _) = read_parquet(&pq);
    assert_eq!(names[..3], ["query_name", "target_name", "target_start"]);
    let ty = |c: &str| types[names.iter().position(|n| n == c).unwrap()].as_str();
    assert_eq!(ty("target_start"), "UInt64");
    assert_eq!(ty("query_len"), "UInt64");
    assert_eq!(ty("target_name"), "Utf8");
    assert_eq!(ty("InternalPriming_TxEndSignedDist"), "Int64");
    assert_eq!(ty("InternalPriming_Flag"), "Boolean");
    assert_eq!(parquet_rows(&pq), tsv); // nulls read back as the TSV's empty cells
}

#[test]
fn parquet_input_round_trip() {
    // a labelled Parquet table as input: only its ID, tp and coordinate columns are read
    let first = tmp("first.parquet");
    let mut a = base("aln.paf");
    a[7] = first.to_string_lossy().into_owned();
    label_internal_priming(&a.iter().map(String::as_str).collect::<Vec<_>>());

    let tracts = data("tracts.bed");
    let cols = [
        "--tp-col", "tp_tag", "--chrom-col", "target_name", "--strand-col", "strand", "--start-col", "target_start",
        "--end-col", "target_end",
    ];
    let again = tmp("again.parquet");
    let mut b = vec!["-a", first.to_str().unwrap(), "-p", &tracts, "-o", again.to_str().unwrap(), "--polya-max-dist", "0"];
    b.extend(cols);
    label_internal_priming(&b);
    let (names, _, _) = read_parquet(&again);
    assert_eq!(names[..12], ["query_name", "target_name", "target_start", "target_end", "strand", "query_len", "query_start", "query_end", "mapq", "tp_tag", "aln_type", "read_3p_clip"]);
    assert_eq!(names[12], "InternalPriming_ReadEnd3p");
    assert!(names.contains(&"InternalPriming_PolyAWithin0bp".to_string()));
    assert!(!names.iter().any(|n| n.starts_with("InternalPriming_TxEnd")), "{names:?}");
    let pr = parquet_rows(&again);
    assert_eq!(col(&pr, "query_name"), ["r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r2"]);
    assert_eq!(col(&pr, "InternalPriming_Flag"), ["True", "False", "False", "False", "True", "False", "False", "", "False"]);
    // its leading columns come back as they were: the table has the record columns under its own names
    let mut c = vec!["-a", first.to_str().unwrap(), "-p", &tracts, "-o", "-"];
    c.extend(cols);
    c.extend(["--query-start-col", "query_start", "--query-end-col", "query_end", "--mapq-col", "mapq"]);
    let tr = rows(&label_internal_priming(&c).0);
    let paf = run(&base("aln.paf"));
    for k in ["query_name", "target_start", "query_len", "query_end", "mapq", "tp_tag", "aln_type", "read_3p_clip"] {
        assert_eq!(col(&tr, k), col(&paf, k), "column {k}");
    }
}

#[test]
fn tsv_input_to_parquet() {
    let pq = tmp("tsv.parquet");
    let mut a = base("aln.tsv");
    a[7] = pq.to_string_lossy().into_owned();
    a.extend(TSV_COLS.map(String::from));
    label_internal_priming(&a.iter().map(String::as_str).collect::<Vec<_>>());
    let (names, types, rows) = read_parquet(&pq);
    assert_eq!(names[..12], ["query_name", "target_name", "target_start", "target_end", "strand", "query_len", "query_start", "query_end", "mapq", "tp_tag", "aln_type", "read_3p_clip"]);
    assert_eq!(types[2], "UInt64");
    assert_eq!(rows.len(), 6);
    assert_eq!(rows[1][0], "r2");
}

#[test]
fn summary_file() {
    let s = tmp("summary.tsv");
    let mut a = base("aln.paf");
    a.extend(["--summary".to_string(), s.to_string_lossy().into_owned()]);
    let (_, err) = label_internal_priming(&a.iter().map(String::as_str).collect::<Vec<_>>());
    // counts are of alignments (PAF lines), so r2's two pieces count twice
    assert!(err.contains("9 alignments (1 unmapped)") && err.contains("unmapped (no labels): 1 (11.1%)"), "{err}");
    let text = std::fs::read_to_string(&s).unwrap();
    // every combination that can occur, zeros included, then the unmapped alignments
    assert_eq!(
        text,
        "InternalPriming_PolyAWithin3bp\tInternalPriming_TxEndWithin50bp\tInternalPriming_Flag\tcount\tpercent\n\
         False\tFalse\tFalse\t3\t33.33\n\
         False\tTrue\tFalse\t0\t0.00\n\
         True\tFalse\tTrue\t2\t22.22\n\
         True\tTrue\tFalse\t3\t33.33\n\
         \t\t\t1\t11.11\n"
    );
    // without an annotation there is no AnnoTranscriptEnd column
    label_internal_priming(&["-a", &data("aln.paf"), "-p", &data("tracts.bed"), "-o", "-", "--summary", s.to_str().unwrap()]);
    let text = std::fs::read_to_string(&s).unwrap();
    assert_eq!(text.lines().next().unwrap(), "InternalPriming_PolyAWithin3bp\tInternalPriming_Flag\tcount\tpercent");
    assert_eq!(text.lines().count(), 4);
    let o = s.to_str().unwrap();
    assert!(label_internal_priming_err(&["-a", &data("aln.paf"), "-p", &data("tracts.bed"), "-o", o, "--summary", o]).contains("same file"));
}

/// Run the subcommand expecting it to fail; returns stderr.
fn label_internal_priming_err(args: &[&str]) -> String {
    let out = subcommand().args(args).output().unwrap();
    assert!(!out.status.success(), "tarotx label-internal-priming {args:?} should have failed");
    String::from_utf8(out.stderr).unwrap()
}

#[test]
fn output_rules() {
    let (paf, tsv, tracts) = (data("aln.paf"), data("aln.tsv"), data("tracts.bed"));
    // -o is always required
    assert!(label_internal_priming_err(&["-a", &paf, "-p", &tracts]).contains("--output"));
    // the table is the only output: there is no tagged or filtered PAF, and no row filtering
    for opt in ["--tagged-out", "--filtered-out", "--exclude-flagged"] {
        assert!(label_internal_priming_err(&["-a", &tsv, "-p", &tracts, "-o", "-", opt, "x.paf"]).contains("unexpected argument"));
    }
    // the GFF3 must have transcript lines: a gene/exon-only GFF3 is refused
    let slim = tmp("slim.gff3");
    let text = std::fs::read_to_string(data("anno.gff3")).unwrap();
    std::fs::write(&slim, text.lines().filter(|l| !l.contains("\ttranscript\t")).collect::<Vec<_>>().join("\n")).unwrap();
    let err = label_internal_priming_err(&["-a", &paf, "-p", &tracts, "-g", slim.to_str().unwrap(), "-o", "-"]);
    assert!(err.contains("has no transcript lines") && err.contains("full GENCODE GFF3"), "{err}");
}
