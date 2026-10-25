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
    assert_eq!(col(&r, "TxAln_3p_end"), ["100", "98", "215", "233", "310", "600", "10", "", "198"]);
    assert_eq!(col(&r, "GenomicPolyA_edge_dist"), ["0", "2", "15", "3", "0", "290", "", "", "2"]);
    // r3 ends 15 bp from both edges of [200,230): the edge beyond its 3' end wins, so the tract is that one
    assert_eq!(col(&r, "GenomicPolyA_tract_len"), ["10", "10", "30", "30", "10", "10", "", "", "30"]);
    assert_eq!(col(&r, "GenomicPolyA_score"), ["90", "90", "100", "100", "80", "80", "", "", "100"]);
    assert_eq!(col(&r, "AnnoTranscriptEnd_signed_dist"), ["-50", "-52", "65", "83", "-30", "-320", "", "", "48"]);
    // TX2 is mRNA_end_NF, so dropped: the '-' reads see TX3
    assert_eq!(col(&r, "AnnoTranscriptEnd_transcript_id"), ["TX1", "TX1", "TX1", "TX1", "TX3", "TX3", "", "", "TX1"]);
    assert_eq!(col(&r, "GenomicPolyA_within_3bp"), ["True", "True", "False", "True", "True", "False", "False", "", "True"]);
    assert_eq!(col(&r, "AnnoTranscriptEnd_within_50bp"), ["True", "False", "False", "False", "True", "False", "False", "", "True"]);
    assert_eq!(col(&r, "internal_priming_flag"), ["False", "True", "False", "True", "False", "False", "False", "", "False"]);
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
fn same_leading_columns_as_label_by_gene() {
    let lead = ["query_name", "tp_tag", "target_name", "target_start", "target_end", "strand"];
    let record = ["query_len", "query_start", "query_end", "mapq", "aln_type", "read_3p_clip"];
    let header = |args: &[String]| -> Vec<String> {
        let (out, _) = label_internal_priming(&args.iter().map(String::as_str).collect::<Vec<_>>());
        out.lines().next().unwrap().split('\t').map(String::from).collect()
    };
    let h = header(&base("aln.paf"));
    assert_eq!(h[..6], lead);
    assert_eq!(h[6..12], record);
    assert_eq!(h[12], "TxAln_3p_end");
    let mut a = base("aln.tsv");
    a.extend(TSV_COLS.map(String::from));
    let h = header(&a);
    assert_eq!(h[..6], lead);
    assert_eq!(h[6], "TxAln_3p_end"); // a table's other columns are not copied, and it has no record columns
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
        for c in ["query_len", "query_start", "query_end", "mapq", "aln_type", "read_3p_clip", "TxAln_3p_end", "internal_priming_flag"] {
            assert_eq!(row[c], paf[i][c], "column {c}, read {}", row["query_name"]);
        }
    }
    assert_eq!(col(&tab, "read_3p_clip"), ["10", "3", ""]);
    assert_eq!(col(&tab, "aln_type"), ["primary", "primary", "unmapped"]);

    // with only some of the four columns, they are left out, and a note says which are missing
    let mut b = a.clone();
    b.extend(["--mapq-col", "no_such_column"].map(String::from));
    let (out, err) = label_internal_priming(&b.iter().map(String::as_str).collect::<Vec<_>>());
    assert!(err.contains("no column [\"no_such_column\"]") && err.contains("--mapq-col"), "{err}");
    assert_eq!(out.lines().next().unwrap().split('\t').nth(6), Some("TxAln_3p_end"));
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
    assert_eq!(col(&r, "internal_priming_flag"), ["False", "True", "False", "False", "", ""]);
    // strand '.' (r9) and '*' with no coordinates (r10) are unmapped: empty coordinates and labels
    assert_eq!(col(&r, "target_name"), ["chrT", "chrT", "chrT", "chrT", "", ""]);
    assert_eq!(col(&r, "TxAln_3p_end"), ["100", "98", "310", "600", "", ""]);

    // --exclude-flagged leaves out r2; unmapped rows are not flagged, so they stay
    a.push("--exclude-flagged".into());
    assert_eq!(col(&run(&a), "query_name"), ["r1", "r5", "r6", "r9", "r10"]);

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
    assert_eq!(r[5]["AnnoTranscriptEnd_signed_dist"], "-200");
    assert_eq!(r[5]["AnnoTranscriptEnd_transcript_id"], "TX2");
    assert_eq!(col(&r, "GenomicPolyA_within_0bp"), ["True", "False", "False", "False", "True", "False", "False", "", "False"]);
    assert_eq!(col(&r, "AnnoTranscriptEnd_within_60bp"), ["True", "True", "False", "False", "True", "False", "False", "", "True"]);

    // basic only: TX3 is not basic and TX2 is mRNA_end_NF, so '-' reads have no annotated end
    let mut a = base("aln.paf");
    a.push("--gencode-basic-only".into());
    let r = run(&a);
    assert_eq!(r[4]["AnnoTranscriptEnd_signed_dist"], "");
    assert_eq!(r[4]["internal_priming_flag"], "True");
}

#[test]
fn without_annotation_flag_is_polya_adjacency() {
    let (out, err) = label_internal_priming(&["-a", &data("aln.paf"), "-p", &data("tracts.bed"), "-o", "-"]);
    let r = rows(&out);
    assert!(!out.lines().next().unwrap().contains("AnnoTranscriptEnd"));
    assert_eq!(col(&r, "internal_priming_flag"), col(&r, "GenomicPolyA_within_3bp"));
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
        "-a", gz.to_str().unwrap(), "--format", "tsv", "-p", &data("tracts.bed"), "-o", "-", "--tp-col", "tp_tag",
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
    assert_eq!(names[..3], ["query_name", "tp_tag", "target_name"]);
    let ty = |c: &str| types[names.iter().position(|n| n == c).unwrap()].as_str();
    assert_eq!(ty("target_start"), "UInt64");
    assert_eq!(ty("query_len"), "UInt64");
    assert_eq!(ty("target_name"), "Utf8");
    assert_eq!(ty("AnnoTranscriptEnd_signed_dist"), "Int64");
    assert_eq!(ty("internal_priming_flag"), "Boolean");
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
    assert_eq!(names[..7], ["query_name", "tp_tag", "target_name", "target_start", "target_end", "strand", "TxAln_3p_end"]);
    assert!(names.contains(&"GenomicPolyA_within_0bp".to_string()));
    assert!(!names.iter().any(|n| n.starts_with("Anno") || n == "query_len"), "{names:?}");
    let pr = parquet_rows(&again);
    assert_eq!(col(&pr, "query_name"), ["r1", "r2", "r3", "r4", "r5", "r6", "r7", "r8", "r2"]);
    assert_eq!(col(&pr, "internal_priming_flag"), ["True", "False", "False", "False", "True", "False", "False", "", "False"]);

    // Parquet in, TSV out, and --exclude-flagged
    let mut c = vec!["-a", first.to_str().unwrap(), "-p", &tracts, "-o", "-", "--exclude-flagged"];
    c.extend(cols);
    let tr = rows(&label_internal_priming(&c).0);
    // no annotation: the flag is poly(A) adjacency within 3 bp; the unmapped r8 stays
    assert_eq!(col(&tr, "query_name"), ["r3", "r6", "r7", "r8"]);
}

#[test]
fn tsv_input_to_parquet() {
    let pq = tmp("tsv.parquet");
    let mut a = base("aln.tsv");
    a[7] = pq.to_string_lossy().into_owned();
    a.extend(TSV_COLS.map(String::from));
    label_internal_priming(&a.iter().map(String::as_str).collect::<Vec<_>>());
    let (names, types, rows) = read_parquet(&pq);
    assert_eq!(names[..6], ["query_name", "tp_tag", "target_name", "target_start", "target_end", "strand"]);
    assert_eq!(types[3], "UInt64");
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
        "GenomicPolyA_within_3bp\tAnnoTranscriptEnd_within_50bp\tinternal_priming_flag\tcount\tpercent\n\
         False\tFalse\tFalse\t3\t33.33\n\
         False\tTrue\tFalse\t0\t0.00\n\
         True\tFalse\tTrue\t2\t22.22\n\
         True\tTrue\tFalse\t3\t33.33\n\
         \t\t\t1\t11.11\n"
    );
    // without an annotation there is no AnnoTranscriptEnd column
    label_internal_priming(&["-a", &data("aln.paf"), "-p", &data("tracts.bed"), "-o", "-", "--summary", s.to_str().unwrap()]);
    let text = std::fs::read_to_string(&s).unwrap();
    assert_eq!(text.lines().next().unwrap(), "GenomicPolyA_within_3bp\tinternal_priming_flag\tcount\tpercent");
    assert_eq!(text.lines().count(), 4);
    let o = s.to_str().unwrap();
    assert!(label_internal_priming_err(&["-a", &data("aln.paf"), "-p", &data("tracts.bed"), "-o", o, "--summary", o]).contains("same file"));
}

#[test]
fn tagged_paf() {
    let tagged = tmp("tagged.paf");
    let t = tagged.to_str().unwrap();
    let mut a = base("aln.paf");
    a.extend(["--tagged-out", t].map(String::from));
    let table = run(&a);

    let input = std::fs::read_to_string(data("aln.paf")).unwrap();
    let out = std::fs::read_to_string(&tagged).unwrap();
    assert_eq!(out.lines().count(), 9);
    for ((inp, line), row) in input.lines().zip(out.lines()).zip(&table) {
        // an unmapped line is copied as it is
        if row["target_name"].is_empty() {
            assert_eq!(line, inp);
            continue;
        }
        // the input line is kept as it is, and the tags are appended
        assert!(line.starts_with(&format!("{inp}\t")), "{line}");
        let tags: HashMap<&str, &str> = line[inp.len() + 1..].split('\t').map(|t| (&t[..2], &t[5..])).collect();
        assert_eq!(tags.get("pd").copied().unwrap_or(""), row["GenomicPolyA_edge_dist"]);
        assert_eq!(tags.get("td").copied().unwrap_or(""), row["AnnoTranscriptEnd_signed_dist"]);
        assert_eq!(tags["ip"], if row["internal_priming_flag"] == "True" { "Y" } else { "N" });
    }
    assert!(out.lines().nth(1).unwrap().ends_with("\tpd:i:2\ttd:i:-52\tip:A:Y"));

    // gzipped when the name ends in .gz
    let gz = tmp("tagged.paf.gz");
    let mut g = base("aln.paf");
    g.extend(["--tagged-out".to_string(), gz.to_string_lossy().into_owned()]);
    run(&g);
    let mut text = String::new();
    std::io::Read::read_to_string(&mut flate2::read::MultiGzDecoder::new(std::fs::File::open(&gz).unwrap()), &mut text)
        .unwrap();
    assert_eq!(text, out);

    // a tagged PAF as input: its old tags are replaced, not duplicated
    let again = tmp("again.paf");
    let mut b = base(t);
    b[1] = t.to_string();
    b.extend(["--tagged-out", again.to_str().unwrap(), "--polya-max-dist", "0"].map(String::from));
    run(&b);
    let again = std::fs::read_to_string(&again).unwrap();
    let second = again.lines().nth(1).unwrap();
    assert_eq!(second.matches("ip:A:").count(), 1);
    assert!(second.ends_with("\tpd:i:2\ttd:i:-52\tip:A:N"), "{second}");
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
    assert!(label_internal_priming_err(&["-a", &paf, "-p", &tracts, "--tagged-out", "x.paf"]).contains("--output"));
    // the tagged output is a PAF, for PAF input, with a path of its own
    assert!(label_internal_priming_err(&["-a", &paf, "-p", &tracts, "-o", "-", "--tagged-out", "x.bam"]).contains("name it .paf or .paf.gz"));
    assert!(label_internal_priming_err(&["-a", &tsv, "-p", &tracts, "-o", "-", "--tagged-out", "x.paf"]).contains("needs PAF input"));
    assert!(label_internal_priming_err(&["-a", &paf, "-p", &tracts, "-o", "x.paf", "--tagged-out", "x.paf"]).contains("are the same file"));
    // there is no filtered output: filter the table, or the tagged PAF on ip:A
    assert!(label_internal_priming_err(&["-a", &paf, "-p", &tracts, "-o", "-", "--filtered-out", "x.paf"]).contains("unexpected argument"));
}
