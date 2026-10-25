# tarotx label-internal-priming: flag likely internal priming in long-read RNA-seq alignments

`tarotx label-internal-priming` labels each long-read alignment by what lies at its 3′ end and flags likely **internal priming**. Internal priming happens when the oligo(dT) primer anneals to an A-rich stretch inside a transcript instead of its poly(A) tail.

For each alignment, label-internal-priming answers two questions:

1. **Genomic poly(A) tract:** is the 3′ end within *X* bp of **either edge** of a same-strand genomic poly(A) tract? The tracts come from [polyscan](https://github.com/maxgmarin/polyscan).
2. **Annotated transcript end** (optional, needs a GFF3): how far is the 3′ end from the nearest same-strand annotated transcript end?

An alignment is flagged when it is next to a poly(A) tract **and not** near an annotated transcript end. Without an annotation, the flag is poly(A) adjacency alone.

- **Input:** a minimap2 PAF, or any TSV or Parquet table with read ID, chromosome, strand, start and end columns, such as a maligno compare table. SAM/BAM/CRAM is not read; convert it to PAF first (see [Inputs](#inputs)).
- **Output:** a per-alignment table, TSV or Parquet, that starts with the same six columns as `label-by-gene`'s table, plus a summary of the counts. Unmapped alignments get a row with empty labels.
- **Rows:** one per input alignment, in input order, as described in [One row per input alignment](../README.md#one-row-per-input-alignment).
- **Tagged PAF:** for PAF input, it can also write the PAF itself, with the labels appended as tags.

Like every tarotx subcommand, it assumes **stranded reads, oriented 5′→3′, with adapters and poly(A) tails trimmed** (e.g. PacBio Iso-Seq/Kinnex reads after `lima` and `isoseq refine`); see [Input reads](../README.md#input-reads-what-tarotx-assumes) for what relies on this. Unstranded or untrimmed reads are not supported yet. Both matter here. Strandedness tells which end of the alignment is the 3′ end. Trimming means that the 3′ end is where the transcript sequence ended, which is the point the poly(A) tracts and annotated ends are measured from.

## Installation

See the [tarotx README](../README.md). This subcommand is `tarotx label-internal-priming`.

---

## Quick start

```bash
# 1. poly(A) tracts for the genome (once)
polyscan -f GRCh38.fa.gz -w 10 -p 80 -n A --mode tracts -o GRCh38.polyA_tracts_W10_P80.bed.gz --index

# 2. Label a PAF: the per-alignment table, the PAF with tags, and the counts
tarotx label-internal-priming -a sample.paf.gz -p GRCh38.polyA_tracts_W10_P80.bed.gz -g gencode.v50.annotation.gff3.gz \
  -o sample.int_priming.tsv.gz --tagged-out sample.int_priming.paf.gz --summary sample.int_priming.summary.tsv
```

A summary goes to stderr. This is the chr21 PAF of the test data:

```
loaded references in 0.1 s
11155 alignments (0 unmapped), 8215 annotated transcript ends, 165921 poly(A) tracts
  AnnoTranscriptEnd_within_50bp:               False                True               total
   GenomicPolyA_within_3bp False        2592 (23.2%)        4566 (40.9%)        7158 (64.2%)
    GenomicPolyA_within_3bp True        3027 (27.1%)          970 (8.7%)        3997 (35.8%)
                           total        5619 (50.4%)        5536 (49.6%)      11155 (100.0%)
internal_priming_flag: 3027 (27.1%)
unmapped (no labels): 0 (0.0%)
done in 0.1 s
```

The counts are of **alignments** (rows), not reads: a split read counts once per piece. With `--summary FILE`, the same counts are also written as a TSV. It has one row per combination of the label columns that can occur, zeros included, then the unmapped alignments. Percent is of all rows:

```
GenomicPolyA_within_3bp	AnnoTranscriptEnd_within_50bp	internal_priming_flag	count	percent
False	False	False	2592	23.24
False	True	False	4566	40.93
True	False	True	3027	27.14
True	True	False	970	8.70
			0	0.00
```

Without `-g`, the `AnnoTranscriptEnd_within` column is left out.

More examples:

```bash
# No annotation: the flag is poly(A) adjacency alone
tarotx label-internal-priming -a sample.paf.gz -p tracts.bed.gz -o sample.int_priming.tsv

# SAM/BAM: convert to PAF first, keeping the unmapped reads
maligno sam2paf -U sample.bam > sample.paf

# A maligno compare table (the default column names are maligno's): Parquet in, Parquet out
tarotx label-internal-priming -a sample.compare.parquet -p tracts.bed.gz -g gencode.gff3.gz -o sample.compare.int_priming.parquet

# Any output can be Parquet: the format follows the -o file name
tarotx label-internal-priming -a sample.paf.gz -p tracts.bed.gz -g gencode.gff3.gz -o sample.int_priming.parquet

# A table with other column names
tarotx label-internal-priming -a aln.tsv -p tracts.bed.gz -o labels.tsv \
  --id-col read --chrom-col chrom --strand-col strand --start-col start --end-col end

# The flagged lines of a tagged PAF
gzip -dc sample.int_priming.paf.gz | awk '$NF == "ip:A:Y"'
```

---

## Inputs

| Input | Option | Format |
|---|---|---|
| Alignments | `-a` | A minimap2 PAF; a TSV table with a header line; or a Parquet table. The format is chosen from the file name (`.paf`, `.tsv`, `.txt`, each optionally `.gz`, or `.parquet`), or set with `--format paf\|tsv\|parquet`. Use `-a -` with `--format` for stdin (not for Parquet). SAM/BAM/CRAM is rejected. |
| poly(A) tracts | `-p` | polyscan `--mode tracts` BED6 (`chrom start end name score strand`), or 5-column `bedtools merge -s -c 5,6 -o max,distinct` output. The layout is detected. `+` = A-run on the plus strand, `-` = T-run (an A-run on the minus strand). |
| Annotation | `-g` (optional) | GENCODE-style GFF3. label-internal-priming uses the `transcript` features and their `transcript_id`, `gene_name` and `tag` attributes. |

All text inputs may be gzip- or BGZF-compressed. Chromosome names must match across files (`chr1` vs `1`); label-internal-priming warns when they don't.

### SAM/BAM input

SAM/BAM is not read. Convert it to PAF first with maligno's `sam2paf`, with `-U` so that unmapped reads are kept as placeholder lines (and so get a row):

```bash
maligno sam2paf -U sample.bam > sample.paf
```

The [SAM/BAM input](../README.md#sambam-input) section of the README covers which records are converted and how. On the HG00099 sample, the converted BAM gives exactly the same rows and labels as minimap2's PAF.

### Rows

Every input alignment gets one row, in input order: every PAF line, and every table row. Nothing is selected among a read's alignments. A maligno compare table already holds one representative alignment per read, so it gives one row per read. The rule and its edge cases are in [One row per input alignment](../README.md#one-row-per-input-alignment).

**Unmapped alignments** get a row with empty coordinates and empty labels. They are never flagged, so `--exclude-flagged` keeps them. An alignment counts as unmapped when:
- **PAF:** the strand is not `+`/`-` (e.g. `*`);
- **tables:** the strand is not `+`/`-`, the chromosome is empty, or the start and end are both empty. A row with a strand whose start and end aren't a 0-based interval is an error that names the row.

---

## Definitions

All coordinates are 0-based, half-open.

**3′ end of an alignment** (`TxAln_3p_end`): where the read stops on the genome. That is the target end for `+` alignments and the target start for `−` alignments.

**Poly(A) adjacency:** each tract contributes two points, its start edge and its end edge. An alignment is adjacent if the nearest same-strand edge is ≤ *X* bp from its 3′ end (`--polya-max-dist`, default 3). This covers both reads that stop at or just before the tract and reads that run through it and stop just past it:

```
                 tract  AAAAAAAAAAAAAA
read ===========>|                         near edge: stops at / just before the tract
read =================================>|   far edge: runs through, stops at / just past it
```

A read that ends in the **middle** of a tract, more than *X* bp from both edges, is **not** adjacent.

**Annotated transcript end:** `end` for `+` transcripts, `start` for `−` transcripts. By default, transcripts tagged `mRNA_end_NF` (3′ end not found) are dropped (`--keep-end-nf` keeps them). `--gencode-basic-only` keeps only transcripts tagged `basic`.

**Signed distance to the annotated end** (`AnnoTranscriptEnd_signed_dist`), in transcript orientation:

- `< 0`: the read stops short of the annotated end.
- `> 0`: the read runs past it.
- `0`: it ends exactly there.

The alignment is near the annotated end if `|signed distance| ≤ D` (`--tes-max-dist`, default 50).

**Ties:** when two points are equally close, the one beyond the read's 3′ end is used. For annotated ends, that gives a negative distance. When several transcripts end at the same position, the first one in the GFF3 is reported.

**Internal priming flag:** `GenomicPolyA_within_{X}bp` and not `AnnoTranscriptEnd_within_{D}bp`. Without an annotation, it is `GenomicPolyA_within_{X}bp` alone. Alignments with no strand (`.`) or on a chromosome without tracts are never flagged.

---

## Output table (`-o`, required)

One row per labelled or unmapped alignment. The format follows the file name:

- `.parquet`: Parquet, zstd-compressed, in row groups of 65,536 rows.
- `.tsv.gz` (any `.gz`): gzipped TSV.
- `-`: TSV on stdout.
- Anything else: plain TSV.

`--exclude-flagged` leaves out flagged rows.

### Leading columns

Every table starts with the same six columns as `label-by-gene`'s, so the two tables can be joined on them:

| Column | Meaning |
|---|---|
| `query_name` | Read name: PAF column 1, or the table's `--id-col` |
| `tp_tag` | `P`, `S`, `I` or `i`: minimap2's `tp:A` tag (PAF), or the table's `--tp-col`; empty when absent |
| `target_name`, `target_start`, `target_end` | Aligned part of the genome; empty when unmapped |
| `strand` | `+` or `-`; empty when unmapped |

For **PAF input**, six columns from the alignment record follow. **Tables** get them too when they have the read length, query span and mapq, as a maligno compare table does (see below):

| Column | Meaning |
|---|---|
| `query_len` | Read length, hard clips included |
| `query_start`, `query_end` | Aligned part of the read, on the read as sequenced (PAF convention) |
| `mapq` | Mapping quality |
| `aln_type` | `primary`, `secondary` or `inversion`, from `tp` (`P`, `S`, `I`/`i`); `unmapped` for an unmapped alignment. minimap2 writes supplementary alignments as `tp:A:P`, so they read `primary` here |
| `read_3p_clip` | Bases of the read past the alignment's 3′ end (`query_len - query_end`). With trimmed reads there should be no tail or adapter left to clip, so a large value points to a split read or an unaligned end. For a split read, the piece with a small value holds the read's real 3′ end. |

**Tables:** every row gets a row, in input order, and the table's other columns are not copied. To bring the labels back to the table, join on the row order, or on `query_name` and the coordinates.

- **Record columns:** a table gets the six record columns when it has all four of these: `--query-len-col` (default `query_len`), `--query-start-col` (`query_start_A`), `--query-end-col` (`query_end_A`) and `--mapq-col` (`mapq_A`). The defaults are maligno's names.
  - The query span must be on the read as sequenced, as in PAF; maligno writes it that way. `aln_type` comes from the tp column.
  - For a maligno compare table, the output is then the same as for the PAF it was made from, `read_3p_clip` included. A compare table holds each read's representative alignment, not every piece of a split read.
- **Partial record columns:** a table with only some of the four gets the six shared columns alone, and a note names the missing ones.

### Label columns

| Column | Meaning |
|---|---|
| `TxAln_3p_end` | 0-based 3′ boundary of the alignment |
| `AnnoTranscriptEnd_signed_dist` | Signed bp to the nearest same-strand annotated end (see above) |
| `AnnoTranscriptEnd_transcript_id` | That annotated transcript |
| `AnnoTranscriptEnd_gene_name` | Its gene |
| `GenomicPolyA_edge_dist` | bp from the 3′ end to the nearest same-strand tract edge |
| `GenomicPolyA_tract_len` | Length of that tract (bp) |
| `GenomicPolyA_score` | That tract's polyscan score (highest window purity, %) |
| `GenomicPolyA_within_{X}bp` | `True` if `GenomicPolyA_edge_dist <= X` |
| `AnnoTranscriptEnd_within_{D}bp` | `True` if `abs(AnnoTranscriptEnd_signed_dist) <= D` |
| `internal_priming_flag` | `True` if poly(A)-adjacent and not near an annotated end |

- `{X}` and `{D}` in the column names follow the thresholds used. The defaults give `GenomicPolyA_within_3bp` and `AnnoTranscriptEnd_within_50bp`.
- The `AnnoTranscriptEnd_*` columns are left out when no annotation is given.
- A value is empty (null in Parquet) when the alignment's chromosome and strand have no tracts or no annotated ends. Every label is empty for an unmapped alignment.
- In Parquet, the distance, length and score columns are `int64`, the ID and name columns are strings, and the `within`/flag columns are booleans. These are the types the Python script writes. The leading and record columns are `uint64`, apart from the name, tp, strand and type columns, which are strings.

---

## Tagged PAF (`--tagged-out`, PAF input)

For PAF input, `--tagged-out` writes every input line, in input order, with these tags appended to the mapped ones. Name the file `.paf`, or `.paf.gz` for gzip; it must differ from the `-o` path.

| Tag | Table column | Meaning |
|---|---|---|
| `pd:i` | `GenomicPolyA_edge_dist` | bp from the 3′ end to the nearest same-strand poly(A) tract edge |
| `td:i` | `AnnoTranscriptEnd_signed_dist` | Signed bp to the nearest same-strand annotated end (`< 0`: stops short). Only with `-g`. |
| `ip:A` | `internal_priming_flag` | `Y` or `N` |

- The tags are appended at the end of each line, and the rest of the line is copied byte for byte. PAF has no header, so this table is where the tags are documented.
- `pd` or `td` is left out when there is no tract or annotated end on that chromosome and strand.
- Any `pd`/`td`/`ip` tags already on a line are replaced, e.g. when a tagged PAF is labelled again.
- Unmapped lines are copied untagged.
- There is no filtered output. To drop the flagged alignments, filter the table, or the tagged PAF on `ip:A:Y`.

The tags are a summary. The full set of per-alignment values is in the `-o` table, which is always written.

---

## Options

| Option | Default | Description |
|---|---|---|
| `-a, --alignments` | required | PAF, TSV or Parquet |
| `-p, --polya-tracts` | required | polyscan poly(A) tracts (BED6 or 5-column) |
| `-g, --gff` | none | GFF3 annotation |
| `--format` | `auto` | `auto`, `paf`, `tsv` or `parquet` |
| `-o, --output` | required | Per-alignment table (`.parquet`, `.tsv`, `.tsv.gz`, or `-` for TSV on stdout) |
| `--tagged-out` | none | PAF input only: every input line, with tags on the mapped ones (`.paf` or `.paf.gz`) |
| `--exclude-flagged` | off | Leave flagged rows out of the table |
| `--summary` | none | Also write the counts as a TSV |
| `--polya-max-dist` | 3 | Max bp from the 3′ end to a tract edge |
| `--tes-max-dist` | 50 | Max bp from the 3′ end to an annotated end |
| `--gencode-basic-only` | off | Use only transcripts tagged `basic` |
| `--keep-end-nf` | off | Keep transcripts tagged `mRNA_end_NF` |
| `--id-col`, `--tp-col`, `--chrom-col`, `--strand-col`, `--start-col`, `--end-col` | `query_name`, `tp_tag_A`, `target_name_A`, `strand_A`, `target_start_A`, `target_end_A` | Table input columns, TSV or Parquet (0-based start, exclusive end). The tp column is optional. |
| `--query-len-col`, `--query-start-col`, `--query-end-col`, `--mapq-col` | `query_len`, `query_start_A`, `query_end_A`, `mapq_A` | Optional table columns; with all four, a table gets the record columns, `read_3p_clip` included |

`-o` is always required, so every run leaves a full record of each call and its supporting values. `--tagged-out` needs PAF input.

---

## Performance

`label-internal-priming` loads the tracts and the annotation fully into memory: about 12 bytes per tract, plus one entry per annotated transcript end. It then streams the alignments one record at a time. Each record gets two binary searches, one over the tract edges and one over the transcript ends of its chromosome and strand, and its row is written at once. Memory therefore depends on the reference files, not on the number of alignments.

Measured on an Apple M5 on mains power, with the GRCh38 tracts (12.3M) and the full GENCODE v50 GFF3, HG00099 Kinnex transcript clusters:

| Input | Rows | Time | Peak memory |
|---|---|---|---|
| PAF, all lines, TSV out | 918,091 | 4.4 s | 282 MB |
| maligno compare Parquet (93 columns), Parquet out | 907,941 | 4.8 s | 331 MB |
| maligno compare Parquet (93 columns), TSV.gz out | 907,941 | 4.7 s | 294 MB |

About 3 s of each run is loading the tracts and the annotation, which are read in parallel. From a table, only the columns used are read. On battery power, expect about twice as long.

The Python script on the same compare Parquet takes 86 s and 10.5 GB.
