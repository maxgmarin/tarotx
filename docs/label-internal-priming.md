# tarotx label-internal-priming: flag likely internal priming in long-read RNA-seq alignments

`tarotx label-internal-priming` labels each long-read alignment by what lies at its 3′ end and flags likely **internal priming**. Internal priming happens when the oligo(dT) primer anneals to an A-rich stretch inside a transcript instead of its poly(A) tail.

For each alignment, it answers two questions:

1. **Genomic poly(A) tract:** is the 3′ end within *X* bp of **either edge** of a same-strand genomic poly(A) tract? The tracts come from [polyscan](https://github.com/maxgmarin/polyscan).
2. **Annotated transcript end** (optional, needs the GENCODE GFF3): how far is the 3′ end from the nearest same-strand annotated transcript end?

An alignment is flagged (`InternalPriming_Flag`) when it is next to a poly(A) tract **and not** near an annotated transcript end. Without an annotation, the flag is poly(A) adjacency alone.

- **Input:** a minimap2 PAF, or a TSV or Parquet table such as a maligno compare table.
- **Output:** one row per input alignment, in input order. Each row has the 12 leading columns every tarotx subcommand writes, then the `InternalPriming_*` labels. A summary of the counts goes to stderr, and to a TSV with `--summary`.
- **Shared behaviour:** the input options, the leading columns and the one-row-per-alignment rule are the same as `label-by-gene`'s; see [What every subcommand shares](../README.md#what-every-subcommand-shares) and [One row per input alignment](../README.md#one-row-per-input-alignment).

Like every tarotx subcommand, it assumes **stranded reads, oriented 5′→3′, with adapters and poly(A) tails trimmed** (e.g. PacBio Iso-Seq/Kinnex reads after `lima` and `isoseq refine`); see [Input reads](../README.md#input-reads-what-tarotx-assumes). Both matter here:
- **Strandedness** tells which end of the alignment is the 3′ end.
- **Trimming** means that the 3′ end is where the transcript sequence ended. That is the point the poly(A) tracts and annotated ends are measured from.

## Installation

See the [tarotx README](../README.md). This subcommand is `tarotx label-internal-priming`.

---

## Quick start

```bash
# 1. poly(A) tracts for the genome (once)
polyscan -f GRCh38.fa.gz -w 10 -p 80 -n A --mode tracts -o GRCh38.polyA_tracts_W10_P80.bed.gz --index

# 2. Label a PAF, with the counts as a TSV
tarotx label-internal-priming -a sample.paf.gz -p GRCh38.polyA_tracts_W10_P80.bed.gz \
  -g gencode.v50.primary_assembly.annotation.gff3.gz -o sample.int_priming.parquet --summary sample.int_priming.summary.tsv

# A maligno compare table (the default column names are maligno's)
tarotx label-internal-priming -a sample.compare.parquet -p tracts.bed.gz -g gencode.gff3.gz -o sample.int_priming.parquet

# No annotation: the flag is poly(A) adjacency alone
tarotx label-internal-priming -a sample.paf.gz -p tracts.bed.gz -o sample.int_priming.tsv.gz

# SAM/BAM: convert to PAF first, keeping the unmapped reads
maligno sam2paf -U sample.bam > sample.paf
```

---

## Inputs

| Input | Option | Format |
|---|---|---|
| Alignments | `-a` | A minimap2 PAF, or a TSV or Parquet table; see [the shared options](../README.md#what-every-subcommand-shares) |
| poly(A) tracts | `-p` | polyscan `--mode tracts` BED6 (`chrom start end name score strand`), or 5-column `bedtools merge -s -c 5,6 -o max,distinct` output. The layout is detected. `+` = A-run on the plus strand, `-` = T-run (an A-run on the minus strand). |
| Annotation | `-g` (optional) | The full GENCODE GFF3. Its `transcript` lines and their `transcript_id`, `gene_name` and `tag` attributes give the annotated 3′ ends. A GFF3 without transcript lines (e.g. a gene/exon-only one), or with none left after the filters, is an error. |

All text inputs may be gzip- or BGZF-compressed. Chromosome names must match across files (`chr1` vs `1`); a warning says when they don't.

---

## How labelling works

All coordinates are 0-based, half-open.

**3′ end of an alignment** (`InternalPriming_ReadEnd3p`): where the read stops on the genome. That is the target end for `+` alignments and the target start for `−` alignments.

**Poly(A) adjacency:** each tract contributes two points, its start edge and its end edge. An alignment is adjacent if the nearest same-strand edge is ≤ *X* bp from its 3′ end (`--polya-max-dist`, default 3). This covers both reads that stop at or just before the tract and reads that run through it and stop just past it:

```
                 tract  AAAAAAAAAAAAAA
read ===========>|                         near edge: stops at / just before the tract
read =================================>|   far edge: runs through, stops at / just past it
```

A read that ends in the **middle** of a tract, more than *X* bp from both edges, is **not** adjacent.

**Annotated transcript end:** `end` for `+` transcripts, `start` for `−` transcripts.
- By default, transcripts tagged `mRNA_end_NF` (3′ end not found) are dropped; `--keep-end-nf` keeps them.
- `--gencode-basic-only` keeps only transcripts tagged `basic`.

**Signed distance to the annotated end** (`InternalPriming_TxEndSignedDist`), in transcript orientation:
- `< 0`: the read stops short of the annotated end.
- `> 0`: the read runs past it.
- `0`: it ends exactly there.

The alignment is near the annotated end if `|signed distance| ≤ D` (`--tes-max-dist`, default 50).

**Ties:** when two points are equally close, the one beyond the read's 3′ end is used. For annotated ends, that gives a negative distance. When several transcripts end at the same position, the first one in the GFF3 is reported.

**The flag:** `InternalPriming_Flag` is `InternalPriming_PolyAWithin{X}bp` and not `InternalPriming_TxEndWithin{D}bp`. Without an annotation, it is `InternalPriming_PolyAWithin{X}bp` alone. An alignment on a chromosome without tracts is never flagged.

**Unmapped alignments** get a row with every label empty, and are never flagged.

---

## Label columns

After the [12 leading columns](../README.md#what-every-subcommand-shares):

| Column | Type | Meaning |
|---|---|---|
| `InternalPriming_ReadEnd3p` | int64 | 0-based 3′ boundary of the alignment |
| `InternalPriming_TxEndSignedDist` | int64 | Signed bp to the nearest same-strand annotated end (see above) |
| `InternalPriming_TxEndTranscriptID` | str | That annotated transcript |
| `InternalPriming_TxEndGeneName` | str | Its gene |
| `InternalPriming_PolyAEdgeDist` | int64 | bp from the 3′ end to the nearest same-strand tract edge |
| `InternalPriming_PolyATractLen` | int64 | Length of that tract (bp) |
| `InternalPriming_PolyAScore` | int64 | That tract's polyscan score (highest window purity, %) |
| `InternalPriming_PolyAWithin{X}bp` | bool | `True` if `InternalPriming_PolyAEdgeDist <= X` |
| `InternalPriming_TxEndWithin{D}bp` | bool | `True` if `abs(InternalPriming_TxEndSignedDist) <= D` |
| `InternalPriming_Flag` | bool | `True` if poly(A)-adjacent and not near an annotated end |

- `{X}` and `{D}` in the column names follow the thresholds used. The defaults give `InternalPriming_PolyAWithin3bp` and `InternalPriming_TxEndWithin50bp`.
- The four `InternalPriming_TxEnd*` columns are left out when no annotation is given.
- A value is empty (null in Parquet) when the alignment's chromosome and strand have no tracts or no annotated ends. Every label is empty for an unmapped alignment.
- In TSV, booleans are written as `True`/`False`, as pandas writes them.

To keep only the alignments that aren't flagged, filter the table, e.g. in pandas: `t[t["InternalPriming_Flag"] != True]`. That keeps unmapped rows too, whose flag is empty.

---

## Summary

At the end of a run, the counts go to stderr. This is the chr21 PAF of the test data:

```
loaded references in 0.0 s
11155 alignments (0 unmapped), 8215 annotated transcript ends, 165921 poly(A) tracts
    InternalPriming_TxEndWithin50bp:               False                True               total
InternalPriming_PolyAWithin3bp False        2592 (23.2%)        4566 (40.9%)        7158 (64.2%)
 InternalPriming_PolyAWithin3bp True        3027 (27.1%)          970 (8.7%)        3997 (35.8%)
                               total        5619 (50.4%)        5536 (49.6%)      11155 (100.0%)
InternalPriming_Flag: 3027 (27.1%)
unmapped (no labels): 0 (0.0%)
done in 0.1 s
```

The counts are of **alignments** (rows), not reads: a split read counts once per piece. With `--summary FILE`, the same counts are also written as a TSV. It has one row per combination of the label columns that can occur, zeros included, then the unmapped alignments. Percent is of all rows:

```
InternalPriming_PolyAWithin3bp	InternalPriming_TxEndWithin50bp	InternalPriming_Flag	count	percent
False	False	False	2592	23.24
False	True	False	4566	40.93
True	False	True	3027	27.14
True	True	False	970	8.70
			0	0.00
```

Without `-g`, the `InternalPriming_TxEndWithin` column is left out.

---

## Options

The input, output and table-column options are [shared with every subcommand](../README.md#what-every-subcommand-shares). Its own options:

| Option | Default | Description |
|---|---|---|
| `-p, --polya-tracts` | required | polyscan poly(A) tracts (BED6 or 5-column) |
| `-g, --gff` | none | The full GENCODE GFF3 |
| `--polya-max-dist` | 3 | Max bp from the 3′ end to a tract edge |
| `--tes-max-dist` | 50 | Max bp from the 3′ end to an annotated end |
| `--gencode-basic-only` | off | Use only transcripts tagged `basic` |
| `--keep-end-nf` | off | Keep transcripts tagged `mRNA_end_NF` |

---

## Performance

The tracts and the annotation are loaded fully into memory: about 12 bytes per tract, plus one entry per annotated transcript end. The alignments are then streamed one at a time. Each gets two binary searches, one over the tract edges and one over the transcript ends of its chromosome and strand, and its row is written at once. Memory therefore depends on the reference files, not on the number of alignments.

Measured on an Apple M5 on mains power, with the GRCh38 tracts (12.3M) and the full GENCODE v50 GFF3, HG00099 Kinnex transcript clusters:

| Input | Rows | Time | Peak memory |
|---|---|---|---|
| PAF, Parquet out | 918,091 | 4.4 s | 324 MB |
| maligno compare Parquet, Parquet out | 907,941 | 4.3 s | 348 MB |
| maligno compare Parquet, TSV.gz out | 907,941 | 5.1 s | 299 MB |

About 3 s of each run is loading the tracts and the annotation, which are read in parallel. From a table, only the columns used are read. On battery power, expect about twice as long.
