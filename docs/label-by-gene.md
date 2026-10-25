# tarotx label-by-gene: assign stranded long-read RNA alignments to annotated genes

`tarotx label-by-gene` writes one classification row per long-read RNA alignment, naming the annotated gene the alignment belongs to. Each row answers two questions, and the answers can be checked from the numbers in the same row:

1. **`GeneAssignment_Status`: which gene, and how clear is that choice?** `unique_gene`, `best_of_multiple`, `tied`, `nearby_gene`, `antisense_only`, `intergenic` or `unmapped`.
2. **`GeneAssignment_SupportLevel`: how is that gene supported?** `known_junctions`, `exonic` or `intronic_only`.

Like every tarotx subcommand, it assumes **stranded reads, oriented 5′→3′, with adapters and poly(A) tails trimmed** (e.g. PacBio Iso-Seq/Kinnex reads after `lima` and `isoseq refine`); see [Input reads](../README.md#input-reads-what-tarotx-assumes) for what relies on this. Unstranded or untrimmed reads are not supported yet. For gene assignment, strandedness is what matters: the alignment strand is taken as the transcript strand, so only genes on that strand can be assigned.

- **Input:** a minimap2 PAF, or a TSV or Parquet table of alignments, such as a maligno compare table.
- **Output:** the classification table, the same 21 columns for every input, plus a summary of the counts.
- **Nothing is filtered out:** every input alignment gets one row, in input order, and you choose what to keep afterwards. See [One row per input alignment](../README.md#one-row-per-input-alignment) for split reads, secondary alignments and unmapped reads.

## Installation

See the [tarotx README](../README.md). This subcommand is `tarotx label-by-gene`.

---

## Quick start

```bash
# A PAF (needs cs:Z or cg:Z tags, i.e. minimap2 --cs or -c)
tarotx label-by-gene -a sample.paf.gz -g gencode.v50.primary_assembly.genes_exons.slim.gff3.gz \
  -o sample.gene_assignment.parquet --summary sample.gene_assignment.summary.tsv

# A maligno compare table (the default column names are maligno's; introns are walked from cs_A)
tarotx label-by-gene -a sample.compare.parquet -g genes.gff3.gz -o sample.gene_assignment.parquet

# Another tool's table, with a CIGAR column
tarotx label-by-gene -a sample.tsv -g genes.gff3.gz -o sample.gene_assignment.tsv --id-col read --chrom-col chrom \
  --strand-col strand --start-col start --end-col end --cigar-col cigar

# Leave out readthrough genes, and keep only protein-coding and lncRNA genes
tarotx label-by-gene -a sample.paf.gz -g genes.gff3.gz -o sample.gene_assignment.tsv.gz --skip-readthrough --gene-types protein_coding,lncRNA
```

A BAM or SAM file is converted to a PAF first, with maligno's `sam2paf` and `-U` to keep the unmapped reads (see [SAM/BAM input](../README.md#sambam-input)):

```bash
maligno sam2paf -U sample.bam > sample.paf
```

### The summary

At the end of a run, the counts are printed to stderr: every status, split by support level, as counts and percent of all alignments. They count alignments (rows), not reads, so a split read in a PAF counts once per piece. With `--summary FILE` the same counts are also written as a TSV, one row per combination that can occur. This is the chr21 PAF of the test data:

```
11155 alignments labelled, 1133 genes, 9101 annotated introns
                     known_junctions            exonic     intronic_only           no gene             total
unique_gene             5535 (49.6%)      2392 (21.4%)      2113 (18.9%)          0 (0.0%)     10040 (90.0%)
best_of_multiple          373 (3.3%)         35 (0.3%)          0 (0.0%)          0 (0.0%)        408 (3.7%)
tied                       44 (0.4%)         19 (0.2%)        222 (2.0%)          0 (0.0%)        285 (2.6%)
nearby_gene                 0 (0.0%)          0 (0.0%)          0 (0.0%)        140 (1.3%)        140 (1.3%)
antisense_only              0 (0.0%)          0 (0.0%)          0 (0.0%)        201 (1.8%)        201 (1.8%)
intergenic                  0 (0.0%)          0 (0.0%)          0 (0.0%)         81 (0.7%)         81 (0.7%)
unmapped                    0 (0.0%)          0 (0.0%)          0 (0.0%)          0 (0.0%)          0 (0.0%)
```

---

## Inputs

| Input | Option | Format |
|---|---|---|
| Alignments | `-a` | A minimap2 PAF (`.paf`, `.paf.gz`), or a table with a header line and a cs or CIGAR column: TSV (`.tsv`, `.txt`, optionally `.gz`) or Parquet (`.parquet`, zstd- or snappy-compressed). The format comes from the file name, or `--format paf\|tsv\|parquet`. Use `-a -` with `--format` for stdin (not for Parquet). |
| Annotation | `-g` | GENCODE-style GFF3, gzipped or not. Only `gene` and `exon` lines are read, so a slimmed GFF works as well as the full one. |

Chromosome names must match between the alignments and the annotation (`chr1` vs `1`); label-by-gene warns when they don't.

### PAF

- **Rows:** every line gets a row, in input order, including secondary and supplementary alignments, and unmapped lines (strand `*`).
- **Introns:** from the `cs:Z` tag's `~` operations, or from the `cg:Z` CIGAR when there's no `cs` tag. A line with neither is an error.
- **Span check:** the `cs` or `cg` must cover exactly the alignment's target span.

### Tables

From a table, only these columns are read:

| Option | Default | Column |
|---|---|---|
| `--id-col` | `query_name` | Read ID |
| `--tp-col` | `tp_tag_A` | Alignment type (`P`, `S`, `I`, `i`), optional |
| `--chrom-col` | `target_name_A` | Chromosome |
| `--strand-col` | `strand_A` | Strand (`+`/`-`; anything else is unmapped) |
| `--start-col` | `target_start_A` | 0-based start |
| `--end-col` | `target_end_A` | Exclusive end |
| `--cs-col` | `cs_A` | The alignment's minimap2 cs tag |
| `--cigar-col` | none | The alignment's CIGAR, used instead of `--cs-col` |

- **Rows:** every row gets a row, in input order. A row whose strand isn't `+`/`-`, or whose start/end is empty, is `unmapped`.
- **Introns:** always walked from the cs tag (`~` operations) or the CIGAR (`N` operations), as for a PAF. Precomputed junction columns, such as maligno's `genomic_junctions_A`, are not read; on the HG00099 compare table they agree with `cs_A` on every row.
- **Span check:** for a mapped row, the cs or CIGAR must cover exactly the span. Anything else is an error that names the row.

---

## How classification works

All coordinates are 0-based, half-open.

### The annotation

- **Genes:** from `gene` lines, with their span, strand, `gene_id`, `gene_name` and `gene_type`.
  - `--gene-types` keeps only the listed types.
  - `--skip-readthrough` leaves out genes tagged `readthrough_gene`. GENCODE v50 has 665 of them, joining two neighbouring genes, e.g. `CENPS-CORT`.
- **Merged exons per gene:** the exons of all the gene's transcripts. Exons that overlap or touch become one interval.
- **Annotated introns per gene:** the gaps between consecutive exons of each transcript, de-duplicated.
- Each exon names its gene through `gene_id` and its transcript through `Parent`, as GENCODE does.

### One alignment

The **aligned blocks** are the alignment span minus its introns:

```
alignment span   |===================================================|
read introns           ^^^^^^^^^              ^^^^^^^^^^^^
aligned blocks   |=====|        |=============|           |==========|
```

1. **Touched genes:** a gene is touched when its body overlaps an aligned block. A gene lying entirely inside one of the read's introns isn't touched.
2. **Two numbers for each touched same-strand gene:**
   - **matched introns:** read introns identical, in both coordinates, to one of the gene's annotated introns;
   - **exonic fraction:** the share of the read's aligned bases on the gene's merged exons.
3. **Support level of each gene:**

   | Support level | Rule |
   |---|---|
   | `known_junctions` | ≥ 1 matched intron |
   | `exonic` | No matched intron, and an exonic fraction ≥ `--min-exonic-frac` (0.2) |
   | `intronic_only` | Neither: the read lies in the gene's body, perhaps grazing an exon |

4. **Ranking.** Same-strand genes are ordered by:
   1. support level (`known_junctions` > `exonic` > `intronic_only`);
   2. more matched introns;
   3. more exonic fraction;
   4. `gene_id`, which only fixes the order.
5. **Competitors:** the genes with `known_junctions` or `exonic` support. If there are none, all the touched genes. A gene the read only grazes doesn't compete with one that has exonic support.
6. **The best gene(s):** the top gene, and every competitor tied with it.
   - With exonic support, a tie means the same number of matched introns and an exonic fraction less than `--tie-margin` (0.1) behind. So the top gene must lead by at least 10% of the read to win on exonic overlap alone. An exact tie is always a tie.
   - At `intronic_only`, all genes are tied.
7. **Status:**

   | Status | Rule |
   |---|---|
   | `unique_gene` | One competitor: only one gene is compatible |
   | `best_of_multiple` | Several competitors, and the top one has no tie |
   | `tied` | The top gene ties with one or more others. The first `GeneAssignment_NumBestGenes` entries of `GeneAssignment_Candidates` are the tied genes. |

8. **No same-strand gene touched.** One of these, checked in order, with the support level empty:

   | Status | Rule |
   |---|---|
   | `nearby_gene` | A same-strand gene within `--near-max-dist` (5,000 bp) of an aligned block, e.g. an extended 3′ UTR; reported with its distance and side. Equal distances go to the smaller `gene_id`. |
   | `antisense_only` | Only opposite-strand genes are touched |
   | `intergenic` | Nothing is touched, and no same-strand gene is near |
   | `unmapped` | No alignment |

A mono-exonic read has no introns, so the best support it can have is `exonic`.

---

## Output table (`-o`, required)

One row per alignment, the same 21 columns whatever the input. The format follows the file name:
- `.parquet`: zstd-compressed, in row groups of 65,536 rows;
- `.tsv.gz` (any `.gz`): gzipped TSV;
- `-`: TSV on stdout;
- anything else: plain TSV.

| # | Column | Type | Meaning |
|---|---|---|---|
| 1 | `query_name` | str | Read name: PAF column 1, or the table's `--id-col` |
| 2 | `tp_tag` | str | `P`, `S`, `I` or `i`: PAF's `tp:A:` tag, or the table's `--tp-col`; empty when absent |
| 3 | `target_name` | str | Chromosome; empty when unmapped |
| 4 | `target_start` | uint64 | 0-based start; empty when unmapped |
| 5 | `target_end` | uint64 | Exclusive end; empty when unmapped |
| 6 | `strand` | str | `+` or `-`; empty when unmapped |
| 7 | `GeneAssignment_Status` | str | Status (above) |
| 8 | `GeneAssignment_SupportLevel` | str | Support level (above); empty when no same-strand gene is touched |
| 9 | `GeneAssignment_GeneID` | str | The chosen gene: the top-ranked one, or the nearest for `nearby_gene` |
| 10 | `GeneAssignment_GeneName` | str | Its name, or its ID when it has none |
| 11 | `GeneAssignment_GeneType` | str | Its `gene_type` |
| 12 | `GeneAssignment_IntronsMatched` | int64 | Read introns matching the chosen gene's annotated introns |
| 13 | `GeneAssignment_ReadIntrons` | int64 | Introns in the read |
| 14 | `GeneAssignment_ExonicFrac` | float64 | Share of the read's aligned bp on the chosen gene's exons (3 decimals) |
| 15 | `GeneAssignment_AlignedBp` | int64 | The read's aligned bp (span minus introns) |
| 16 | `GeneAssignment_NumSameStrandGenes` | int64 | Same-strand genes the read touches |
| 17 | `GeneAssignment_NumBestGenes` | int64 | Genes tied for the best call: 1 for `unique_gene` and `best_of_multiple`, 2+ for `tied`, 0 otherwise |
| 18 | `GeneAssignment_Candidates` | str | Every touched same-strand gene in rank order, as `NAME (ID) matched/introns exonic_frac`, separated by `; ` |
| 19 | `GeneAssignment_AntisenseGenes` | str | Opposite-strand genes touched, as `NAME (ID)`, separated by `; ` |
| 20 | `GeneAssignment_NearbyGeneDist` | int64 | For `nearby_gene`: bp from the nearest aligned block to the gene |
| 21 | `GeneAssignment_NearbyGeneSide` | str | For `nearby_gene`: `downstream_of_3p` (past the gene's 3′ end) or `upstream_of_5p` |

Empty means not applicable (null in Parquet).

### Examples (chr21 test data)

| | DONSON read (…/1001086) | IFNAR2 read (…/1224638) |
|---|---|---|
| `GeneAssignment_Status` | `best_of_multiple` | `tied` |
| `GeneAssignment_SupportLevel` | `known_junctions` | `known_junctions` |
| `GeneAssignment_GeneName` | DONSON | IFNAR2 |
| `GeneAssignment_IntronsMatched` / `ReadIntrons` | 8 / 8 | 8 / 8 |
| `GeneAssignment_NumSameStrandGenes` / `NumBestGenes` | 2 / 1 | 2 / 2 |
| `GeneAssignment_Candidates` | `DONSON (ENSG00000159147.20) 8/8 1.0; ENSG00000249209 (ENSG00000249209.2) 1/8 0.138` | `IFNAR2 (ENSG00000159110.24) 8/8 1.0; IFNAR2-IL10RB (ENSG00000249624.10) 8/8 1.0` |

The IFNAR2 read fits IFNAR2 and the IFNAR2-IL10RB readthrough gene equally well. With `--skip-readthrough` it would be `unique_gene`.

### Joining the labels back, and filtering

The table holds only each alignment's identity and its labels.
- **To a compare table:** merge on `query_name`, since it has one row per read.
- **To a PAF:** use `query_name`, `target_name`, `target_start`, `target_end` and `strand` together.

```python
import pandas as pd
cmp = pd.read_parquet("sample.compare.parquet")
lab = pd.read_parquet("sample.gene_assignment.parquet", dtype_backend="numpy_nullable")
cmp = cmp.merge(lab.drop(columns=["tp_tag", "target_name", "target_start", "target_end", "strand"]),
                on="query_name", how="left", validate="one_to_one")

# one clearly assigned gene, supported by exons
confident = lab[lab["GeneAssignment_Status"].isin(["unique_gene", "best_of_multiple"])
                & lab["GeneAssignment_SupportLevel"].isin(["known_junctions", "exonic"])]
```

`dtype_backend="numpy_nullable"` keeps the integer columns as integers where they have nulls.

Gene-level read counts that don't filter on `GeneAssignment_SupportLevel` include `intronic_only` reads, which inflates genes with long introns. Intronic reads are also where artifacts such as internal priming and genomic DNA collect.

---

## Options

| Option | Default | Description |
|---|---|---|
| `-a, --alignments` | required | PAF, TSV or Parquet |
| `-g, --gff` | required | GENCODE-style GFF3 |
| `-o, --output` | required | The classification table |
| `--summary` | none | Also write the status × support-level counts as a TSV |
| `--format` | `auto` | `auto`, `paf`, `tsv` or `parquet` |
| `--min-exonic-frac` | 0.2 | Exonic fraction that makes a gene's support `exonic` |
| `--tie-margin` | 0.1 | Lead in exonic fraction needed to beat a gene with the same matched introns |
| `--near-max-dist` | 5000 | Max bp from an aligned block to a same-strand gene for `nearby_gene` |
| `--skip-readthrough` | off | Leave out genes tagged `readthrough_gene` |
| `--gene-types` | all | Gene types to use, comma-separated |
| `--id-col`, `--tp-col`, `--chrom-col`, `--strand-col`, `--start-col`, `--end-col`, `--cs-col`, `--cigar-col` | see Tables | Table columns |

---

## Performance

- **The annotation** is read into memory once: about 79k genes, 460k merged exons and 712k annotated introns for GENCODE v50, built in about 2 s.
- **The alignments** are streamed. A typical read touches 1–3 genes, so each one needs a binary search for its genes, then walks over those genes' sorted exon and intron lists. Memory doesn't grow with the number of alignments.

On an Apple M5, with HPRC HG00099 Kinnex transcript clusters and the slim GENCODE v50 GFF:
- **Whole compare table:** 907,941 rows, Parquet in and out, in 3.6 s with 355 MB peak memory (5.6 s with the full GFF). On battery power, expect about twice as long.
- **PAF:** the 918,091-line PAF took the same time as the compare table when run side by side.
- **Output size:** the table is about 34 MB as Parquet.
