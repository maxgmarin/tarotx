# tarotx

`tarotx` is a toolkit for QC and processing of long-read RNA-seq alignments: its subcommands label every alignment with information for filtering and quality control. The *tx* is for transcripts.

`tarotx` is one program with a subcommand for each job. The subcommands label alignments of **stranded long-read RNA-seq reads whose adapters and poly(A) tails have been trimmed**, and write one row per labelled alignment, as TSV or Parquet. PacBio Iso-Seq/Kinnex reads after `lima` and `isoseq refine` are reads of this kind. See [Input reads](#input-reads-what-tarotx-assumes) before using other data.

| Subcommand | What it does | Documentation |
|---|---|---|
| `label-by-gene` | Assigns each alignment to the annotated gene it most likely comes from. It reports how clear that choice is and how the gene is supported (known junctions, exonic or intronic only). | [docs/label-by-gene.md](docs/label-by-gene.md) |
| `label-internal-priming` | Labels each alignment's 3′ end by its distance to the nearest genomic poly(A) tract and annotated transcript end, and flags likely internal priming. | [docs/label-internal-priming.md](docs/label-internal-priming.md) |

`tarotx <subcommand> --help` lists a subcommand's options.

The poly(A) tracts that `label-internal-priming` needs come from [polyscan](https://github.com/maxgmarin/polyscan), a separate tool.

---

## Input reads: what tarotx assumes

tarotx is currently designed for one kind of data, and its rules are built on two assumptions about the reads:

| Assumption | What relies on it |
|---|---|
| **Stranded:** each read is oriented 5′→3′ along its transcript, so the alignment strand is the transcript strand. | **Both subcommands.** `label-by-gene` compares the alignment only with genes on its own strand, and calls `antisense_only` from the other strand. `label-internal-priming` takes the 3′ end as the alignment's end on `+` and its start on `−`, and looks only at same-strand poly(A) tracts and transcript ends. |
| **Trimmed:** the adapters/primers and the poly(A) tail have been removed, so the read ends where the transcript sequence ends. | **Mainly `label-internal-priming`.** It treats the alignment's 3′ end as the site where the oligo(dT) primer bound. If the tail or an adapter is left on, it may be soft-clipped (`read_3p_clip` grows) or align into a downstream genomic A-run, and either way the 3′ end no longer marks where the transcript ended. `label-by-gene` works from aligned blocks and introns, so it depends mainly on strandedness. |

**Neither assumption is checked.** tarotx can't tell an unstranded or untrimmed read from a stranded, trimmed one: it labels them all the same way, and their labels can be wrong. Examples of input that isn't supported yet:
- unstranded or unoriented reads, such as ONT cDNA reads before orientation, or any data where a read can come from either strand;
- reads that still carry their adapters or poly(A) tail.

Support for unstranded data and untrimmed reads may be added in a later version.

---

## Install

### Installing prebuilt binaries

Prebuilt binaries for Linux and macOS are attached to each [GitHub release](https://github.com/maxgmarin/tarotx/releases). Each tarball holds the `tarotx` binary and this README.

```bash
# Go to the directory where you want to place the tarotx binary
cd <directory_on_PATH>

# Download the latest release (Linux x86_64; static, no runtime dependencies)
wget https://github.com/maxgmarin/tarotx/releases/download/v0.2.0/tarotx-v0.2.0-x86_64-unknown-linux-musl.tar.gz
# or macOS (Apple Silicon):
# wget https://github.com/maxgmarin/tarotx/releases/download/v0.2.0/tarotx-v0.2.0-aarch64-apple-darwin.tar.gz

# Uncompress the .tar.gz to get the binary
tar -xzf tarotx-v0.2.0-x86_64-unknown-linux-musl.tar.gz

# Test that the binary runs
./tarotx --help
./tarotx --version

# Then add the directory with the binary to your $PATH variable
```

`checksums-sha256.txt` on the release page verifies the download (`shasum -a 256 -c checksums-sha256.txt`).

### Building from source

Requires a [Rust toolchain](https://rustup.rs) and a C compiler.

```bash
cargo build --release              # binary at target/release/tarotx
cargo install --locked --path .    # or install to ~/.cargo/bin (usually on PATH)
```

A static Linux (musl) build, e.g. for an HPC cluster, is described in the [reference](docs/REFERENCE.md#build).

---

## Quick start

```bash
# Gene assignment of a minimap2 PAF (cs or cg tags), with the counts as a TSV
tarotx label-by-gene -a sample.paf.gz -g gencode.v50.primary_assembly.annotation.gff3.gz \
  -o sample.gene_assignment.parquet --summary sample.gene_assignment.summary.tsv

# Internal priming of the same PAF
tarotx label-internal-priming -a sample.paf.gz -p GRCh38.polyA_tracts_W10_P80.bed.gz \
  -g gencode.v50.primary_assembly.annotation.gff3.gz -o sample.int_priming.parquet --summary sample.int_priming.summary.tsv

# A BAM: convert it to PAF first, keeping the unmapped reads (see docs/REFERENCE.md#sambam-input)
maligno sam2paf -U sample.bam > sample.paf

# Both on a maligno compare table (the default column names are maligno's)
tarotx label-by-gene -a sample.compare.parquet -g gencode.gff3.gz -o sample.gene_assignment.parquet
tarotx label-internal-priming -a sample.compare.parquet -p tracts.bed.gz -g gencode.gff3.gz -o sample.int_priming.parquet
```

---

## Extended documentation

- **[`docs/REFERENCE.md`](docs/REFERENCE.md):**
  - what every subcommand shares: the input and output options, the 12 leading columns, the summary and the annotation;
  - the one-row-per-input-alignment rule, with split reads, secondary alignments and unmapped reads;
  - SAM/BAM input;
  - the static Linux build;
  - the code layout.
- **[`docs/label-by-gene.md`](docs/label-by-gene.md):** how alignments are assigned to genes, and the `GeneAssignment_*` columns.
- **[`docs/label-internal-priming.md`](docs/label-internal-priming.md):** how internal priming is labelled, and the `InternalPriming_*` columns.

---

## Scripts

`scripts/annotate_alncompare_table.WithTarotxInfo.py` appends the key tarotx labels to a maligno alignment comparison table. It takes the compare table plus the `label-by-gene` and `label-internal-priming` tables made from it, and writes the compare table with 12 label columns added at the end:
- 7 `GeneAssignment_*` columns: Status, SupportLevel, GeneID, GeneName, GeneType, NumSameStrandGenes, NumBestGenes;
- 5 `InternalPriming_*` columns: Flag, PolyAWithin3bp, TxEndWithin50bp, PolyAEdgeDist, TxEndSignedDist.

Rows are joined on `query_name`. The script stops unless every read is present once in all three tables, at the same A-side alignment. It needs pandas and pyarrow, and expects `label-internal-priming` to have been run with `-g` and the default thresholds.

```bash
python scripts/annotate_alncompare_table.WithTarotxInfo.py --compare sample.compare.parquet \
  --gene-labels sample.gene_labels.parquet --internal-priming sample.internal_priming.parquet \
  -o sample.compare.annotated.parquet
```

---

## License

GNU General Public License v3.0; see [LICENSE](LICENSE).
