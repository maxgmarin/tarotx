# tarotx

***Every read has a fortune.***

A toolkit for stranded, trimmed long-read RNA-seq alignments. `tarotx` reads each alignment's fortune: its subcommands label every alignment with information for filtering and quality control. The *tx* is for transcripts.

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

## Installation

Building needs the [Rust toolchain](https://www.rust-lang.org/tools/install) and a C compiler, because zstd (for Parquet) is compiled from source. On macOS, the Xcode command line tools provide one (`xcode-select --install`).

```bash
cargo install --locked --path .
```

`--locked` builds with the exact dependency versions in `Cargo.lock`.

### Static Linux build

Cross-compile from macOS with the musl toolchain:
- `brew install FiloSottile/musl-cross/musl-cross`
- `rustup target add x86_64-unknown-linux-musl`

The linker is set in `.cargo/config.toml`. The `CC`/`AR` variables are needed for zstd's C code:

```bash
CC_x86_64_unknown_linux_musl=x86_64-linux-musl-gcc \
AR_x86_64_unknown_linux_musl=x86_64-linux-musl-ar \
cargo build --release --target x86_64-unknown-linux-musl
```

---

## Quick start

```bash
# Gene assignment of a minimap2 PAF (cs or cg tags), with the counts as a TSV
tarotx label-by-gene -a sample.paf.gz -g gencode.v50.primary_assembly.genes_exons.slim.gff3.gz \
  -o sample.gene_assignment.parquet --summary sample.gene_assignment.summary.tsv

# Internal priming of a PAF: the label table, and the PAF with the labels as tags
tarotx label-internal-priming -a sample.paf.gz -p GRCh38.polyA_tracts_W10_P80.bed.gz \
  -g gencode.v50.annotation.gff3.gz -o sample.int_priming.tsv.gz --tagged-out sample.int_priming.paf.gz

# A BAM: convert it to PAF first, keeping the unmapped reads (see "SAM/BAM input")
maligno sam2paf -U sample.bam > sample.paf

# Both on a maligno compare table (the default column names are maligno's)
tarotx label-by-gene -a sample.compare.parquet -g genes.gff3.gz -o sample.gene_assignment.parquet
tarotx label-internal-priming -a sample.compare.parquet -p tracts.bed.gz -g gencode.gff3.gz -o sample.int_priming.parquet
```

Each subcommand writes its own table, and both tables start with the same six columns: `query_name`, `tp_tag`, `target_name`, `target_start`, `target_end` and `strand`. To combine the tables, join on those columns, for example in pandas or polars. Both subcommands take `--summary FILE` to save their counts as a TSV.

**Inputs by subcommand:**

| Input | `label-by-gene` | `label-internal-priming` |
|---|---|---|
| minimap2 PAF | Yes (needs `cs:Z` or `cg:Z`) | Yes |
| TSV / Parquet table | Yes (needs a cs or CIGAR column) | Yes |
| SAM/BAM | No: convert to PAF with `maligno sam2paf -U` | No: the same |

**Provenance:** every Parquet table records `tarotx_version` and `tarotx_command` (the subcommand) in its metadata.

---

## One row per input alignment

Both subcommands label **every input alignment, in input order, and write one row for each**: every PAF line, and every row of a TSV or Parquet table. They never choose among a read's alignments. The input decides what is labelled, and the output says which alignment each row is (`query_name`, `tp_tag`, the coordinates and strand). What to keep is decided downstream.

How that plays out:

| Case | What you get |
|---|---|
| **One alignment per read** | One row per read |
| **maligno compare table** | One row per read: maligno has already reduced each read to one representative alignment, and that alignment is the one labelled |
| **Split / chimeric read in a PAF** (a primary plus supplementary pieces covering different parts of the read) | One row per piece. minimap2 writes supplementary pieces as `tp:A:P`, the same as the primary, so `tp_tag` can't tell them apart. In `label-internal-priming`, `read_3p_clip` shows which piece holds the read's 3′ end. |
| **Secondary alignments** (other placements of the same part of the read) | One row each, with `tp_tag` `S`, when the aligner wrote them (minimap2's `--secondary=no` leaves them out) |
| **Unmapped read** | One row with empty coordinates and labels, from a PAF line with strand `*` (written by `maligno sam2paf -U`, or by minimap2 with `--paf-no-hit`) or from a table row without strand |
| **A BAM** | Not read. Convert it with `maligno sam2paf -U` (below); each record is then a PAF line, as above |

**Counts are of alignments, not reads.** The summaries (stderr and `--summary`) count rows, so a split read counts once per piece. For a compare table, alignments and reads are the same thing. In the HG00099 test sample, 9,708 of 907,941 reads (1.1%) have more than one alignment in the PAF.

To keep one row per read from a PAF, filter the table afterwards, e.g. on `tp_tag`, `read_3p_clip` or the coordinates, or label a table that already has one alignment per read.

### SAM/BAM input

Neither subcommand reads SAM/BAM. Convert it to PAF with maligno's `sam2paf`, with `-U` so that unmapped reads are kept:

```bash
maligno sam2paf -U sample.bam > sample.paf
```

- **Unmapped reads:** with `-U`, each one becomes a placeholder line with strand `*` that keeps the read length, so it gets a row. Without `-U`, unmapped reads are dropped and get no row.
- **Which records:** by default every record is converted: primary, supplementary and secondary. Add `-p` to skip secondary records, or `-P` to keep only the primary ones.
- **Supplementary records** are written as `tp:A:P`, like the primary, as minimap2 does in PAF.
- **`cs` tags** come from the BAM's `cs` tag, or are rebuilt from its `MD` tag. `label-by-gene` needs them, for the introns.
- **Row order** follows the BAM.

On the HG00099 test sample, the BAM converted this way gives exactly the same rows and labels as minimap2's own PAF, in both subcommands: 918,091 rows, the 1,217 unmapped reads included.

---

## Code layout

Each subcommand is one folder under `src/`. Its `mod.rs` holds the command-line options (`Args`) and the entry point (`run`), and starts with a list of the folder's files.

```
src/
  main.rs                       the command line: the subcommand list, and error printing
  common/                       shared by every subcommand
    io.rs                         text input and output, with transparent gzip/BGZF
    strand.rs                     Strand (+/-)
    mod.rs                        Parquet provenance metadata
  label_by_gene/                tarotx label-by-gene
    annot.rs  blocks.rs  label.rs  run.rs  parquet_in.rs  table.rs
  label_internal_priming/       tarotx label-internal-priming
    points.rs  label.rs  run.rs  parquet_in.rs  table.rs
tests/
  label_by_gene.rs              end-to-end tests, one file per subcommand
  label_internal_priming.rs
  data/<subcommand>/            their hand-built fixtures
```

The two subcommands share only code that was already identical in both. Each keeps its own table writer and runners, so changing one can't change the other's output.

**To add a subcommand:**
1. Add a folder `src/<name>/` whose `mod.rs` has `pub struct Args` (`#[derive(clap::Args)]`) and `pub fn run(args: Args) -> anyhow::Result<()>`.
2. Add a variant to `Command` in `src/main.rs`, with its name and its dispatch.
3. Add `tests/<name>.rs` with fixtures in `tests/data/<name>/`, and `docs/<name>.md`.

`cargo test --release` runs the unit tests in each module and the end-to-end tests in `tests/`.

---

## License

GNU General Public License v3.0; see [LICENSE](LICENSE).
