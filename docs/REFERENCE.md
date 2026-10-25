# tarotx — Reference

The reference for **tarotx**: what every subcommand shares, how alignments become rows, building from source, and how the code is laid out. For installation and a quick start, see the [README](../README.md). Each subcommand's rules, label columns and options are in its own document: [`label-by-gene`](label-by-gene.md) and [`label-internal-priming`](label-internal-priming.md).

**Contents**

- [What every subcommand shares](#what-every-subcommand-shares): the input and output options, the 12 leading columns, the summary, the annotation
- [One row per input alignment](#one-row-per-input-alignment): split reads, secondary alignments, unmapped reads
  - [SAM/BAM input](#sambam-input): converting with `maligno sam2paf -U`
- [Build](#build): building from source, and the static Linux build
- [Code layout](#code-layout): the source tree, and how to add a subcommand

---

## What every subcommand shares

**The same options for the input and output:**

| Option | Default | Description |
|---|---|---|
| `-a, --alignments` | required | A minimap2 PAF (`.paf`, `.paf.gz`), or a TSV (`.tsv`, `.txt`, optionally `.gz`) or Parquet table with a header line. `-a -` reads stdin with `--in-format` (not for Parquet) |
| `--in-format` | `auto` | `auto` (from the file name), `paf`, `tsv` or `parquet` |
| `-o, --output` | required | The per-alignment table: `.parquet` (zstd, row groups of 65,536 rows), `.tsv.gz` (any `.gz`), `-` (TSV on stdout), or plain TSV |
| `--summary` | none | Also write the counts printed at the end as a TSV |
| `--id-col`, `--chrom-col`, `--strand-col`, `--start-col`, `--end-col` | `query_name`, `target_name_A`, `strand_A`, `target_start_A`, `target_end_A` | Table input: the read ID, chromosome, strand (`+`/`-`; anything else is unmapped), 0-based start and exclusive end. Required |
| `--query-len-col`, `--query-start-col`, `--query-end-col`, `--mapq-col`, `--tp-col` | `query_len`, `query_start_A`, `query_end_A`, `mapq_A`, `tp_tag_A` | Table input: the read length, the query span (on the read as sequenced, as in PAF), mapq and tp tag. Optional; the matching leading columns are empty without them |

The defaults are the column names of a maligno compare table, so it needs no column options.

**The same 12 leading columns.** Every table starts with them, so the tables of different subcommands can be joined on them, for example in pandas or polars:

| # | Column | Type | From a PAF | From a table |
|---|---|---|---|---|
| 1 | `query_name` | str | column 1 | `--id-col` |
| 2–4 | `target_name`, `target_start`, `target_end` | str, uint64, uint64 | columns 6, 8, 9 | `--chrom-col`, `--start-col`, `--end-col` |
| 5 | `strand` | str | column 5 | `--strand-col` |
| 6 | `query_len` | uint64 | column 2 | `--query-len-col` |
| 7–8 | `query_start`, `query_end` | uint64 | columns 3, 4 | `--query-start-col`, `--query-end-col` |
| 9 | `mapq` | uint64 | column 12 | `--mapq-col` |
| 10 | `tp_tag` | str | `tp:A` tag (`P`, `S`, `I`, `i`) | `--tp-col` |
| 11 | `aln_type` | str | `primary`, `secondary` or `inversion`, from the tp tag; `unmapped` for an unmapped alignment | the same |
| 12 | `read_3p_clip` | uint64 | bases of the read past the alignment's 3′ end: `query_len - query_end` | the same, when both are known |

- **Unmapped alignments:** the coordinates, strand, query span and mapq are empty, and `aln_type` is `unmapped`.
- **Missing values:** a value the input doesn't have is empty in TSV, and null in Parquet. For example, a table without a mapq column leaves `mapq` empty.
- **`read_3p_clip`:** with trimmed reads, there's no tail or adapter left to clip, so a large value points to a split read or an unaligned end. For a split read, the piece with a small value holds the read's real 3′ end.

**The same summary shape.** On stderr: the loading time, then `N alignments (U unmapped), …`, the subcommand's counts, and the total time. `--summary FILE` writes the counts as a TSV.

**Provenance:** every Parquet table records `tarotx_version` and `tarotx_command` (the subcommand) in its metadata.

**Annotation:** use the full GENCODE GFF3 (e.g. `gencode.v50.primary_assembly.annotation.gff3.gz`) for both subcommands.
- `label-by-gene` reads its `gene` and `exon` lines, and `label-internal-priming` its `transcript` lines.
- A GFF3 without the lines a subcommand needs is an error. A gene/exon-only GFF3, for instance, is refused by `label-internal-priming`.

**Inputs by subcommand:**

| Input | `label-by-gene` | `label-internal-priming` |
|---|---|---|
| minimap2 PAF | Yes (needs `cs:Z` or `cg:Z`) | Yes |
| TSV / Parquet table | Yes (needs a cs or CIGAR column: `--cs-col`, default `cs_A`, or `--cigar-col`) | Yes |
| SAM/BAM | No: convert to PAF with `maligno sam2paf -U` | No: the same |

---

## One row per input alignment

Both subcommands label **every input alignment, in input order, and write one row for each**: every PAF line, and every row of a TSV or Parquet table. They never choose among a read's alignments. The input decides what is labelled, and the 12 leading columns say which alignment each row is. Nothing is filtered: what to keep is decided downstream.

How that plays out:

| Case | What you get |
|---|---|
| **One alignment per read** | One row per read |
| **maligno compare table** | One row per read: maligno has already reduced each read to one representative alignment, and that alignment is the one labelled |
| **Split / chimeric read in a PAF** (a primary plus supplementary pieces covering different parts of the read) | One row per piece. minimap2 writes supplementary pieces as `tp:A:P`, the same as the primary, so `tp_tag` can't tell them apart. `read_3p_clip` shows which piece holds the read's 3′ end. |
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

## Build

Building needs the [Rust toolchain](https://www.rust-lang.org/tools/install) and a C compiler, because zstd (for Parquet) is compiled from source. On macOS, the Xcode command line tools provide one (`xcode-select --install`).

```bash
cargo build --release              # binary at target/release/tarotx
cargo install --locked --path .    # or install to ~/.cargo/bin (usually on PATH)
```

`--locked` builds with the exact dependency versions in `Cargo.lock`.

### Static Linux build

A static Linux binary (no runtime dependencies, e.g. for an HPC cluster) is cross-compiled from macOS with the musl toolchain:
- `brew install FiloSottile/musl-cross/musl-cross`
- `rustup target add x86_64-unknown-linux-musl`

The linker is set in `.cargo/config.toml`. The `CC`/`AR` variables are needed for zstd's C code:

```bash
CC_x86_64_unknown_linux_musl=x86_64-linux-musl-gcc \
AR_x86_64_unknown_linux_musl=x86_64-linux-musl-ar \
cargo build --release --target x86_64-unknown-linux-musl
```

The binary is at `target/x86_64-unknown-linux-musl/release/tarotx`, and `file` reports it as `static-pie linked`.

---

## Code layout

Each subcommand is one folder under `src/`. Its `mod.rs` holds the command-line options (`Args`) and the entry point (`run`), and starts with a list of the folder's files. Everything about the input and the output table is shared, in `src/common/`, so a subcommand holds only its labelling logic.

```
src/
  main.rs                       the command line: the subcommand list, and error printing
  common/                       shared by every subcommand
    input.rs                      reading PAF, TSV and Parquet input; the 12 leading columns
    args.rs                       the shared input, output and table-column options
    table.rs                      the output table (TSV or Parquet)
    io.rs                         text input and output, with transparent gzip/BGZF
    strand.rs                     Strand (+/-)
    mod.rs                        Parquet provenance metadata
  label_by_gene/                tarotx label-by-gene
    annot.rs  blocks.rs  label.rs  run.rs
  label_internal_priming/       tarotx label-internal-priming
    points.rs  label.rs  run.rs
tests/
  label_by_gene.rs              end-to-end tests, one file per subcommand
  label_internal_priming.rs
  data/<subcommand>/            their hand-built fixtures
```

**To add a subcommand:**
1. Add a folder `src/<name>/` whose `mod.rs` has `pub struct Args` (`#[derive(clap::Args)]`, flattening `common::args::InputArgs` and `TableColumnArgs`, plus its own `-o/--output` and `--summary`) and `pub fn run(args: Args) -> anyhow::Result<()>`.
2. In `run`, open the table with `common::table::Sink::open(output, &common::input::table_columns(label_columns), name)`, then call `common::input::for_each_alignment` and write `aln.lead_cells()` plus the labels for each alignment.
3. Add a variant to `Command` in `src/main.rs`, with its name and its dispatch.
4. Add `tests/<name>.rs` with fixtures in `tests/data/<name>/`, and `docs/<name>.md`.

`cargo test --release` runs the unit tests in each module and the end-to-end tests in `tests/`.
