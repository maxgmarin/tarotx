#!/usr/bin/env python3
"""Append the key tarotx labels to a maligno alignment comparison table.

Inputs:
  --compare           the maligno compare table (Parquet), one row per read
  --gene-labels       `tarotx label-by-gene` output for that table (Parquet)
  --internal-priming  `tarotx label-internal-priming` output for that table (Parquet), run with -g and the default
                      thresholds (--polya-max-dist 3, --tes-max-dist 50)
Output (-o): the compare table with every original column, in the original row order, plus these columns at the end,
under tarotx's own names:
  GeneAssignment_Status, GeneAssignment_SupportLevel, GeneAssignment_GeneID, GeneAssignment_GeneName,
  GeneAssignment_GeneType, GeneAssignment_NumSameStrandGenes, GeneAssignment_NumBestGenes,
  InternalPriming_Flag, InternalPriming_PolyAWithin3bp, InternalPriming_TxEndWithin50bp,
  InternalPriming_PolyAEdgeDist, InternalPriming_TxEndSignedDist

tarotx labels the compare table's A side (its default columns are the _A ones), so these labels describe the A-side
alignment. Rows are joined on query_name, after strict checks (see check_labels); unmapped reads keep empty labels.

Usage: python annotate_compare_table.py --compare C.parquet --gene-labels G.parquet --internal-priming I.parquet -o OUT.parquet
Needs pandas and pyarrow.
"""
import argparse
import sys

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

GENE_COLUMNS = [
    "GeneAssignment_Status", "GeneAssignment_SupportLevel", "GeneAssignment_GeneID", "GeneAssignment_GeneName",
    "GeneAssignment_GeneType", "GeneAssignment_NumSameStrandGenes", "GeneAssignment_NumBestGenes",
]
INTERNAL_PRIMING_COLUMNS = [
    "InternalPriming_Flag", "InternalPriming_PolyAWithin3bp", "InternalPriming_TxEndWithin50bp",
    "InternalPriming_PolyAEdgeDist", "InternalPriming_TxEndSignedDist",
]
# tarotx's leading columns, and the compare-table columns (A side) they were read from
COORDS = {"target_name": "target_name_A", "target_start": "target_start_A", "target_end": "target_end_A",
          "strand": "strand_A"}


def read(path, columns=None):
    """A Parquet table with Arrow-backed types, so integer and boolean columns with nulls keep their types."""
    return pd.read_parquet(path, columns=columns, dtype_backend="pyarrow")


def check_labels(name, labels, compare, columns):
    """Stop unless `labels` (a tarotx table) holds exactly the compare table's reads, once each, with the same A-side
    alignment, and has the columns to transfer."""
    missing = [c for c in columns if c not in labels.columns]
    if missing:
        sys.exit(f"{name}: no column {missing}. label-internal-priming needs -g and the default --polya-max-dist 3 "
                 "and --tes-max-dist 50" if "InternalPriming" in missing[0] else f"{name}: no column {missing}")
    clash = [c for c in columns if c in compare.columns]
    if clash:
        sys.exit(f"the compare table already has column(s) {clash}")
    if labels["query_name"].duplicated().any():
        sys.exit(f"{name}: query_name is not unique ({labels['query_name'].duplicated().sum()} repeats)")
    ours, theirs = set(labels["query_name"]), set(compare["query_name"])
    if ours != theirs:
        sys.exit(f"{name}: {len(theirs - ours)} compare-table reads are missing, {len(ours - theirs)} reads are extra")
    # the same alignment: mapped in both or neither, and at the same place when mapped
    both = compare[["query_name", *COORDS.values()]].merge(labels[["query_name", *COORDS]], on="query_name",
                                                           validate="one_to_one")
    mapped_compare = both["strand_A"].isin(["+", "-"])
    mapped_labels = both["target_name"].notna()
    if (mapped_compare != mapped_labels).any():
        sys.exit(f"{name}: {(mapped_compare != mapped_labels).sum()} reads are mapped in one table and not the other")
    m = both[mapped_compare]
    for ours_col, theirs_col in COORDS.items():
        differ = (m[ours_col].astype(str) != m[theirs_col].astype(str)).sum()
        if differ:
            sys.exit(f"{name}: {ours_col} differs from the compare table's {theirs_col} on {differ} reads")


def main():
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--compare", required=True, help="maligno compare table (Parquet)")
    a.add_argument("--gene-labels", required=True, help="tarotx label-by-gene output (Parquet)")
    a.add_argument("--internal-priming", required=True, help="tarotx label-internal-priming output (Parquet)")
    a.add_argument("-o", "--output", required=True, help="annotated compare table (Parquet)")
    args = a.parse_args()

    compare = read(args.compare)
    if compare["query_name"].duplicated().any():
        sys.exit("the compare table's query_name is not unique")
    keys = ["query_name", *COORDS]
    gene = read(args.gene_labels, keys + GENE_COLUMNS)
    ip = read(args.internal_priming, keys + INTERNAL_PRIMING_COLUMNS)
    check_labels(args.gene_labels, gene, compare, GENE_COLUMNS)
    check_labels(args.internal_priming, ip, compare, INTERNAL_PRIMING_COLUMNS)

    # left joins on query_name keep the compare table's rows and their order
    out = (compare
           .merge(gene[["query_name", *GENE_COLUMNS]], on="query_name", how="left", validate="one_to_one")
           .merge(ip[["query_name", *INTERNAL_PRIMING_COLUMNS]], on="query_name", how="left", validate="one_to_one"))
    # keep the compare table's schema exactly (pandas would drop its not-null flags), then the new columns
    table = pa.Table.from_pandas(out, preserve_index=False)
    added_fields = [table.schema.field(c) for c in GENE_COLUMNS + INTERNAL_PRIMING_COLUMNS]
    schema = pa.schema(list(pq.read_schema(args.compare)) + added_fields)
    pq.write_table(table.cast(schema), args.output, compression="zstd")
    added = len(GENE_COLUMNS) + len(INTERNAL_PRIMING_COLUMNS)
    print(f"{args.output}: {len(out)} rows, {len(compare.columns)} compare-table columns + {added} label columns",
          file=sys.stderr)


if __name__ == "__main__":
    main()
