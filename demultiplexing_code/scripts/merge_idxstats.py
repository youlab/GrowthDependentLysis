#!/usr/bin/env python3
"""
merge_idxstats.py — Combine SAMtools idxstats outputs into a
sample × amplicon read count matrix.

Galaxy equivalent:
    After running SAMtools idxstats on each BAM, use
    "Join, Subtract and Group" → "Join two files" iteratively,
    or upload all idxstats to a Collection and use
    "Collapse Collection" → then this script's logic.

idxstats columns: contig  length  mapped_reads  unmapped_reads

Output count_matrix.tsv:
    rows    = amplicons (Seq_ID)
    columns = samples
    values  = mapped reads

Usage:
    python3 merge_idxstats.py \
        --idxstats-dir results/06_counts \
        --primers-csv  validated_primers_k20-seed18_with_adaptors.csv \
        --out-dir      results/06_counts
"""

import argparse
import os
from pathlib import Path

import pandas as pd


def merge(idxstats_dir: str, primers_csv: str, out_dir: str):
    primers = pd.read_csv(primers_csv)
    # Normalise contig names to match BAM reference header
    primers["contig"] = (
        primers["Seq_ID"]
        .str.replace(" ", "_")
        .str.replace("/", "-")
    )

    frames = {}
    for f in sorted(Path(idxstats_dir).glob("*.idxstats")):
        # Filename pattern: {sample}_{amplicon_id}.idxstats
        # We want the sample part (everything before the last underscore block)
        sample = f.stem  # e.g. "1-1_ES55"
        df = pd.read_csv(
            f, sep="\t",
            names=["contig", "length", "mapped", "unmapped"],
        )
        # Drop the '*' unmapped row
        df = df[df["contig"] != "*"]
        frames[sample] = df.set_index("contig")["mapped"]

    if not frames:
        print("No .idxstats files found.")
        return

    # Build matrix: index = contig, columns = samples
    matrix = pd.DataFrame(frames).fillna(0).astype(int)

    # Add human-readable Seq_ID as index
    contig_to_seqid = dict(zip(primers["contig"], primers["Seq_ID"]))
    matrix.index = [contig_to_seqid.get(c, c) for c in matrix.index]
    matrix.index.name = "Seq_ID"

    # Sort by Seq_ID
    matrix = matrix.sort_index()

    # ── Outputs ───────────────────────────────────────────────────────────────

    # 1. Full count matrix (sample × amplicon)
    matrix_path = os.path.join(out_dir, "count_matrix.tsv")
    matrix.to_csv(matrix_path, sep="\t")
    print(f"Count matrix ({matrix.shape[0]} amplicons × "
          f"{matrix.shape[1]} samples) → {matrix_path}")

    # 2. Long format (easier to import into R/pandas/Galaxy)
    long = matrix.reset_index().melt(
        id_vars="Seq_ID", var_name="Sample", value_name="Mapped_Reads"
    )
    long_path = os.path.join(out_dir, "count_matrix_long.tsv")
    long.to_csv(long_path, sep="\t", index=False)
    print(f"Long format → {long_path}")

    # 3. Summary stats per amplicon
    summary = matrix.agg(["sum", "mean", "median", "min", "max"]).T
    summary.columns = ["Total_Reads", "Mean", "Median", "Min", "Max"]
    summary.index.name = "Seq_ID"

    # Merge primer metadata
    meta = primers[["Seq_ID", "Product_Size", "Fwd_Tm", "Rev_Tm",
                     "Forward_Primer", "Reverse_Primer"]].set_index("Seq_ID")
    summary = summary.join(meta, how="left")

    summary_path = os.path.join(out_dir, "amplicon_summary.tsv")
    summary.to_csv(summary_path, sep="\t")
    print(f"Amplicon summary → {summary_path}")

    # 4. Per-sample total mapped reads
    sample_totals = matrix.sum(axis=0).reset_index()
    sample_totals.columns = ["Sample", "Total_Mapped_Reads"]
    sample_totals_path = os.path.join(out_dir, "sample_totals.tsv")
    sample_totals.to_csv(sample_totals_path, sep="\t", index=False)
    print(f"Sample totals → {sample_totals_path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--idxstats-dir", required=True)
    p.add_argument("--primers-csv",  required=True)
    p.add_argument("--out-dir",      required=True)
    args = p.parse_args()
    merge(args.idxstats_dir, args.primers_csv, args.out_dir)
