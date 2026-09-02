#!/usr/bin/env python3
"""
primers_to_cutadapt_barcodes.py — Build a FASTA of target barcodes for
cutadapt --demultiplex mode.

Galaxy equivalent:
    Cutadapt → "5' adapter" → Source: File
    The uploaded FASTA is exactly what this script produces.

Each entry is a linked adapter:
    >Seq_ID
    BARCODE...FORWARD_PRIMER

This tells cutadapt:
  1. The read must start with BARCODE (within --error-rate)
  2. Immediately followed by the forward primer sequence
  → Only true on-target reads pass; chimeras and off-targets are discarded

Usage:
    python3 primers_to_cutadapt_barcodes.py \
        --primers-csv validated_primers_k20-seed18_with_adaptors.csv \
        --out-fasta   results/04_demux_barcode/barcodes.fasta \
        --mismatches  1
"""

import argparse
import pandas as pd


def build_fasta(primers_csv: str, out_fasta: str, mismatches: int):
    df = pd.read_csv(primers_csv)

    with open(out_fasta, "w") as f:
        for _, row in df.iterrows():
            seq_id  = row["Seq_ID"].replace(" ", "_").replace("/", "-")
            barcode = row["Target_Barcode"].strip().upper()
            fwd_primer = row["Forward_Primer"].strip().upper()

            # Linked adapter: barcode anchored at 5' end, followed by primer
            # Equivalent to cutadapt: -g ^BARCODE...PRIMER
            f.write(f">{seq_id}\n")
            f.write(f"^{barcode}...{fwd_primer}\n")

    bc_len = len(df["Target_Barcode"].iloc[0])
    error_rate = round(mismatches / bc_len, 3)
    print(f"Wrote {len(df)} barcode entries → {out_fasta}")
    print(f"Suggested cutadapt --error-rate: {error_rate}  "
          f"(={mismatches} mismatch / {bc_len}bp barcode)")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--primers-csv", required=True)
    p.add_argument("--out-fasta",   required=True)
    p.add_argument("--mismatches",  type=int, default=1)
    args = p.parse_args()
    build_fasta(args.primers_csv, args.out_fasta, args.mismatches)
