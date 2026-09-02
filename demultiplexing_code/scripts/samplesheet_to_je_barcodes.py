#!/usr/bin/env python3
"""
samplesheet_to_je_barcodes.py — Convert Plate1_SampleSheet.txt to
Je-demultiplex barcode file format.

Je barcode format (one sample per line):
    Sample_ID:I7_seq:I5_seq

Galaxy equivalent:
    Je Demultiplex → "Barcode file" input (tab or colon-delimited)

Usage:
    python3 samplesheet_to_je_barcodes.py \
        --sample-sheet Plate1_SampleSheet.txt \
        --out-file results/03_demux_index/je_barcodes.txt
"""

import argparse
import pandas as pd


def convert(sample_sheet: str, out_file: str):
    df = pd.read_csv(sample_sheet, sep="\t")
    with open(out_file, "w") as f:
        f.write("# Sample_ID:I7_index:I5_index\n")
        for _, row in df.iterrows():
            f.write(f"{row['Sample_ID']}:{row['index']}:{row['index2']}\n")
    print(f"Wrote {len(df)} barcodes → {out_file}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--out-file",     required=True)
    args = p.parse_args()
    convert(args.sample_sheet, args.out_file)
