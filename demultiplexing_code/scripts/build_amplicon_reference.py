#!/usr/bin/env python3
"""
build_amplicon_reference.py — Create a BWA-indexable FASTA from amplicon sequences.

If you do not have a whole-genome reference, use this script to build a
"mini-reference" consisting of the amplicon sequences. Alignment is then
done against this custom reference instead.

Usage:
    python3 build_amplicon_reference.py \
        --primers-csv validated_primers_k20-seed18_with_adaptors.csv \
        --out-fasta   reference/amplicons.fa

    # Then index:
    bwa index reference/amplicons.fa
    samtools faidx reference/amplicons.fa
"""

import argparse
import os
import subprocess
import logging

import pandas as pd

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s")
log = logging.getLogger(__name__)


def build_reference(primers_csv: str, out_fasta: str):
    df = pd.read_csv(primers_csv)
    os.makedirs(os.path.dirname(out_fasta) or ".", exist_ok=True)

    written = 0
    with open(out_fasta, "w") as f:
        for _, row in df.iterrows():
            seq_id = row["Seq_ID"].replace(" ", "_").replace("/", "-")
            seq    = row["Amplicon_Sequence"].strip().upper()
            if not seq:
                log.warning("Empty sequence for %s — skipping", seq_id)
                continue
            # Write 60-char wrapped FASTA
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
            written += 1

    log.info("Wrote %d amplicon sequences → %s", written, out_fasta)

    # Index
    log.info("Indexing with BWA...")
    subprocess.run(["bwa", "index", out_fasta], check=True)
    log.info("Indexing with samtools faidx...")
    subprocess.run(["samtools", "faidx", out_fasta], check=True)

    # Create sequence dictionary for GATK/bcftools
    dict_path = out_fasta.replace(".fa", ".dict").replace(".fasta", ".dict")
    subprocess.run(
        ["samtools", "dict", out_fasta, "-o", dict_path], check=True
    )
    log.info("Reference ready: %s", out_fasta)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--primers-csv", required=True)
    p.add_argument("--out-fasta",   default="reference/amplicons.fa")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_reference(args.primers_csv, args.out_fasta)
