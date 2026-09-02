#!/usr/bin/env python3
"""
demux_by_index.py — Dual-index demultiplexing using I1/I2 index FASTQ files.

Matching rules (confirmed from data):
  - I1 (I7): full 8bp exact match (no mismatch)
  - I2 (I5): strip first base of I2 read → remaining 7bp matches SampleSheet index2 (7bp)

The I2 read is 8bp but has a leading extra base (T), so we use I2[1:] = 7bp.
SampleSheet uses the correct 7bp A5xx index sequences.

Reads R1/R2/I1/I2 in lockstep (all have identical read counts).
Outputs per-sample R1/R2 FASTQ pairs for downstream fastp merging.

Usage:
    python3 demux_by_index.py \
        --r1           data/raw/Undetermined_S0_L001_R1_001.fastq.gz \
        --r2           data/raw/Undetermined_S0_L001_R2_001.fastq.gz \
        --i1           data/raw/Undetermined_S0_L001_I1_001.fastq.gz \
        --i2           data/raw/Undetermined_S0_L001_I2_001.fastq.gz \
        --sample-sheet Plate1_SampleSheet2.txt \
        --out-dir      results/02_demux_index \
        --mismatches   0
"""

import argparse
import gzip
import logging
import os
from collections import defaultdict

import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(levelname)s %(message)s")
log = logging.getLogger(__name__)

REPORT_INTERVAL = 500_000


def load_sample_sheet(path):
    df = pd.read_csv(path, sep="\t")
    df["index"]  = df["index"].str.upper()
    #df["index2"] = df["index2"].str.upper()  
    df["index2"] = df["index2"].str.upper().str[:7] # use positions 1-7 from Plate1_SampleSheet index2 

    i7_map = defaultdict(dict)
    for _, row in df.iterrows():
        i7_map[row["index"]][row["index2"]] = row["Sample_ID"]

    i7_list = df["index"].unique().tolist()
    log.info("Loaded %d samples | %d unique I7 | %d unique I5",
             len(df), df["index"].nunique(), df["index2"].nunique())
    log.info("I7 sequences : %s", ", ".join(sorted(i7_list)))
    log.info("I5 sequences : %s", ", ".join(sorted(df["index2"].unique())))
    return dict(i7_map), i7_list


def open_fq(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_records(fh):
    while True:
        h = fh.readline().rstrip()
        if not h:
            break
        yield h, fh.readline().rstrip(), fh.readline().rstrip(), fh.readline().rstrip()


def demux(r1_path, r2_path, i1_path, i2_path, sample_sheet, out_dir, mismatches):
    os.makedirs(out_dir, exist_ok=True)
    i7_map, i7_list = load_sample_sheet(sample_sheet)

    stats        = defaultdict(int)
    fail_reasons = defaultdict(int)
    i7_counts    = defaultdict(int)
    i5_counts    = defaultdict(int)
    handles      = {}

    def get_handles(sid):
        if sid not in handles:
            handles[sid] = (
                gzip.open(os.path.join(out_dir, f"{sid}_R1.fastq.gz"), "wt"),
                gzip.open(os.path.join(out_dir, f"{sid}_R2.fastq.gz"), "wt"),
            )
        return handles[sid]

    undet_r1 = gzip.open(os.path.join(out_dir, "undetermined_R1.fastq.gz"), "wt")
    undet_r2 = gzip.open(os.path.join(out_dir, "undetermined_R2.fastq.gz"), "wt")

    log.info("Starting demux (lockstep R1/R2/I1/I2)...")
    log.info("  R1: %s", r1_path)
    log.info("  R2: %s", r2_path)
    log.info("  I1 (I7): %s", i1_path)
    log.info("  I2 (I5): %s", i2_path)
    log.info("  Mismatches: %d (I7 exact match only)", mismatches)

    total = 0
    with open_fq(r1_path) as f1, open_fq(r2_path) as f2, \
         open_fq(i1_path) as fi1, open_fq(i2_path) as fi2:

        for (h1, s1, p1, q1), (h2, s2, p2, q2), \
            (_, si1, _, _), (_, si2, _, _) in zip(
                read_records(f1), read_records(f2),
                read_records(fi1), read_records(fi2)):

            total += 1
            obs_i7 = si1.strip().upper()        # full 8bp from I1
            obs_i5 = si2.strip().upper()[1:]    # strip leading base, use 7bp in Undetermined_S0_L001_I2_001.fastq.gz

            i7_counts[obs_i7] += 1
            i5_counts[obs_i5] += 1

            sid = None

            if obs_i7 in i7_map:
                if obs_i5 in i7_map[obs_i7]:
                    sid = i7_map[obs_i7][obs_i5]
                else:
                    fail_reasons["I7_ok_I5_fail"] += 1
            else:
                fail_reasons["I7_fail"] += 1

            if sid:
                fh1, fh2 = get_handles(sid)
                fh1.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                fh2.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                stats[sid] += 1
            else:
                undet_r1.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                undet_r2.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
                stats["undetermined"] += 1

            if total % REPORT_INTERVAL == 0:
                assigned = total - stats["undetermined"]
                log.info("  Progress: %d reads | assigned %d (%.1f%%)",
                         total, assigned, 100 * assigned / total)

    for fh1, fh2 in handles.values():
        fh1.close(); fh2.close()
    undet_r1.close(); undet_r2.close()

    # ── Final report ──────────────────────────────────────────────────────────
    assigned = total - stats["undetermined"]
    log.info("=" * 60)
    log.info("DEMUX COMPLETE")
    log.info("  Total read pairs : %d", total)
    log.info("  Assigned         : %d (%.1f%%)", assigned, 100 * assigned / max(total, 1))
    log.info("  Undetermined     : %d (%.1f%%)",
             stats["undetermined"], 100 * stats["undetermined"] / max(total, 1))
    log.info("")
    log.info("Fail reasons:")
    for reason, n in sorted(fail_reasons.items(), key=lambda x: -x[1]):
        log.info("  %-30s: %d (%.1f%%)", reason, n, 100 * n / max(total, 1))
    log.info("")
    log.info("Per-sample read counts:")
    for sid, n in sorted(stats.items(), key=lambda x: -x[1]):
        if sid == "undetermined":
            continue
        log.info("  %-20s: %d (%.1f%%)", sid, n, 100 * n / max(total, 1))

    stats_path = os.path.join(out_dir, "demux_index_stats.tsv")
    with open(stats_path, "w") as f:
        f.write("Sample_ID\tRead_Pairs\tPct_Total\n")
        for sid, n in sorted(stats.items()):
            pct = round(100 * n / max(total, 1), 2)
            f.write(f"{sid}\t{n}\t{pct}\n")
    log.info("")
    log.info("Stats → %s", stats_path)
    log.info("=" * 60)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--r1",           required=True)
    p.add_argument("--r2",           required=True)
    p.add_argument("--i1",           required=True)
    p.add_argument("--i2",           required=True)
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--out-dir",      required=True)
    p.add_argument("--mismatches",   type=int, default=0)
    args = p.parse_args()
    demux(args.r1, args.r2, args.i1, args.i2,
          args.sample_sheet, args.out_dir, args.mismatches)