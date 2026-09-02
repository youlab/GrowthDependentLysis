# =============================================================================
# Pipeline Configuration
# =============================================================================

RAW_DIR="data/raw"
SAMPLE_SHEET="Plate1_SampleSheet.txt"
PRIMERS_CSV="sink_dual_amplicons.csv"

# Build with: python3 scripts/build_amplicon_reference.py
REFERENCE="reference/amplicons.fa"

OUT_DIR="results"
THREADS=8

# Quality
MIN_QUALITY=20
MIN_LENGTH=50
OVERLAP_MIN=50           # fastp min overlap for PE merging

# Demux
MAX_MISMATCH_INDEX=0     # per index (I7 or I5)
MAX_MISMATCH_BARCODE=1   # in Target_Barcode

# Illumina TruSeq adaptors (from your primers CSV)
ADAPTER_FWD="AACACTCTTTCCCTACACGACGCTCTTCCGATCT"
ADAPTER_REV="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
