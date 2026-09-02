#!/usr/bin/env bash
# =============================================================================
# Amplicon Sequencing Pipeline — Galaxy-compatible tool stack
# =============================================================================
#
# Each step uses the same tool Galaxy uses, so this script and a Galaxy
# workflow are directly interchangeable.
#
# Galaxy tool mapping:
#   Step 1  FastQC              → toolshed: fastqc
#   Step 2  fastp               → toolshed: fastp
#   Step 3  Index demux (I7+I5) → toolshed: je_demultiplex (custom fallback)
#   Step 4  BWA-MEM             → toolshed: bwa_mem
#   Step 5  SAMtools flagstat   → toolshed: samtools_flagstat
#           SAMtools idxstats   → toolshed: samtools_idxstats   ← read count matrix
#           SAMtools sort/index → toolshed: samtools_sort / samtools_index
#   Step 6  MultiQC             → toolshed: multiqc
#
# NOTE: No barcode demux step — each amplicon has a unique sequence (barcode),
# BWA-MEM alignment to the amplicon reference handles species identification.
# samtools idxstats counts reads per contig = reads per species per sample.
#
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.sh"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE"

# ── Defaults ──────────────────────────────────────────────────────────────────
: "${RAW_DIR:=data/raw}"
: "${OUT_DIR:=results}"
: "${SAMPLE_SHEET:=Plate1_SampleSheet.txt}"
: "${PRIMERS_CSV:=sink_dual_amplicons.csv}"
: "${REFERENCE:=reference/amplicons.fa}"
: "${THREADS:=8}"
: "${MIN_QUALITY:=20}"
: "${MIN_LENGTH:=50}"
: "${OVERLAP_MIN:=50}"
: "${MAX_MISMATCH_INDEX:=0}"
: "${ADAPTER_FWD:=AACACTCTTTCCCTACACGACGCTCTTCCGATCT}"
: "${ADAPTER_REV:=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT}"

# ── Output dirs ───────────────────────────────────────────────────────────────
QC_RAW_DIR="${OUT_DIR}/01_fastqc_raw"
DEMUX_INDEX_DIR="${OUT_DIR}/02_demux_index"
FASTP_DIR="${OUT_DIR}/03_fastp"
ALIGN_DIR="${OUT_DIR}/04_aligned"
COUNTS_DIR="${OUT_DIR}/05_counts"
REPORT_DIR="${OUT_DIR}/06_report"

# ── Helpers ───────────────────────────────────────────────────────────────────
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
check_tool() { command -v "$1" &>/dev/null || die "Tool not found: $1"; }

check_tools() {
    log "Checking tools..."
    for t in fastqc fastp bwa samtools multiqc python3; do
        check_tool "$t"
    done
    python3 -c "import pandas" 2>/dev/null || die "pip install pandas"
    log "All tools OK."
}

# =============================================================================
# STEP 1 — FastQC (raw)
# Galaxy tool: FastQC (Galaxy Version 0.74+galaxy0)
# =============================================================================
run_fastqc_raw() {
    log "=== STEP 1: FastQC (raw reads) ==="
    mkdir -p "$QC_RAW_DIR"

    fastqc \
        --threads "$THREADS" \
        --outdir  "$QC_RAW_DIR" \
        "${RAW_DIR}"/*.fastq.gz

    log "FastQC raw → $QC_RAW_DIR"
}

# =============================================================================
# STEP 2 — Index demultiplexing using I1/I2 FASTQ files
# Galaxy tool: Je Demultiplex (je_demultiplex)
#
# Matching rules:
#   I1 (I7): full 8bp exact match (no mismatch)
#   I2 (I5): last 7bp of I2 matched against last 7bp of SampleSheet index2
# Input:  raw R1/R2 + I1/I2 (lockstep)
# Output: per-sample R1/R2 FASTQ pairs
# =============================================================================
run_demux_index() {
    log "=== STEP 2: Dual-index demultiplexing (I1/I2 files) ==="
    mkdir -p "$DEMUX_INDEX_DIR"

    python3 "${SCRIPT_DIR}/scripts/demux_by_index.py" \
        --r1           "${RAW_DIR}/Undetermined_S0_L001_R1_001.fastq.gz" \
        --r2           "${RAW_DIR}/Undetermined_S0_L001_R2_001.fastq.gz" \
        --i1           "${RAW_DIR}/Undetermined_S0_L001_I1_001.fastq.gz" \
        --i2           "${RAW_DIR}/Undetermined_S0_L001_I2_001.fastq.gz" \
        --sample-sheet "$SAMPLE_SHEET" \
        --out-dir      "$DEMUX_INDEX_DIR" \
        --mismatches   0

    log "Index demux → $DEMUX_INDEX_DIR"
}

# =============================================================================
# STEP 3 — fastp: trim + merge per-sample paired-end reads
# Galaxy tool: fastp (Galaxy Version 0.23.4+galaxy0)
#
# Input:  per-sample R1/R2 from Step 2
# Output: per-sample merged FASTQ
# =============================================================================
run_fastp() {
    log "=== STEP 3: fastp (trim + PE merge, per sample) ==="
    mkdir -p "$FASTP_DIR"

    for r1 in "${DEMUX_INDEX_DIR}"/*_R1.fastq.gz; do
        [[ -f "$r1" ]] || continue
        r2="${r1/_R1./_R2.}"
        [[ -f "$r2" ]] || { log "WARNING: no R2 for $r1 — skipping"; continue; }

        base=$(basename "$r1" _R1.fastq.gz)
        log "  fastp: $base"

        fastp \
            --in1  "$r1" \
            --in2  "$r2" \
            --adapter_sequence    "$ADAPTER_FWD" \
            --adapter_sequence_r2 "$ADAPTER_REV" \
            --merge \
            --merged_out          "${FASTP_DIR}/${base}_merged.fastq.gz" \
            --out1                "${FASTP_DIR}/${base}_unmerged_R1.fastq.gz" \
            --out2                "${FASTP_DIR}/${base}_unmerged_R2.fastq.gz" \
            --overlap_len_require     "$OVERLAP_MIN" \
            --qualified_quality_phred "$MIN_QUALITY" \
            --length_required     "$MIN_LENGTH" \
            --thread              "$THREADS" \
            --html "${FASTP_DIR}/${base}_fastp.html" \
            --json "${FASTP_DIR}/${base}_fastp.json"
    done

    log "fastp → $FASTP_DIR"
}

# =============================================================================
# STEP 4 — Alignment: BWA-MEM
# Galaxy tool: Map with BWA-MEM (bwa_mem, Galaxy Version 0.7.18)
#
# Input: per-sample merged FASTQs from Step 3 (I7+I5 demux)
# Each read is aligned to the amplicon reference (89 contigs = 89 species).
# BWA automatically assigns reads to the correct contig → species identification.
# =============================================================================
run_alignment() {
    log "=== STEP 4: BWA-MEM alignment ==="
    mkdir -p "$ALIGN_DIR"

    [[ -f "${REFERENCE}.bwt" ]] || {
        log "  Indexing reference..."
        bwa index "$REFERENCE"
        samtools faidx "$REFERENCE"
    }

    for merged_fq in "${FASTP_DIR}"/*_merged.fastq.gz; do
        [[ -f "$merged_fq" ]] || continue
        sample=$(basename "$merged_fq" _merged.fastq.gz)
        out_bam="${ALIGN_DIR}/${sample}.bam"
        log "  Aligning $sample..."

        bwa mem \
            -t "$THREADS" \
            -B 4 \
            -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
            "$REFERENCE" \
            "$merged_fq" | \
        samtools sort -@ "$THREADS" -o "$out_bam"

        samtools index "$out_bam"
    done

    log "Alignment → $ALIGN_DIR"
}

# =============================================================================
# STEP 5 — Read count matrix via SAMtools idxstats
# Galaxy tool: SAMtools idxstats (samtools_idxstats, Galaxy Version 2.0.5)
#
# idxstats columns: contig  length  mapped_reads  unmapped_reads
# One BAM per sample, one row per species contig → count matrix directly.
# =============================================================================
run_counts() {
    log "=== STEP 5: Read count matrix (samtools idxstats) ==="
    mkdir -p "$COUNTS_DIR"

    for bam in "${ALIGN_DIR}"/*.bam; do
        [[ -f "$bam" ]] || continue
        base=$(basename "$bam" .bam)

        samtools idxstats "$bam" > "${COUNTS_DIR}/${base}.idxstats"
        samtools flagstat "$bam" > "${COUNTS_DIR}/${base}.flagstat"
    done

    python3 "${SCRIPT_DIR}/scripts/merge_idxstats.py" \
        --idxstats-dir "$COUNTS_DIR" \
        --primers-csv  "$PRIMERS_CSV" \
        --out-dir      "$COUNTS_DIR"

    log "Count matrix → ${COUNTS_DIR}/count_matrix.tsv"
}

# =============================================================================
# STEP 6 — MultiQC report
# Galaxy tool: MultiQC (multiqc, Galaxy Version 1.11+galaxy1)
# =============================================================================
run_report() {
    log "=== STEP 6: MultiQC report ==="
    mkdir -p "$REPORT_DIR"

    multiqc \
        "$QC_RAW_DIR" \
        "$FASTP_DIR" \
        "$COUNTS_DIR" \
        --outdir   "$REPORT_DIR" \
        --filename "multiqc_report" \
        --title    "Amplicon Sequencing — $(basename "$SAMPLE_SHEET" .txt)"

    log "Report → ${REPORT_DIR}/multiqc_report.html"
}

# =============================================================================
# MAIN
# =============================================================================
usage() {
    cat <<EOF
Usage: $(basename "$0") [STEPS...]

Steps:
  all         Full pipeline (default)
  fastqc      Step 1: FastQC raw
  fastp       Step 2: Trim + PE merge
  demux_idx   Step 3: Index demux (I7+I5 → per-sample FASTQ)
  align       Step 4: BWA-MEM (per-sample → BAM aligned to amplicon ref)
  counts      Step 5: samtools idxstats → count matrix (sample × species)
  report      Step 6: MultiQC

Galaxy tool equivalents are documented inline in each step.
EOF
}

main() {
    local steps=("$@")
    [[ ${#steps[@]} -eq 0 ]] && steps=("all")
    check_tools
    for step in "${steps[@]}"; do
        case "$step" in
            all)       run_fastqc_raw; run_demux_index; run_fastp;
                       run_alignment; run_counts; run_report ;;
            fastqc)    run_fastqc_raw ;;
            fastp)     run_fastp ;;
            demux_idx) run_demux_index ;;
            align)     run_alignment ;;
            counts)    run_counts ;;
            report)    run_report ;;
            help|-h)   usage; exit 0 ;;
            *)         die "Unknown step: $step" ;;
        esac
    done
    log "Done."
}

main "$@"