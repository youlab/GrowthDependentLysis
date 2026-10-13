# Amplicon Sequencing Demo Dataset

This directory contains a small demonstration dataset for testing the amplicon sequencing and barcode-processing pipeline used in the GrowthDependentLysis study.

## Dataset description

The demo dataset contains 1,000 synchronized Illumina read groups subsampled from the original sequencing dataset, including paired-end reads (R1/R2) and corresponding index reads (I1/I2).

The dataset is intended to demonstrate pipeline functionality rather than reproduce the quantitative results of the full sequencing analysis.

## Input files

| File | Description |
|---|---|
| `data/raw/*R1*.fastq.gz` | Forward sequencing reads |
| `data/raw/*R2*.fastq.gz` | Reverse sequencing reads |
| `data/raw/*I1*.fastq.gz` | Index 1 reads |
| `data/raw/*I2*.fastq.gz` | Index 2 reads |
| `Plate1_SampleSheet.txt` | Sample-specific i5/i7 indices |
| `sink_dual_amplicons.csv` | Target barcodes, primer sequences, and amplicon references |
| `reference/amplicons.fa` | Reference sequences for alignment |

## Running the demo

The demonstration uses the same processing workflow as the full amplicon sequencing analysis.

1. Install the dependencies described in the parent directory's README.
2. Configure the pipeline to use the files in this demo directory.
3. Run the pipeline using the instructions provided in the parent directory.
4. Examine the generated count matrices and quality-control reports.

## Expected output

The pipeline generates a sample-by-amplicon read-count matrix and associated quality-control summaries.

Because the dataset contains only 1,000 read groups, some samples and amplicons may have zero counts.

## Runtime

The demonstration is expected to run substantially faster than the full sequencing dataset. The actual runtime should be measured after testing.

## Reproducibility

The original sequencing analysis requires the complete dataset and corresponding processing parameters. This reduced dataset is intended for testing the computational workflow only.
