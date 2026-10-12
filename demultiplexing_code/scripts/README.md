# Amplicon Sequencing Analysis Pipeline

A reproducible, Galaxy-compatible pipeline for processing paired-end Illumina amplicon sequencing data, demultiplexing samples and strain-specific barcode targets, aligning reads to amplicon references, and generating read-count matrices for downstream analysis.

## Overview

This pipeline processes paired-end Illumina sequencing data through seven main steps:

1. Quality assessment of raw sequencing reads.
2. Adapter trimming and paired-end read merging.
3. Sample demultiplexing using Illumina i5/i7 indices.
4. Amplicon demultiplexing using target-specific barcode sequences.
5. Alignment to amplicon reference sequences.
6. Read counting and count-matrix generation.
7. Quality-control report aggregation.

The pipeline supports both command-line execution and equivalent analysis steps using Galaxy tools. Custom Python scripts are provided for input preparation, demultiplexing, reference construction, and count-matrix generation.

## 1. Software requirements

### Dependencies

| Step | Command-line tool | Galaxy equivalent | Purpose |
|---|---|---|---|
| 1 | FastQC | FastQC | Raw read quality assessment |
| 2 | fastp | fastp | Adapter trimming and paired-end merging |
| 3 | `demux_by_index.py` | Je Demultiplex | Sample demultiplexing using i5/i7 indices |
| 4 | Cutadapt | Cutadapt | Target-specific barcode demultiplexing |
| 5 | BWA-MEM | Map with BWA-MEM | Alignment to amplicon references |
| 6 | SAMtools | SAMtools idxstats | Read counting |
| 7 | MultiQC | MultiQC | Quality-control report aggregation |

Additional Python dependencies:

- Python 3
- pandas
- pysam

### Operating system and hardware

The pipeline is designed for execution in a Unix-like environment using Conda-compatible software packages.

A standard multicore CPU system is sufficient for the analysis. Memory and storage requirements depend on the sequencing depth, number of samples, and number of target amplicons.

No GPU is required.

**Tested environment:** [Specify the operating system, Python version, and software versions used for the final analysis.]

**Installation time:** [Specify the approximate time required to create the Conda environment.]

## 2. Installation

### 2.1 Clone the repository

```bash
git clone https://github.com/youlab/GrowthDependentLysis.git
cd GrowthDependentLysis
```

Navigate to the directory containing the amplicon sequencing pipeline.

### 2.2 Create the Conda environment

```bash
conda create -n amplicon_pipeline \
  -c conda-forge -c bioconda \
  python=3 \
  fastqc fastp cutadapt bwa samtools multiqc \
  pandas pysam -y
```

Activate the environment:

```bash
conda activate amplicon_pipeline
```

The versions used in the final analysis should be recorded for reproducibility. Where possible, use an exported Conda environment file to recreate the original analysis environment.

### 2.3 Verify installation

```bash
fastqc --version
fastp --version
cutadapt --version
bwa
samtools --version
multiqc --version
python --version
```

## 3. Input data

The pipeline requires paired-end FASTQ files, an Illumina sample sheet, and a primer/amplicon reference table.

### 3.1 Raw sequencing reads

Input files should contain paired-end Illumina reads in FASTQ or compressed FASTQ format.

Example:

```text
reads_R1.fastq.gz
reads_R2.fastq.gz
```

The paired-end reads must originate from the same sequencing library and maintain consistent read identifiers.

### 3.2 Illumina sample sheet

Example: `Plate1_SampleSheet.txt`

```text
Sample_ID   I7_Index_ID   index      I5_Index_ID   index2
1-1         D701          CGAGTAAT   D501          TATAGCCT
```

The sample sheet specifies the i7 and i5 index sequences used to assign reads to individual samples.

The custom script `samplesheet_to_je_barcodes.py` converts the sample sheet into a barcode file compatible with Je Demultiplex.

Example converted entry:

```text
1-1:CGAGTAAT:TATAGCCT
```

### 3.3 Primer and amplicon reference table

The primer CSV file contains the sequences required for target-specific barcode identification and reference construction.

Required fields include:

```text
Seq_ID,Target_Barcode,Forward_Primer,Reverse_Primer,Amplicon_Sequence
```

The following scripts process this information:

- `primers_to_cutadapt_barcodes.py`: generates barcode definitions for Cutadapt.
- `build_amplicon_reference.py`: generates the reference FASTA file for BWA-MEM alignment.

The `Target_Barcode` and `Forward_Primer` fields are used to construct linked-adapter definitions for barcode demultiplexing.

The `Amplicon_Sequence` field is used to construct the alignment reference.

## 4. Pipeline configuration

The pipeline uses a shell-based configuration and execution workflow.

Before running the analysis, verify that the configuration contains the correct paths to:

- Raw paired-end FASTQ files.
- Illumina sample sheet.
- Primer and amplicon CSV file.
- Reference FASTA file.
- Output directories.
- Custom Python scripts.

Review the configuration in `amplicon_pipeline.sh` and any associated configuration files before execution.

All paths should be adjusted to the local computational environment.

## 5. Running the pipeline

### Step 1. Activate the environment

```bash
conda activate amplicon_pipeline
```

### Step 2. Generate the amplicon reference

From the analysis directory containing the primer CSV file:

```bash
python3 scripts/build_amplicon_reference.py \
  --primers-csv sink_dual_amplicons.csv \
  --out-fasta reference/amplicons.fa
```

This step generates the reference sequences used for downstream read alignment.

### Step 3. Execute the pipeline

```bash
chmod +x scripts/amplicon_pipeline.sh
bash scripts/amplicon_pipeline.sh all
```

The `all` argument runs the configured workflow, including quality control, read processing, demultiplexing, alignment, counting, and report generation.

Ensure that the script paths and input directories are configured correctly before execution.

## 6. Pipeline workflow

```text
Paired-end FASTQ files (R1 + R2)
              |
              v
Step 1. FastQC
        Raw read quality assessment
              |
              v
Step 2. fastp
        Adapter trimming and paired-end merging
              |
              v
Step 3. Sample demultiplexing
        i5/i7 index-based sample assignment
        Je Demultiplex / demux_by_index.py
              |
              v
Step 4. Amplicon demultiplexing
        Target-specific barcode identification
        Cutadapt linked adapters
              |
              v
Step 5. BWA-MEM
        Alignment to amplicon reference sequences
              |
              v
Step 6. SAMtools idxstats
        Read counting and matrix generation
              |
              v
Step 7. MultiQC
        Aggregated quality-control reporting
```

The number of samples and amplicon targets depends on the sample sheet and primer reference table supplied for each sequencing run.

## 7. Galaxy-compatible workflow

The pipeline can also be implemented using Galaxy tools.

### Sample demultiplexing

1. Upload `Plate1_SampleSheet.txt`.
2. Convert the sample sheet using `samplesheet_to_je_barcodes.py`.
3. Provide the generated barcode file to Je Demultiplex.
4. Demultiplex reads according to sample-specific i5/i7 indices.

### Amplicon demultiplexing

1. Upload the primer CSV file.
2. Generate barcode definitions using `primers_to_cutadapt_barcodes.py`.
3. Provide the resulting `barcodes.fasta` to Cutadapt.
4. Demultiplex reads according to target-specific barcode sequences.

### Reference-based alignment

1. Generate `amplicons.fa` using `build_amplicon_reference.py`.
2. Upload the reference FASTA file to Galaxy.
3. Align amplicon reads using BWA-MEM.
4. Obtain alignment statistics using SAMtools idxstats.
5. Aggregate counts across samples using the provided count-merging workflow.

The Galaxy workflow and command-line pipeline use corresponding analysis stages, although implementation-specific parameters should be checked for consistency.

## 8. Output files

The pipeline generates the following primary outputs:

| Output file | Description |
|---|---|
| `results/06_counts/count_matrix.tsv` | Primary amplicon read-count matrix |
| `results/06_counts/count_matrix_long.tsv` | Long-format read-count table |
| `results/06_counts/amplicon_summary.tsv` | Summary statistics for each amplicon |
| `results/06_counts/sample_totals.tsv` | Total mapped reads per sample |
| `results/07_report/multiqc_report.html` | Aggregated quality-control report |

### Primary count matrix

The primary output is a tab-separated matrix containing read counts for each amplicon across samples.

Example:

```text
Seq_ID          1-1    1-2    1-3
A. baylyi ADP1  1523   1847   902
DA28102_1       3201   2984   3456
```

Rows represent target amplicons, and columns represent sequencing samples.

The resulting matrix can be imported into R or Python for downstream analyses, including barcode quantification, relative-abundance estimation, and growth-rate inference.

## 9. Demonstration and expected results

A small example dataset is recommended for validating installation and demonstrating pipeline functionality.

**Demo dataset:** [Provide the location of a small example FASTQ dataset, sample sheet, and primer CSV.]

**Demo command:** [Provide the command used to execute the demonstration.]

**Expected output:** The demonstration should generate a read-count matrix, associated summary tables, and a MultiQC report.

**Expected runtime:** [Record the runtime measured on a standard desktop or workstation.]

The demo should allow users to verify that the pipeline executes successfully without requiring the complete experimental sequencing dataset.

## 10. Reproducibility

To support reproducibility, users should record:

- Software and package versions.
- Input sequencing files and associated sample sheets.
- Primer and reference sequences.
- Demultiplexing parameters.
- Alignment and read-counting parameters.
- Pipeline configuration files.
- Analysis outputs and quality-control reports.

For analyses associated with the manuscript, the resulting read-count matrices are used in downstream strain-resolved abundance and growth–lysis analyses.

The mathematical framework and downstream analysis methods are described in the manuscript and Supplementary Information.

## 11. Troubleshooting

| Problem | Possible cause | Suggested action |
|---|---|---|
| Conda installation fails | Package conflicts or channel configuration | Check channel priorities and dependency compatibility |
| No reads assigned to samples | Incorrect index sequences or index orientation | Verify sample-sheet indices and demultiplexing settings |
| Low amplicon assignment | Barcode mismatch or adapter configuration | Check barcode sequences and Cutadapt parameters |
| Low alignment rate | Incorrect reference sequences or incomplete trimming | Verify amplicon references and inspect read quality |
| Missing count-matrix entries | Missing or incomplete alignment files | Check BAM files and SAMtools idxstats outputs |
| MultiQC report is incomplete | Missing QC output files | Confirm that upstream QC steps completed successfully |

## 12. Code availability

The pipeline is part of the code repository associated with the study:

**GrowthDependentLysis**

https://github.com/youlab/GrowthDependentLysis

The repository contains code and processed data associated with microbial community experiments, barcode sequencing, growth–lysis analysis, and ecological simulations.

## 13. Contact

For questions regarding the pipeline, code, or associated analyses, please contact:

Lingchong You  
Department of Biomedical Engineering  
Duke University  
lingchong.you@duke.edu
