# GrowthDependentLysis

**Code and processed data associated with the study:**

*Growth-dependent killing can preserve microbial diversity during antibiotic treatment*

## Overview

This repository contains analysis code and processed data associated with a study investigating how growth-dependent β-lactam killing influences microbial population dynamics and community diversity.

The study integrates clinical metagenomic analysis, strain-resolved dual-barcoding measurements, synthetic microbial community experiments, and ecological modeling.

Two defined microbial communities were examined:

- A 57-member *Escherichia coli* Keio community consisting of genetically barcoded strains.
- An 11-member community assembled from bacterial isolates collected from hospital sink environments.

The computational workflows include sequencing demultiplexing, barcode quantification, strain-resolved abundance estimation, peak-to-trough ratio (PTR) calibration, effective growth- and lysis-rate inference, Bayesian regression, statistical analyses, and ecological simulations.

## 1. Repository structure

| Directory | Description |
|---|---|
| `comm_Keio/` | Processed data and analysis code for the Keio synthetic community experiments |
| `comm_sink/` | Processed data and analysis code for the sink-isolate synthetic community experiments |
| `demultiplexing_code/` | Amplicon sequencing demultiplexing and barcode-processing workflows |
| `growth_lysis_analysis/` | Effective growth-rate estimation, lysis-rate inference, community dynamics, and Bayesian regression |
| `simulation_code/` | Generalized Lotka–Volterra modeling and ecological simulations |

The repository contains independent analysis workflows. Additional scripts, notebooks, processed datasets, and documentation are organized within the corresponding directories.

## 2. System requirements

### Software

The analyses use R, Python, Galaxy, and established bioinformatics tools.

| Software | Version | Purpose |
|---|---|---|
| R | 4.2.2 | Statistical analysis, data processing, and visualization |
| Python | 3.9 | Data processing, Bayesian regression, and simulations |
| Galaxy Main | 26.1.2.dev0 (server version on October 12, 2026) | Barcode sequencing and demultiplexing |
| BBMap | 39.52 | Metagenomic read processing |
| Kraken2 | 2.17.1 | Metagenomic taxonomic classification |
| Bracken | 3.1 | Species-level abundance estimation |
| FastQC | See pipeline environment | Sequencing quality assessment |
| fastp | See pipeline environment | Adapter trimming and paired-end read merging |
| Cutadapt | See pipeline environment | Barcode demultiplexing |
| BWA-MEM | See pipeline environment | Amplicon read alignment |
| SAMtools | See pipeline environment | Alignment statistics and read counting |
| MultiQC | See pipeline environment | Quality-control reporting |
| NumPyro / JAX | See Bayesian pipeline dependencies | Bayesian regression and MCMC sampling |
| Jupyter | See simulation environment | Simulation notebook execution |

The Galaxy Main version listed above refers to the server version observed on October 12, 2026, and may differ from the version used during the original analyses. The exported Galaxy workflow files contain tool-specific version information.

### Operating system and hardware

**Operating system:** macOS for personal-computer analyses; Linux-based environments for applicable sequencing workflows.

**Specialized hardware:** No specialized hardware is required for standard downstream analyses. A personal computer with sufficient memory and storage can be used.

**Typical installation time:** A few minutes for most analysis environments, depending on internet speed and package availability.

### Expected runtime

| Analysis | Approximate runtime |
|---|---|
| Standard R/Python data processing and statistical analysis | A few minutes |
| Ecological simulations | A few minutes under typical configurations |
| Amplicon sequencing processing and demultiplexing | Approximately 2 hours |
| Full Bayesian regression pipeline | Approximately 4–5 hours |

Runtime depends on dataset size, sequencing depth, computational resources, and model settings.

The Bayesian runtime refers to the complete growth–lysis regression pipeline using the default MCMC settings.

## 3. Installation

### Download the repository

Clone the repository:

```bash
git clone https://github.com/youlab/GrowthDependentLysis.git
cd GrowthDependentLysis
```

Alternatively, download the repository as a ZIP archive from GitHub.

### Install analysis-specific dependencies

Different analysis components have different software requirements.

**R analyses**

Install the R packages required by the corresponding scripts. R 4.2.2 was used for the study analyses.

**Python and Bayesian regression**

The Bayesian pipeline can be installed using a Python virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Run these commands from the Bayesian pipeline directory containing its Python package configuration file.

**Amplicon sequencing analysis**

Create a Conda environment containing the required tools:

```bash
conda create -n amplicon_pipeline \
  -c conda-forge -c bioconda \
  python=3 fastqc fastp cutadapt bwa samtools multiqc \
  pandas pysam -y

conda activate amplicon_pipeline
```

The Galaxy workflows can be imported into a compatible Galaxy instance, such as [Galaxy Main](https://usegalaxy.org/).

**Ecological simulations**

The simulation notebooks require Python, Jupyter, and the dependencies imported by the individual notebooks.

Consult the corresponding subdirectory documentation for additional requirements.

## 4. Running the analyses

The repository is organized into independent analysis components. Users can run individual workflows without executing the entire repository.

### 4.1 Amplicon sequencing and barcode processing

**Directory:** `demultiplexing_code/`

The pipeline processes paired-end Illumina sequencing data through the following stages:

1. Raw sequencing quality assessment using FastQC.
2. Adapter trimming and paired-end merging using fastp.
3. Sample demultiplexing using i5/i7 index sequences.
4. Amplicon demultiplexing using target-specific barcodes and Cutadapt.
5. Alignment to amplicon references using BWA-MEM.
6. Read counting using SAMtools idxstats.
7. Quality-control aggregation using MultiQC.

The primary output is a sample-by-amplicon read-count matrix:

`results/06_counts/count_matrix.tsv`

Additional outputs include long-format count tables, amplicon summaries, sample-level read totals, and quality-control reports.

**Expected runtime:** Approximately 2 hours for a typical sequencing dataset.

Detailed configuration and execution instructions are provided in `demultiplexing_code/`.

### 4.2 Keio synthetic community analysis

**Directory:** `comm_Keio/`

This directory contains code and processed data for the 57-member Keio synthetic community.

The analyses include barcode-based community composition estimation, strain-resolved abundance measurements, and microbial community responses to β-lactam treatment.

The Galaxy workflow used for Keio barcode processing is provided at:

`comm_Keio/Galaxy-Workflow-Chromosome_BarcodeProcessing_StrainBarcodes_Paired_End_w__counting_(counting_R1).ga`

### 4.3 Sink-isolate community analysis

**Directory:** `comm_sink/`

This directory contains code and processed data for the 11-member sink-isolate community.

The analyses include strain-specific barcode quantification, community composition estimation, and antibiotic-response measurements.

### 4.4 Growth- and lysis-rate inference

**Directory:** `growth_lysis_analysis/`

Strain-specific effective growth rates are estimated using calibrated PTR measurements. Net population rates are obtained from fitted time-resolved biomass trajectories, and lysis rates are inferred from the difference between effective growth and net population rates.

The directory also contains a Bayesian regression pipeline that propagates uncertainty through the growth- and lysis-rate inference procedure.

The Bayesian pipeline consists of three stages:

1. **Model 1:** Bayesian linear regression of effective growth rate against ln(PTR).
2. **Model 2:** Bayesian cubic polynomial regression of log-transformed biomass against time.
3. **Model 4:** Stacked Bayesian regression combining posterior information from Models 1 and 2 to infer lysis rates and estimate the growth–lysis relationship.

Model 4 in the code corresponds to Model 3 described in Supplementary Note 4.

The Bayesian pipeline includes convergence diagnostics based on effective sample size (ESS) and split R-hat.

**Expected runtime:** Approximately 4–5 hours for the complete analysis of 57 strains using the default MCMC settings.

Detailed installation, execution, model parameters, output structure, and convergence-diagnostic instructions are provided in the Bayesian pipeline README.

### 4.5 Ecological simulations

**Directory:** `simulation_code/`

The simulations use a generalized Lotka–Volterra framework to examine how growth-dependent lysis affects microbial community dynamics and diversity.

The directory contains three Jupyter notebooks:

| Notebook | Description | Related results |
|---|---|---|
| `gLV_simulation_Shannon.ipynb` | Effects of growth–lysis coupling on Shannon diversity | Supplementary Fig. S15 |
| `gLV_simulation_debris.ipynb` | Effects of biomass residue on lysis-rate inference | Supplementary Fig. S9 |
| `gLV_simulation_net_rate_fitting.ipynb` | Accuracy of net-rate estimation using polynomial fitting | Supplementary Fig. S7 |

Model equations, parameter definitions, and simulation assumptions are described in Supplementary Notes 1–3.

**Expected runtime:** Typically a few minutes under standard simulation settings.

## 5. Input data and expected outputs

Different workflows require different input data formats.

| Analysis | Input | Expected output |
|---|---|---|
| Amplicon sequencing | Paired-end FASTQ files, sample sheet, primer/amplicon references | Barcode read-count matrices and QC reports |
| Keio community | Barcode counts and experimental measurements | Strain abundance and community-response analyses |
| Sink community | Barcode counts and experimental measurements | Species abundance and community-response analyses |
| Growth–lysis analysis | PTR calibration data and time-resolved biomass measurements | Growth-rate, net-rate, and lysis-rate estimates |
| Bayesian regression | PTR–growth calibration data and strain-resolved biomass trajectories | Posterior estimates, diagnostic plots, and summary tables |
| Ecological simulations | Model parameters and generated initial conditions | Simulated population trajectories, regression results, and diversity measurements |

Processed data and analysis scripts are organized within the corresponding directories.

## 6. Reproducing manuscript analyses

The repository contains computational workflows associated with the main and supplementary analyses.

| Analysis | Code location | Manuscript reference |
|---|---|---|
| Barcode sequencing and quantification | `comm_Keio/`, `demultiplexing_code/` | Methods: Amplicon sequencing and data processing |
| PTR calibration and abundance estimation | `comm_Keio/`, `comm_sink/` | Fig. 3 |
| Growth- and lysis-rate inference | `growth_lysis_analysis/` | Figs. 4–6 |
| Bayesian growth–lysis regression | `growth_lysis_analysis/` | Supplementary Fig. S10; Note 4 |
| Net-rate fitting validation | `simulation_code/gLV_simulation_net_rate_fitting.ipynb` | Supplementary Fig. S7; Note 2 |
| Biomass residue simulations | `simulation_code/gLV_simulation_debris.ipynb` | Supplementary Fig. S9; Note 3 |
| Shannon diversity simulations | `simulation_code/gLV_simulation_Shannon.ipynb` | Supplementary Fig. S15; Note 1 |

The full Bayesian regression pipeline can be executed using the input datasets and commands provided in its README.

For computational simulations, the model equations and parameter settings are described in Supplementary Notes 1–3. Stochastic simulations may exhibit numerical variation across executions depending on random seeds.

## 7. Demonstration and reproducibility

The repository contains analysis-specific input files and workflows that can be used to demonstrate the computational methods.

For the Bayesian regression pipeline, the default input datasets and execution commands are documented in the corresponding README. The full analysis requires approximately 4–5 hours.

For the amplicon sequencing pipeline, example input formats, processing steps, and expected count-matrix outputs are documented in `demultiplexing_code/`.

For ecological simulations, the Jupyter notebooks provide executable implementations of the mathematical models.

Users should consult the corresponding subdirectory documentation for input requirements, software dependencies, execution instructions, and expected outputs.

For reproducibility, the original software versions, model parameters, MCMC settings, and random seeds should be retained whenever available.

## 8. Code contributions and acknowledgments

The sequencing demultiplexing and barcode-processing workflow for the Keio community was implemented in Galaxy and is available at:

`comm_Keio/Galaxy-Workflow-Chromosome_BarcodeProcessing_StrainBarcodes_Paired_End_w__counting_(counting_R1).ga`

The dual-barcoding design workflow for the sink isolates was adapted from the [xiaoli_PCR_design](https://github.com/jaeminkim518/xiaoli_PCR_design) workflow developed by Jaemin Kim.

The version used for the final sink-isolate primer design includes modifications to the ori-ter amplicon design workflow:

- Fixed amplicon length: 105 bp
- `seed_len = 4`
- `max_mismatches = 3`

The sequencing demultiplexing and barcode-processing workflow for the sink community was developed in collaboration with **Jiesi Lei** and is available in `demultiplexing_code/`.

## 9. Data and code availability

Custom analysis code and processed data associated with the experimental and computational analyses are available in this repository:

https://github.com/youlab/GrowthDependentLysis

The clinical metagenomic dataset was provided by Nicholas A. Turner. Availability of clinical sequencing data and associated patient-level metadata is subject to confirmation of the applicable institutional and ethical data-sharing requirements.

Information regarding clinical data access will be provided in the manuscript's Data Availability Statement.

## 10. License

No specific software license has been assigned to this repository. The source code is publicly accessible for inspection and reproducibility purposes.

## 11. Contact

For questions regarding the study, code, or data, please contact:

**Lingchong You**  
Department of Biomedical Engineering  
Duke University  
lingchong.you@duke.edu
