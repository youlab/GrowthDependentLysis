# GrowthDependentLysis

Code and processed data associated with the study:

**Growth-dependent killing can preserve microbial diversity during antibiotic treatment**

## Overview

This repository contains analysis code and processed data used to investigate how growth-dependent β-lactam killing shapes microbial community dynamics and diversity.

The study combines strain-resolved dual-barcoding measurements and synthetic microbial community experiments to quantify population-specific growth and lysis responses under β-lactam treatment.

## Repository structure

- `comm_Keio/` — Analysis of the Keio synthetic community experiments
- `comm_sink/` — Analysis of the sink-isolate synthetic community experiments
- `demultiplexing_code/` — Sequencing demultiplexing and barcode-processing workflows
- `growth_lysis_analysis/` — Analysis of effective growth rates, lysis rates, and community dynamics
- `simulation_code/` — Mathematical modeling and simulation analyses

## Code contributions and acknowledgments

The sequencing demultiplexing and barcode-processing workflow for the Keio community was implemented in Galaxy and is available at `comm_Keio/Galaxy-Workflow-Chromosome_BarcodeProcessing_StrainBarcodes_Paired_End_w__counting_(counting_R1).ga`.

The dual-barcoding design workflow for the sink isolates was adapted from the
[`xiaoli_PCR_design`](https://github.com/jaeminkim518/xiaoli_PCR_design)
workflow developed by Jaemin Kim. The version used for the final sink-isolate
primer design includes modifications to the ori-ter amplicon design workflow,
including a fixed amplicon length of 105 bp, `seed_len = 4`, and
`max_mismatches = 3`.

The sequencing demultiplexing and barcode-processing workflow for the sink community was developed in collaboration with **Jiesi Lei** and is available in `demultiplexing_code/`.

## Requirements

Software dependencies and package requirements are provided within the corresponding analysis directories where applicable.

## Contact

For questions regarding the study, code, or data, please contact:

**Lingchong You**  
Department of Biomedical Engineering  
Duke University  
lingchong.you@duke.edu
