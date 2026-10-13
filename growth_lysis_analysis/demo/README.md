# Bayesian Growth–Lysis Analysis: Example Dataset

## Overview

This directory documents example input datasets for the Bayesian growth–lysis regression pipeline used in the study *Growth-dependent killing can preserve microbial diversity during antibiotic treatment*.

The datasets are stored in the parent `data/` directory and are used to demonstrate PTR-based growth-rate calibration, biomass trajectory fitting, and growth–lysis regression.

## Example datasets

| File | Description |
|---|---|
| `../data/first_data-lnPTR-growth_rate-plot_57strains.csv` | PTR–growth-rate calibration data for 57 bacterial strains |
| `../data/od_3rep_250511_LB_sensitive_Carbenicillin70_37_1.txt` | Time-resolved biomass and barcode measurements for three biological replicates under one experimental condition |

These files represent example inputs from the study rather than the complete experimental dataset.

## Analysis workflow

The Bayesian pipeline consists of three stages:

1. **Model 1:** Estimate effective growth rates from ln(PTR).
2. **Model 2:** Fit cubic polynomial models to log-transformed biomass trajectories and estimate net population rates.
3. **Model 4:** Combine posterior samples from Models 1 and 2 to estimate lysis rates and their relationship with effective growth rates.

Model 4 in the code corresponds to Model 3 in Supplementary Note 4.

## Running the example

From the `growth_lysis_analysis/` directory, activate the Python environment described in the main README.

### Model 1: PTR–growth-rate calibration

```bash
python strain_results/fit_all_strains.py \
  --data data/first_data-lnPTR-growth_rate-plot_57strains.csv \
  --name example

python strain_results/analyze_strains.py \
  --data data/first_data-lnPTR-growth_rate-plot_57strains.csv \
  --name example
```

### Model 2: Biomass trajectory fitting

```bash
python cubic_growth_rate/fit_all.py \
  --data data/od_3rep_250511_LB_sensitive_Carbenicillin70_37_1.txt \
  --name example

python cubic_growth_rate/analyze_all.py \
  --data data/od_3rep_250511_LB_sensitive_Carbenicillin70_37_1.txt \
  --name example
```

### Model 4: Growth–lysis regression

```bash
python stacked_model/fit_model4.py \
  --data data/od_3rep_250511_LB_sensitive_Carbenicillin70_37_1.txt \
  --name example \
  --model1-run example \
  --model2-run example \
  --time 6.0

python stacked_model/analyze_model4.py \
  --name example
```

### Convergence diagnostics

```bash
python diagnostics/check_convergence.py --name example
```

## Expected outputs

The analysis produces posterior model files, diagnostic plots, and summary tables.

| Directory | Expected output |
|---|---|
| `strain_results/results/example/` | Model 1 posterior results and calibration summaries |
| `cubic_growth_rate/results/example/` | Model 2 fitted trajectories and net-rate estimates |
| `stacked_model/results/example/` | Model 4 growth–lysis regression results |
| `diagnostics/results/example/` | MCMC convergence diagnostics |

## Expected runtime

The complete Bayesian pipeline typically requires approximately 4–5 hours using the default MCMC settings. Runtime depends on hardware, dataset size, and sampling parameters.

These example commands use the standard analysis settings and are not intended as a reduced-runtime demonstration.

## Reproducibility

The pipeline uses Bayesian MCMC sampling. Convergence should be assessed using effective sample size and split R-hat diagnostics.

For detailed model definitions, prior distributions, and uncertainty propagation, see Supplementary Note 4 and the parent README.

The example datasets are intended to demonstrate the analysis workflow. Reproduction of the full manuscript results requires the corresponding complete datasets and analysis settings.
