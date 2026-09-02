# GrowthLysis

Estimated runtime for a dataset (growth rates v.s. ln(PTR), ln(OD) v.s. time, and lysis rate v.s. growth rate): 4-5 Hours

Bayesian regression pipeline for microbial growth data. Predicts growth rate from ln(PTR) across 57 bacterial strains using three stacked models.

## Setup

Create and activate a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate        # Linux/macOS
.venv/Scripts/activate           # Windows
```

Install dependencies (dev mode):

```bash
pip install -e .
```

## Pipeline Overview

```
Model 1 (Linear)                Model 2 (Cubic Polynomial)
growth_rate ~ ln_ptr            ln_od ~ t + t^2 + t^3
strain_results/                 cubic_growth_rate/
        |                               |
        |  G_hat = predicted             |  dy/dt = d(ln_od)/dt
        |  growth rate                   |  at time t
        v                               v
       +-----------------------------------+
       |  Model 4 (Stacked Regression)     |
       |  L = G_hat - dy/dt               |
       |  L ~ gamma_0 + gamma_1 * G_hat   |
       |  stacked_model/                   |
       +-----------------------------------+
```

- **Model 1**: Per-strain Bayesian linear regression of growth rate on ln(PTR).
- **Model 2**: Per-(strain, rep) cubic polynomial fit of ln(OD) over time; growth rate is the analytical derivative.
- **Model 4**: Stacked regression combining Models 1 and 2 via composition sampling.

## Usage

All scripts support `--name` to organize results into named runs under `results/<name>/`. The default run name is `"default"`.

### Model 1: Linear Regression

```bash
# Fit
python strain_results/fit_all_strains.py [--data PATH] [--name RUN] \
    [--mcmc-warmup N] [--mcmc-samples N]

# Analyze
python strain_results/analyze_strains.py [--data PATH] [--name RUN]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--data` | `data/first_data-lnPTR-growth_rate-plot_57strains.csv` | Input CSV |
| `--name` | `default` | Run name |
| `--mcmc-warmup` | `5000` | MCMC warmup iterations |
| `--mcmc-samples` | `5000` | MCMC posterior samples |

### Model 2: Cubic Polynomial

```bash
# Fit
python cubic_growth_rate/fit_all.py [--data PATH] [--name RUN] \
    [--mcmc-warmup N] [--mcmc-samples N]

# Analyze
python cubic_growth_rate/analyze_all.py [--data PATH] [--name RUN]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--data` | `data/OD_barcode_counts_3replicates.txt` | Input tab-separated file |
| `--name` | `default` | Run name |
| `--mcmc-warmup` | `10000` | MCMC warmup iterations |
| `--mcmc-samples` | `10000` | MCMC posterior samples |

### Model 4: Stacked Regression

```bash
# Fit (depends on Models 1 and 2)
python stacked_model/fit_model4.py [--data PATH] [--name RUN] \
    [--model1-run RUN] [--model2-run RUN] [--time T]

# Analyze
python stacked_model/analyze_model4.py [--name RUN]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--data` | `data/OD_barcode_counts_3replicates.txt` | Input OD data |
| `--name` | `default` | Run name for Model 4 output |
| `--model1-run` | `default` | Which Model 1 run to use |
| `--model2-run` | `default` | Which Model 2 run to use |
| `--time` | `6.0` | Time point to filter OD data |

### Example: Full Pipeline with Custom Data

```bash
python strain_results/fit_all_strains.py --data data/new_data.csv --name experiment1
python strain_results/analyze_strains.py --data data/new_data.csv --name experiment1

python cubic_growth_rate/fit_all.py --data data/new_od.txt --name experiment1
python cubic_growth_rate/analyze_all.py --data data/new_od.txt --name experiment1

python stacked_model/fit_model4.py --data data/new_od.txt --name experiment1 \
    --model1-run experiment1 --model2-run experiment1 --time 6.0
python stacked_model/analyze_model4.py --name experiment1
```

## Output Structure

Each run produces results under `results/<name>/` within its model directory:

```
strain_results/results/<name>/
    models/          # k1.pkl, k3.pkl, ...
    plots/           # k1.png, k3.png, ...
    summary.csv

cubic_growth_rate/results/<name>/
    models/          # k1_rep1.pkl, k1_rep2.pkl, ...
    plots/           # k1.png, k3.png, ...
    summary.csv

stacked_model/results/<name>/
    model4_rep1.pkl, model4_rep2.pkl, model4_rep3.pkl
    plots/           # model4_rep1.png, ...
    summary.csv
```

## Convergence Diagnostics

After fitting, check that the MCMC chains have converged:

```bash
python diagnostics/check_convergence.py [--name RUN]
```

This loads all saved posterior samples from Models 1, 2, and 4, computes two standard diagnostics for every parameter, and reports any warnings. Results are saved to `diagnostics/results/<name>/convergence.csv`.

| Argument | Default | Description |
|----------|---------|-------------|
| `--name` | `default` | Run name |
| `--fix` | off | Auto-refit failing Model 1/2 models with doubled iterations |
| `--max-retries` | `3` | Max refit attempts (iterations multiply by 2x, 4x, 8x) |
| `--data1` | `data/first_data-lnPTR-growth_rate-plot_57strains.csv` | Model 1 data path (used with `--fix`) |
| `--data2` | `data/OD_barcode_counts_3replicates.txt` | Model 2 data path (used with `--fix`) |

### Diagnostics Explained

**Effective Sample Size (ESS)** measures how many independent samples your MCMC chain is equivalent to. MCMC draws are autocorrelated (each sample resembles the previous one), so 2,000 draws may carry the information content of far fewer independent samples. ESS < 100 means posterior summaries, especially tail quantiles like the 2.5th and 97.5th percentiles used for 95% credible intervals, may be unreliable. The fix is to increase `--mcmc-samples`.

**Split R-hat** checks whether the chain has reached its stationary distribution. It splits the chain in half and compares the two halves: if they look like they came from the same distribution, R-hat is close to 1.0. R-hat > 1.05 indicates the chain may still be drifting and has not converged. The fix is to increase `--mcmc-warmup` (giving the sampler more time to adapt) and/or `--mcmc-samples`.

### Thresholds

| Diagnostic | Passing | Warning | Meaning |
|------------|---------|---------|---------|
| ESS | >= 100 | < 100 | Too few effective samples for reliable posterior summaries |
| R-hat | < 1.05 | >= 1.05 | Chain has not converged to the stationary distribution |

### Fixing Convergence Warnings

The easiest way to fix warnings is the `--fix` flag, which automatically refits only the failing Model 1/2 models with progressively more MCMC iterations:

```bash
python diagnostics/check_convergence.py --fix [--max-retries 3]
```

This doubles warmup and samples on each retry (2x, 4x, 8x) and re-checks after each round, stopping once all models pass or retries are exhausted. Model 4 is skipped since it depends on Models 1 and 2.

Alternatively, you can manually re-fit all models with more iterations:

```bash
# Model 1: increase from default 5000
python strain_results/fit_all_strains.py --name refit --mcmc-warmup 10000 --mcmc-samples 10000

# Model 2: increase from default 10000 (most common source of warnings due to n=5)
python cubic_growth_rate/fit_all.py --name refit --mcmc-warmup 20000 --mcmc-samples 20000

# Re-check
python diagnostics/check_convergence.py --name refit
```

Model 2 (cubic polynomial) is the most likely to show convergence issues because it fits 4 parameters to only 5 data points per strain-replicate, leaving the posterior poorly constrained for some combinations.
