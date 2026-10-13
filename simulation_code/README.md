# Ecological simulation notebooks and demos

This directory contains three full-analysis Jupyter notebooks and three reduced, executable demonstration notebooks based on the same simulation functions.

## Full-analysis notebooks

| File | Purpose |
|---|---|
| `gLV_simulation_Shannon.ipynb` | Community diversity under growth-dependent lysis (Supplementary Note 1; Fig. S15). |
| `gLV_simulation_debris.ipynb` | Biomass residue effects on inferred lysis (Supplementary Note 3; Fig. S9). |
| `gLV_simulation_net_rate_fitting.ipynb` | Accuracy of net-rate inference from sampled biomass trajectories (Supplementary Note 2; Fig. S7). |

## Demonstration notebooks

The `demo/` directory contains `demo_Shannon.ipynb`, `demo_debris.ipynb`, and `demo_net_rate_fitting.ipynb`. Each runs a **reduced configuration with 8 populations and 2 random seeds**. The Shannon demo tests four parameter combinations; the other two test two combinations each. They are **functional demonstrations only**, not reproductions of the full manuscript figures or numerical estimates.

### Requirements

Python with Jupyter, NumPy, SciPy, Matplotlib, and scikit-learn. The full-analysis notebooks may additionally require pandas and seaborn. Install in a suitable Python environment and launch Jupyter:

```bash
python -m pip install jupyter numpy scipy matplotlib scikit-learn pandas seaborn
jupyter lab
```

### How to run

Open each `demo/demo_*.ipynb` notebook and select **Run All Cells**, with the working directory set to `demo/`. Each notebook generates an `.npz` output in the current working directory:

| Demo | Generated file |
|---|---|
| Shannon | `demo_shannon_results.npz` |
| Debris | `simulation_results_2seeds_with_debris.npz` |
| Net-rate fitting | `simulation_results_2seeds.npz` |

Reference `.npz` outputs and corresponding numerical summaries are included in `demo/expected_output/`. The demos were executed successfully in a Python environment with the dependencies listed above. Observed notebook execution times were approximately 2 seconds each in that test environment; execution time may vary across machines and does not include environment installation.

### Notes

- Run full-analysis notebooks in a separate working directory to avoid overwriting generated data.
- Stochastic simulation results may vary across software environments.
- The full notebooks retain the original manuscript-scale settings and are not intended as quick demos.
- These notebooks generate their own simulation inputs; no external experimental FASTQ or clinical data are needed.
