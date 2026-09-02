"""
Fit a BayesLinearRegress model per strain from the 57-strain CSV
and save each model under strain_results/results/<name>/models/.
"""
import sys
import os
import csv
import argparse
import numpy as np

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(project_root, "src"))

from BayesLinearRegress import BayesLinearRegress

DEFAULT_DATA = os.path.join(project_root, "data", "first_data-lnPTR-growth_rate-plot_57strains.csv")


def parse_args():
    parser = argparse.ArgumentParser(description="Fit linear models per strain")
    parser.add_argument("--data", default=DEFAULT_DATA, help="Path to CSV data file")
    parser.add_argument("--name", default="default", help="Run name (results stored under results/<name>/)")
    parser.add_argument("--mcmc-warmup", type=int, default=5000, help="MCMC warmup iterations (default: 3000)")
    parser.add_argument("--mcmc-samples", type=int, default=5000, help="MCMC posterior samples (default: 3000)")
    return parser.parse_args()


def load_data(path):
    """Read CSV and group by strain. Returns {strain: (X, y)}."""
    strains = {}
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            strain = row["strain"].strip().strip('"')
            gr = float(row["growth_rate"])
            lnptr = float(row["ln_ptr"])
            strains.setdefault(strain, ([], []))
            strains[strain][0].append(lnptr)
            strains[strain][1].append(gr)
    return {s: (np.array(x).reshape(-1, 1), np.array(y))
            for s, (x, y) in strains.items()}


def main():
    args = parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    models_dir = os.path.join(script_dir, "results", args.name, "models")

    os.makedirs(models_dir, exist_ok=True)
    data = load_data(args.data)
    strain_names = sorted(data.keys(), key=lambda s: int(s[1:]))

    print(f"Found {len(strain_names)} strains")
    print(f"Saving models to {models_dir}/")

    for strain in strain_names:
        X, y = data[strain]
        n = len(y)
        print(f"  {strain} ({n} obs) ... ", end="", flush=True)

        model = BayesLinearRegress(
            X, y,
            use_intercept=True,
            rescale_data=True,
            MCMC_warmup=args.mcmc_warmup,
            MCMC_samples=args.mcmc_samples,
        )
        model.fit(rng_seed=42)

        save_path = os.path.join(models_dir, f"{strain}.pkl")
        model.save(save_path)
        print("done")

    print(f"\nAll models saved to {models_dir}/")


if __name__ == "__main__":
    main()
