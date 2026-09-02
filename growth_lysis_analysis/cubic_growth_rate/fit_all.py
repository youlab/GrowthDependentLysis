"""
Fit a cubic polynomial (order 3) to ln_od_ind as a function of time for
each (species, rep) combination, using BayesLinearRegress.

Model: ln_od_ind = beta0 + beta1*t + beta2*t^2 + beta3*t^3 + noise

Uses rescale_data=True so the wide HalfCauchy(100) priors are appropriate
in standardized space. Without rescaling, 5 data points can't constrain the
posterior under such diffuse priors. The model handles all standardization
internally — no manual centering needed.

Growth rate: mu(t) = d(ln_od)/dt (computed in analyze_all.py)
"""
import sys
import os
import csv
import argparse
import numpy as np

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(project_root, "src"))

from BayesLinearRegress import BayesLinearRegress

DEFAULT_DATA = os.path.join(project_root, "data", "OD_barcode_counts_3replicates.txt")

def parse_args():
    parser = argparse.ArgumentParser(description="Fit cubic polynomial models per (species, rep)")
    parser.add_argument("--data", default=DEFAULT_DATA, help="Path to OD data file (tab-separated)")
    parser.add_argument("--name", default="default", help="Run name (results stored under results/<name>/)")
    parser.add_argument("--mcmc-warmup", type=int, default=20000, help="MCMC warmup iterations (default: 4000)")
    parser.add_argument("--mcmc-samples", type=int, default=20000, help="MCMC posterior samples (default: 4000)")
    return parser.parse_args()

def load_data(path):
    """Read tab-separated data, group by (species, rep).

    Returns {(species, rep): (time_array, ln_od_ind_array)}.
    """
    groups = {}
    with open(path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            species = row["name"].strip()
            rep = int(row["rep"])
            t = float(row["time"])
            ln_od = float(row["ln_od_ind"])
            groups.setdefault((species, rep), ([], []))
            groups[(species, rep)][0].append(t)
            groups[(species, rep)][1].append(ln_od)
    return {k: (np.array(ts), np.array(ys)) for k, (ts, ys) in groups.items()}


def build_poly_features(t_raw):
    """Build polynomial design matrix [t, t^2, t^3].

    The intercept column is NOT included — BayesLinearRegress adds it
    via use_intercept=True. Standardization is handled by rescale_data=True.
    """
    return np.column_stack([t_raw, t_raw**2, t_raw**3])


def main():
    args = parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    models_dir = os.path.join(script_dir, "results", args.name, "models")

    os.makedirs(models_dir, exist_ok=True)
    data = load_data(args.data)
    keys = sorted(data.keys(), key=lambda k: (int(k[0][1:]), k[1]))

    print(f"Found {len(keys)} (species, rep) combinations")
    print(f"Saving models to {models_dir}/")

    for species, rep in keys:
        t_raw, y = data[(species, rep)]
        order = np.argsort(t_raw)
        t_raw, y = t_raw[order], y[order]

        X = build_poly_features(t_raw)

        print(f"  {species} rep{rep} ({len(y)} obs) ... ", end="", flush=True)

        model = BayesLinearRegress(
            X, y,
            use_intercept=True,
            rescale_data=True,
            MCMC_warmup=args.mcmc_warmup,
            MCMC_samples=args.mcmc_samples,
        )
        model.fit(rng_seed=42)

        save_path = os.path.join(models_dir, f"{species}_rep{rep}.pkl")
        model.save(save_path)
        print("done")

    print(f"\nAll {len(keys)} models saved to {models_dir}/")


if __name__ == "__main__":
    main()
