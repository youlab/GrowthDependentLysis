"""
Model 4: Stacked Regression (G_hat -> L) via Composition Sampling.

Minimal edits to compute multiple timepoints' G/L summaries while keeping
the stacked regression (composition sampling) unchanged and still using the
single timepoint specified by --time for the regression itself.

For each posterior draw s from Models 1 and 2:
  1. G_hat^(s) = Model 1 prediction of growth rate (at OD data's ln_ptr)
  2. dy/dt^(s) = Model 2 time derivative of ln_od (at OD data's time)
  3. L^(s) = G_hat^(s) - dy/dt^(s)
  4. Fit BLR: L ~ N(gamma_0 + gamma_1 * G_hat, sigma_eta^2)
  5. Draw one sample (gamma^(s), sigma_eta^(s)) from this posterior

Run per rep, pooling all 57 strains (n = 57 * 5 = 285 observations per fit).
"""
import sys
import os
import csv
import pickle
import argparse
import time as timer
import numpy as np

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(project_root, "src"))

import jax
jax.config.update("jax_platform_name", "cpu")
import jax.random as random
import numpyro

numpyro.set_host_device_count(1)

from BayesLinearRegress import BayesLinearRegress

DEFAULT_DATA = os.path.join(project_root, "data", "OD_barcode_counts_3replicates.txt")

WARMUP = 2000

def parse_args():
    parser = argparse.ArgumentParser(description="Fit Model 4 stacked regression via composition sampling")
    parser.add_argument("--data", default=DEFAULT_DATA, help="Path to OD data file (tab-separated)")
    parser.add_argument("--name", default="default", help="Run name (results stored under results/<name>/)")
    parser.add_argument("--model1-run", default="default", help="Model 1 run name (strain_results/results/<name>/)")
    parser.add_argument("--model2-run", default="default", help="Model 2 run name (cubic_growth_rate/results/<name>/)")
    parser.add_argument("--time", type=float, default=6.0, help="Time point to filter OD data (default: 6.0)")
    return parser.parse_args()


def load_od_data(path):
    """Load OD data, group by (species, rep)."""
    groups = {}
    with open(path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            species = row["name"].strip()
            rep = int(row["rep"])
            key = (species, rep)
            groups.setdefault(key, {"time": [], "ln_ptr": [], "ln_od_ind": []})
            groups[key]["time"].append(float(row["time"]))
            groups[key]["ln_ptr"].append(float(row["ln_ptr"]))
            groups[key]["ln_od_ind"].append(float(row["ln_od_ind"]))
    return {k: {c: np.array(v) for c, v in d.items()} for k, d in groups.items()}


def predict_G_hat(model1, ln_ptr):
    """Predict growth rate using Model 1 through rescaled-space path.

    Args:
        model1: fitted BayesLinearRegress (use_intercept=True, rescale_data=True)
        ln_ptr: (n,) ln_ptr values from OD data

    Returns:
        (S, n) predictions of growth rate in original space
    """
    n = len(ln_ptr)
    X_full = np.column_stack([np.ones(n), ln_ptr])
    X_r = (X_full - model1.X_mean) / model1.X_std
    y_r = X_r @ model1.beta.T  # (n, S)
    return (y_r * model1.y_std + model1.y_mean).T  # (S, n)


def compute_dydt(model2, time):
    """Compute d(ln_od)/dt from Model 2 posterior samples.

    Args:
        model2: fitted BayesLinearRegress for cubic poly (rescale_data=True)
        time: (n,) raw time values

    Returns:
        (S, n) growth rates d(ln_od)/dt
    """
    X_std = model2.X_std
    b1 = model2.beta[:, 1]
    b2 = model2.beta[:, 2]
    b3 = model2.beta[:, 3]
    return model2.y_std * (
        b1[:, None] / X_std[1]
        + 2 * b2[:, None] * time[None, :] / X_std[2]
        + 3 * b3[:, None] * time[None, :]**2 / X_std[3]
    )
    
    
### ADD (helper writers) ###
def write_table_tsv(path, rows, fieldnames):
    """Write rows (list of dicts) to TSV."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
### END ADD ###


def main():
    args = parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "results", args.name)
    model1_dir = os.path.join(project_root, "strain_results", "results", args.model1_run, "models")
    model2_dir = os.path.join(project_root, "cubic_growth_rate", "results", args.model2_run, "models")

    os.makedirs(output_dir, exist_ok=True)

    od_data = load_od_data(args.data)
    species_list = sorted(
        set(k[0] for k in od_data.keys()),
        key=lambda s: int(s[1:]),
    )
    print(f"Loaded OD data for {len(species_list)} species")

    # list of species to exclude from the regression / Bayesian refit
    data_filename = os.path.basename(args.data).lower()
    if "resistant" in data_filename:
        EXCLUDE_SPECIES = {"k30"}
    else:
        EXCLUDE_SPECIES = set()

    print(f"Data file: {data_filename}")
    print(f"Excluded species: {sorted(EXCLUDE_SPECIES)}")
    # filter species_list BEFORE loading models
    species_list = [s for s in species_list if s not in EXCLUDE_SPECIES]
    print(f"Using {len(species_list)} species for modeling.")

    # Load all Model 1s (one per strain)
    model1s = {}
    for species in species_list:
        model1s[species] = BayesLinearRegress.load(
            os.path.join(model1_dir, f"{species}.pkl")
        )
    S1 = model1s[species_list[0]].beta.shape[0]
    print(f"Model 1: {len(model1s)} models, {S1} samples each")
    
    
    # === MODIFIED: prepare scatter rows collector (to be written at end) ===
    scatter_rows = []  # will collect dicts with keys: rep,species,time,G_mean,L_mean,G_lo,G_hi,L_lo,L_hi
    # === END MODIFIED ===

    for rep in [1, 2, 3]:
        print(f"\n{'='*50}")
        print(f"  Rep {rep}")
        print(f"{'='*50}")

        # Pre-compute G_hat and dy/dt for all strains (for composition sampling)
        obs_info = []
        G_hat_parts = []
        dydt_parts = []

        for species in species_list:
            d = od_data[(species, rep)]
            order = np.argsort(d["time"])
            time_vals_all = d["time"][order]  # === MODIFIED: all times available for this (species,rep)
            ln_ptr_vals_all = d["ln_ptr"][order]

            # === MODIFIED: compute G_hat and dydt for ALL timepoints for scatter output ===
            if len(time_vals_all) == 0:
                print(f"    WARNING: {species} rep {rep} has no time points, skipping")
                continue

            # Model 1: predict growth rate at ALL ln_ptr values
            G_hat_all = predict_G_hat(model1s[species], ln_ptr_vals_all)  # (S, n_all)

            # Model 2: compute d(ln_od)/dt at ALL times
            model2_all = BayesLinearRegress.load(
                os.path.join(model2_dir, f"{species}_rep{rep}.pkl")
            )
            dydt_all = compute_dydt(model2_all, time_vals_all)  # (S, n_all)

            # For each timepoint at this species/rep compute summary stats and append to scatter_rows
            # note: we compute G_mean/G_lo/G_hi and L_mean/L_lo/L_hi for each (species,rep,time)
            S_all = G_hat_all.shape[0]
            for j, t_j in enumerate(time_vals_all):
                # samples from model1 and model2 for this timepoint
                G_samples = G_hat_all[:, j]       # shape (S_g,)
                d_samples = dydt_all[:, j]        # shape (S_d,)

                # === ALIGN SAMPLE SIZES (take last min samples) ===
                S_g = G_samples.shape[0]
                S_d = d_samples.shape[0]
                S_common = min(S_g, S_d)
                if S_common == 0:
                    # nothing to do
                    continue

                # take the most recent S_common draws from each (mimics earlier behavior of using [-S:])
                if S_g != S_common:
                    G_samples = G_samples[-S_common:]
                if S_d != S_common:
                    d_samples = d_samples[-S_common:]
                # ===================================================

                # now compute L and summaries
                G_samples = np.maximum(G_samples, 0.0)
                L_samples = np.maximum(G_samples - d_samples, 0.0)

                G_mean = float(np.mean(G_samples))
                G_lo = float(np.percentile(G_samples, 2.5))
                G_hi = float(np.percentile(G_samples, 97.5))

                L_mean = float(np.mean(L_samples))
                L_lo = float(np.percentile(L_samples, 2.5))
                L_hi = float(np.percentile(L_samples, 97.5))

                scatter_rows.append({
                    "rep": rep,
                    "species": species,
                    "time": float(t_j),
                    "G_mean": G_mean,
                    "L_mean": L_mean,
                    "G_lo": G_lo,
                    "G_hi": G_hi,
                    "L_lo": L_lo,
                    "L_hi": L_hi,
                })

            # Now continue original flow: select the single timepoint args.time for composition/regression
            # Filter to selected time point (exact match)
            t_mask = time_vals_all == args.time
            if not np.any(t_mask):
                print(f"    WARNING: {species} rep {rep} has no t={args.time} data, skipping composition sampling for this species")
                continue
            time_vals = time_vals_all[t_mask]
            ln_ptr_vals = ln_ptr_vals_all[t_mask]

            # Model 1: predict growth rate at OD data's ln_ptr (for the regression/composition sampling)
            G_hat = predict_G_hat(model1s[species], ln_ptr_vals)

            # Model 2: compute d(ln_od)/dt at OD data's time (for the regression/composition sampling)
            model2 = model2_all 
            dydt = compute_dydt(model2, time_vals)

            S = min(S1, dydt.shape[0])
            G_hat_parts.append(G_hat[-S:])
            dydt_parts.append(dydt[-S:])

            for t in time_vals:
                obs_info.append({
                    "species": species,
                    "time": float(t),
                    "rep": rep,
                })

        # Pool across strains: (S, n_obs) for composition/regression (unchanged)
        if len(G_hat_parts) == 0:
            print("  No observations for composition sampling in this rep, skipping.")
            continue

        
        G_hat_pooled = np.concatenate(G_hat_parts, axis=1)
        dydt_pooled = np.concatenate(dydt_parts, axis=1)

        # truncate G first
        G_hat_pooled = np.maximum(G_hat_pooled, 0.0)
        # then define L consistently from truncated G
        L_pooled = np.maximum(G_hat_pooled - dydt_pooled, 0.0)     

        # Filter out observations where mean G_hat across draws is negative (unchanged)
        mean_G_hat = G_hat_pooled.mean(axis=0)
        pos_mask = mean_G_hat > 0
        n_dropped = (~pos_mask).sum()
        if n_dropped > 0:
            print(f"  Dropping {n_dropped} obs with negative mean G_hat")
            G_hat_pooled = G_hat_pooled[:, pos_mask]
            dydt_pooled = dydt_pooled[:, pos_mask]
            L_pooled = L_pooled[:, pos_mask]
            obs_info = [obs_info[i] for i in range(len(obs_info)) if pos_mask[i]]

        S = G_hat_pooled.shape[0]
        n_obs = G_hat_pooled.shape[1]
        print(f"  {n_obs} observations, {S} composition samples")

        # Create BLR from mean data (for normalization stats)
        G_mean = G_hat_pooled.mean(axis=0)  # (n_obs,)
        L_mean = L_pooled.mean(axis=0)      # (n_obs,)
        blr = BayesLinearRegress(
            X=G_mean.reshape(-1, 1), y=L_mean,
            use_intercept=True, rescale_data=True,
            MCMC_warmup=WARMUP,
        )
        kernel = blr.create_kernel()

        # Composition sampling (unchanged)
        gamma_0_out = np.zeros(S)
        gamma_1_out = np.zeros(S)
        sigma_out = np.zeros(S)

        t0 = timer.time()
        print(f"  Composition sampling (warmup={WARMUP})...", flush=True)

        for s in range(S):
            if s % 100 == 0:
                elapsed = timer.time() - t0
                rate = s / elapsed if elapsed > 0 and s > 0 else 0
                eta = (S - s) / rate / 60 if rate > 0 else 0
                print(
                    f"    s={s}/{S}  ({elapsed:.0f}s elapsed, "
                    f"~{eta:.1f}min remaining)",
                    flush=True,
                )

            samp = blr.run_single(
                kernel,
                X_raw=G_hat_pooled[s].reshape(-1, 1),
                y_raw=L_pooled[s],
                rng_key=random.PRNGKey(s),
            )

            # Convert rescaled-space samples to original-space gamma/sigma
            beta_r = samp["beta"][0]      # (2,) rescaled
            sigma_r = float(samp["sigma"][0])
            gamma_1_out[s] = blr.y_std * beta_r[1] / blr.X_std[1]
            gamma_0_out[s] = blr.y_mean - gamma_1_out[s] * blr.X_mean[1]
            sigma_out[s] = sigma_r * blr.y_std

        elapsed = timer.time() - t0
        print(f"  Done in {elapsed:.0f}s ({elapsed/60:.1f}min)")

        # Save results (unchanged)
        results = {
            "G_hat": G_hat_pooled,       # (S, n_obs)
            "L": L_pooled,               # (S, n_obs)
            "dydt": dydt_pooled,         # (S, n_obs)
            "gamma_0": gamma_0_out,      # (S,) intercept
            "gamma_1": gamma_1_out,      # (S,) slope
            "sigma_eta": sigma_out,      # (S,) noise
            "obs_info": obs_info,        # list of dicts
        }

        path = os.path.join(output_dir, f"model4_rep{rep}.pkl")
        with open(path, "wb") as f:
            pickle.dump(results, f)

        g1_lo = np.percentile(gamma_1_out, 2.5)
        g1_hi = np.percentile(gamma_1_out, 97.5)
        print(f"  gamma_1 (slope): {gamma_1_out.mean():.4f} "
              f"[{g1_lo:.4f}, {g1_hi:.4f}]")
        print(f"  gamma_0 (intercept): {gamma_0_out.mean():.4f}")
        print(f"  sigma_eta: {sigma_out.mean():.4f}")
        print(f"  Saved -> {path}")

    # === MODIFIED: write scatter files with all collected timepoints ===
    if len(scatter_rows) > 0:
        scatter_tsv = os.path.join(output_dir, "scatter.txt")
        scatter_csv = os.path.join(output_dir, "scatter.csv")
        fieldnames = ["rep", "species", "time", "G_mean", "L_mean", "G_lo", "G_hi", "L_lo", "L_hi"]
        write_table_tsv(scatter_tsv, scatter_rows, fieldnames)
        write_csv(scatter_csv, scatter_rows, fieldnames)
        print(f"\nWrote scatter summaries for all timepoints -> {scatter_tsv} and {scatter_csv}")
    else:
        print("\nNo scatter rows were collected; no scatter files written.")
    # === END MODIFIED ===

    print("\nAll reps complete!")


if __name__ == "__main__":
    main()
    