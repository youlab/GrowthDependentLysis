"""
MCMC Convergence Diagnostics for All Models.

Loads saved posterior samples from Models 1, 2, and 4, computes split R-hat
and effective sample size (ESS) for each parameter, and reports warnings
for R-hat > 1.05 or ESS < 100.

With --fix, automatically refits failing Model 1 and Model 2 models with
doubled MCMC iterations (up to --max-retries times) until they pass.

Usage:
    python diagnostics/check_convergence.py [--name RUN_NAME]
    python diagnostics/check_convergence.py --fix [--max-retries 3]
"""
import sys
import os
import csv
import pickle
import argparse
import numpy as np

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(project_root, "src"))

from numpyro.diagnostics import summary as numpyro_summary

RHAT_THRESHOLD = 1.05
ESS_THRESHOLD = 100

DEFAULT_DATA1 = os.path.join(project_root, "data", "first_data-lnPTR-growth_rate-plot_57strains.csv")
DEFAULT_DATA2 = os.path.join(project_root, "data", "OD_barcode_counts_3replicates.txt")


def parse_args():
    parser = argparse.ArgumentParser(description="Check MCMC convergence diagnostics")
    parser.add_argument("--name", default="default", help="Run name")
    parser.add_argument("--fix", action="store_true",
                        help="Auto-refit failing Model 1/2 models with more iterations")
    parser.add_argument("--max-retries", type=int, default=3,
                        help="Max refit attempts per model (default: 3)")
    parser.add_argument("--data1", default=DEFAULT_DATA1,
                        help="Model 1 data path")
    parser.add_argument("--data2", default=DEFAULT_DATA2,
                        help="Model 2 data path")
    return parser.parse_args()


def diagnose_samples(samples_dict, label, skip=None):
    """Compute ESS and split R-hat for a dict of posterior sample arrays.

    Args:
        samples_dict: {param_name: array of shape (S,) or (S, K)}
        label: string label for this model (e.g. "Model1/k1")
        skip: set of (param_name, index_or_None) tuples to exclude from
              warnings (still recorded in CSV but flag is cleared)

    Returns:
        list of dicts with diagnostic info, one per parameter (or per
        parameter component for vector params)
    """
    skip = skip or set()
    diag = numpyro_summary(samples_dict, prob=0.9, group_by_chain=False)
    rows = []
    for param, stats in diag.items():
        n_eff = np.atleast_1d(stats["n_eff"])
        r_hat = np.atleast_1d(stats["r_hat"])
        for i in range(len(n_eff)):
            suffix = f"[{i}]" if len(n_eff) > 1 else ""
            # Check if this param/index should be skipped
            is_skipped = (param, i) in skip or (param, None) in skip
            flag = ""
            if not is_skipped:
                if r_hat[i] > RHAT_THRESHOLD:
                    flag += "R-hat>{:.2f} ".format(RHAT_THRESHOLD)
                if n_eff[i] < ESS_THRESHOLD:
                    flag += "ESS<{} ".format(ESS_THRESHOLD)
            rows.append({
                "model": label,
                "param": f"{param}{suffix}",
                "ess": f"{n_eff[i]:.0f}",
                "r_hat": f"{r_hat[i]:.4f}",
                "flag": flag.strip(),
            })
    return rows


def load_blr_samples(filepath):
    """Load BayesLinearRegress pickle and return posterior samples dict.

    Also returns a set of (param, index) tuples to skip in diagnostics.
    When rescale_data=True and use_intercept=True, the intercept column
    (all ones) gets standardized to near-zero, making beta[0] unidentifiable.
    tau inherits this problem since it governs the shared prior scale.
    These are known artifacts, not real convergence failures.
    """
    with open(filepath, "rb") as f:
        state = pickle.load(f)

    skip = set()
    if state.get("use_intercept") and state.get("rescale_data"):
        skip.add(("beta", 0))
        skip.add(("tau", None))

    return {
        "beta": state["beta"],
        "sigma": state["sigma"],
        "tau": state["tau"],
    }, skip


def load_model4_samples(filepath):
    """Load Model 4 pickle and return posterior samples dict."""
    with open(filepath, "rb") as f:
        r = pickle.load(f)
    return {
        "gamma_0": r["gamma_0"],
        "gamma_1": r["gamma_1"],
        "sigma_eta": r["sigma_eta"],
    }


def load_model1_data(path):
    """Read Model 1 CSV and group by strain. Returns {strain: (X, y)}.

    X is (n, 1) ln_ptr, y is (n,) growth_rate.
    """
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


def load_model2_data(path):
    """Read Model 2 tab-separated data, group by (species, rep).

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
    """Build polynomial design matrix [t, t^2, t^3]."""
    return np.column_stack([t_raw, t_raw**2, t_raw**3])


def refit_model(filepath, X, y, multiplier=2):
    """Refit a BayesLinearRegress model with increased MCMC iterations.

    Loads the existing pickle to read current settings, creates a new model
    with warmup and samples multiplied by `multiplier`, fits it, and
    overwrites the saved file.

    Returns the new posterior samples dict for re-diagnosis.
    """
    from BayesLinearRegress import BayesLinearRegress

    with open(filepath, "rb") as f:
        state = pickle.load(f)

    new_warmup = state["MCMC_warmup"] * multiplier
    new_samples = state["MCMC_samples"] * multiplier

    model = BayesLinearRegress(
        X, y,
        use_intercept=state["use_intercept"],
        rescale_data=state["rescale_data"],
        MCMC_warmup=new_warmup,
        MCMC_samples=new_samples,
        param_hyperprior=state["param_hyperprior"],
        measurement_noise_hyperprior=state["measurement_noise_hyperprior"],
    )
    model.fit(rng_seed=42)
    model.save(filepath)

    return {
        "beta": model.beta,
        "sigma": model.sigma,
        "tau": model.tau,
    }


def run_diagnostics(model1_dir, model2_dir, model4_dir, only_models=None):
    """Run diagnostics on all models and return (all_rows, failing_models).

    Args:
        model1_dir, model2_dir, model4_dir: paths to model directories
        only_models: if set, only diagnose these model labels (for re-checking)

    Returns:
        all_rows: list of diagnostic row dicts
        failing_models: set of model labels with at least one non-skipped flag
    """
    all_rows = []
    failing_models = set()

    # --- Model 1: Linear per strain ---
    if os.path.isdir(model1_dir):
        files = sorted(
            [f for f in os.listdir(model1_dir) if f.endswith(".pkl")],
            key=lambda f: int(f.replace(".pkl", "")[1:]),
        )
        if only_models is None:
            print(f"Model 1 (Linear): {len(files)} models from {model1_dir}")
        for fname in files:
            strain = fname.replace(".pkl", "")
            label = f"Model1/{strain}"
            if only_models is not None and label not in only_models:
                continue
            samples, skip = load_blr_samples(os.path.join(model1_dir, fname))
            rows = diagnose_samples(samples, label, skip)
            all_rows.extend(rows)
            if any(r["flag"] for r in rows):
                failing_models.add(label)
    elif only_models is None:
        print(f"Model 1: directory not found ({model1_dir})")

    # --- Model 2: Cubic per (strain, rep) ---
    if os.path.isdir(model2_dir):
        files = sorted(
            [f for f in os.listdir(model2_dir) if f.endswith(".pkl")],
            key=lambda f: (
                int(f.split("_rep")[0][1:]),
                int(f.split("_rep")[1].replace(".pkl", "")),
            ),
        )
        if only_models is None:
            print(f"Model 2 (Cubic):  {len(files)} models from {model2_dir}")
        for fname in files:
            label_name = fname.replace(".pkl", "")
            label = f"Model2/{label_name}"
            if only_models is not None and label not in only_models:
                continue
            samples, skip = load_blr_samples(os.path.join(model2_dir, fname))
            rows = diagnose_samples(samples, label, skip)
            all_rows.extend(rows)
            if any(r["flag"] for r in rows):
                failing_models.add(label)
    elif only_models is None:
        print(f"Model 2: directory not found ({model2_dir})")

    # --- Model 4: Stacked per rep ---
    if os.path.isdir(model4_dir):
        rep_files = sorted(
            [f for f in os.listdir(model4_dir) if f.startswith("model4_rep") and f.endswith(".pkl")]
        )
        if only_models is None:
            print(f"Model 4 (Stacked): {len(rep_files)} reps from {model4_dir}")
        for fname in rep_files:
            rep = fname.replace("model4_rep", "").replace(".pkl", "")
            label = f"Model4/rep{rep}"
            if only_models is not None and label not in only_models:
                continue
            samples = load_model4_samples(os.path.join(model4_dir, fname))
            rows = diagnose_samples(samples, label)
            all_rows.extend(rows)
            if any(r["flag"] for r in rows):
                failing_models.add(label)
    elif only_models is None:
        print(f"Model 4: directory not found ({model4_dir})")

    return all_rows, failing_models


def print_summary(all_rows):
    """Print the diagnostics summary table."""
    print(f"\n{'='*70}")
    print(f"  Convergence Diagnostics Summary")
    print(f"{'='*70}")
    print(f"  Thresholds: R-hat < {RHAT_THRESHOLD}, ESS > {ESS_THRESHOLD}")
    print(f"  Total parameters checked: {len(all_rows)}")

    flagged = [r for r in all_rows if r["flag"]]
    n_warnings = len(flagged)

    if n_warnings == 0:
        print(f"\n  All parameters passed. No convergence issues detected.")
    else:
        print(f"\n  {n_warnings} WARNING(S):")
        print(f"  {'Model':<25} {'Param':<12} {'ESS':>6} {'R-hat':>8}  Flag")
        print(f"  {'-'*25} {'-'*12} {'-'*6} {'-'*8}  {'-'*20}")
        for r in flagged:
            print(f"  {r['model']:<25} {r['param']:<12} {r['ess']:>6} {r['r_hat']:>8}  {r['flag']}")

    # Print per-model summary stats
    print(f"\n  Per-model R-hat and ESS ranges:")
    for prefix, label in [("Model1/", "Model 1 (Linear)"),
                          ("Model2/", "Model 2 (Cubic)"),
                          ("Model4/", "Model 4 (Stacked)")]:
        model_rows = [r for r in all_rows if r["model"].startswith(prefix)]
        if not model_rows:
            continue
        rhats = [float(r["r_hat"]) for r in model_rows]
        esses = [float(r["ess"]) for r in model_rows]
        n_flag = sum(1 for r in model_rows if r["flag"])
        print(f"    {label}:")
        print(f"      R-hat: [{min(rhats):.4f}, {max(rhats):.4f}]  "
              f"ESS: [{min(esses):.0f}, {max(esses):.0f}]  "
              f"Warnings: {n_flag}/{len(model_rows)}")

    return n_warnings


def main():
    args = parse_args()

    model1_dir = os.path.join(project_root, "strain_results", "results", args.name, "models")
    model2_dir = os.path.join(project_root, "cubic_growth_rate", "results", args.name, "models")
    model4_dir = os.path.join(project_root, "stacked_model", "results", args.name)

    # --- Initial diagnostics pass ---
    all_rows, failing_models = run_diagnostics(model1_dir, model2_dir, model4_dir)
    print_summary(all_rows)

    # Save CSV
    output_dir = os.path.join(project_root, "diagnostics", "results", args.name)
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, "convergence.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["model", "param", "ess", "r_hat", "flag"])
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\n  Full results -> {csv_path}")

    # --- Auto-fix loop ---
    if not args.fix:
        return

    # Filter to only Model 1 and Model 2 failures (skip Model 4)
    fixable = {m for m in failing_models if m.startswith("Model1/") or m.startswith("Model2/")}
    if not fixable:
        print("\n  --fix: No fixable failures (Model 1/2). Nothing to do.")
        return

    # Load raw data lazily (only if needed)
    model1_data = None
    model2_data = None

    retry = 0
    while retry < args.max_retries and fixable:
        multiplier = 2 ** (retry + 1)
        print(f"\n{'='*70}")
        print(f"  Auto-fix attempt {retry + 1}/{args.max_retries} "
              f"({len(fixable)} model(s), {multiplier}x iterations)")
        print(f"{'='*70}")

        for label in sorted(fixable):
            if label.startswith("Model1/"):
                strain = label.split("/")[1]
                filepath = os.path.join(model1_dir, f"{strain}.pkl")
                if model1_data is None:
                    print(f"  Loading Model 1 data from {args.data1}")
                    model1_data = load_model1_data(args.data1)
                X, y = model1_data[strain]
                print(f"  Refitting {label} ({multiplier}x) ... ", end="", flush=True)
                refit_model(filepath, X, y, multiplier=multiplier)
                print("done")

            elif label.startswith("Model2/"):
                name = label.split("/")[1]  # e.g. "k9_rep1"
                species, rep_str = name.rsplit("_rep", 1)
                rep = int(rep_str)
                filepath = os.path.join(model2_dir, f"{name}.pkl")
                if model2_data is None:
                    print(f"  Loading Model 2 data from {args.data2}")
                    model2_data = load_model2_data(args.data2)
                t_raw, ln_od = model2_data[(species, rep)]
                order = np.argsort(t_raw)
                t_raw, ln_od = t_raw[order], ln_od[order]
                X = build_poly_features(t_raw)
                print(f"  Refitting {label} ({multiplier}x) ... ", end="", flush=True)
                refit_model(filepath, X, ln_od, multiplier=multiplier)
                print("done")

        # Re-diagnose only the models we just refitted
        print(f"\n  Re-checking refitted models...")
        recheck_rows, still_failing = run_diagnostics(
            model1_dir, model2_dir, model4_dir, only_models=fixable
        )

        # Update all_rows: remove old entries for refitted models, add new ones
        all_rows = [r for r in all_rows if r["model"] not in fixable]
        all_rows.extend(recheck_rows)
        print_summary(all_rows)

        # Save updated CSV
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["model", "param", "ess", "r_hat", "flag"])
            writer.writeheader()
            writer.writerows(all_rows)
        print(f"\n  Updated results -> {csv_path}")

        # Update fixable set to only models that still fail
        fixable = {m for m in still_failing if m.startswith("Model1/") or m.startswith("Model2/")}
        retry += 1

    if fixable:
        print(f"\n  {len(fixable)} model(s) still failing after {args.max_retries} retries:")
        for m in sorted(fixable):
            print(f"    {m}")
    else:
        print(f"\n  All models now pass convergence checks.")


if __name__ == "__main__":
    main()
