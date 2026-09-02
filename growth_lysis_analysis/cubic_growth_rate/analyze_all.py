"""
Analyze fitted cubic polynomial models to extract growth rates.

Growth rate mu(t) = d(ln_od_ind)/dt is computed analytically from
the posterior samples of the polynomial coefficients.

Since rescale_data=True, predictions and derivatives go through
the rescaled-space path. The constant (intercept) column drops out
of the derivative, so the known rescaling bug doesn't affect growth rates.

Design matrix is [1, t, t^2, t^3] (raw time, no manual centering).
The model standardizes internally.

Derivative formula (original space):
    mu(t) = y_std * (beta1/X_std[1] + 2*beta2*t/X_std[2] + 3*beta3*t^2/X_std[3])

Produces:
  - Per-species plots (ln OD fit + growth rate, 3 reps overlaid)
  - Summary CSV with growth rate at each time point
"""
import sys
import os
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(project_root, "src"))

from BayesLinearRegress import BayesLinearRegress

DEFAULT_DATA = os.path.join(project_root, "data", "OD_barcode_counts_3replicates.txt")


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze fitted cubic polynomial models")
    parser.add_argument("--data", default=DEFAULT_DATA, help="Path to OD data file (tab-separated)")
    parser.add_argument("--name", default="default", help="Run name (results under results/<name>/)")
    return parser.parse_args()

TIME_POINTS = np.array([4.0, 5.0, 6.0, 8.0, 10.0])

CI_LEVEL = 0.95
CI_LO = (1 - CI_LEVEL) / 2
CI_HI = 1 - CI_LO

REP_COLORS = {1: "tab:blue", 2: "tab:orange", 3: "tab:green"}


def load_data(path):
    """Read tab-separated data, group by (species, rep)."""
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


def predict_curve(model, t_eval):
    """Posterior samples of ln_od at given times via rescaled-space path.

    Goes through the model's own rescaled coordinates so the CI band
    is not inflated by the unidentifiable intercept column.

    Args:
        model: fitted BayesLinearRegress with rescale_data=True
        t_eval: (n_t,) time in original coordinates

    Returns:
        (n_mcmc, n_t) posterior samples of ln_od_ind in original space
    """
    X_full = np.column_stack([np.ones_like(t_eval), t_eval, t_eval**2, t_eval**3])
    X_r = (X_full - model.X_mean) / model.X_std
    y_r = (X_r @ model.beta.T).T  # (n_mcmc, n_t)
    return y_r * model.y_std + model.y_mean


def compute_growth_rate(model, t_eval):
    """Posterior samples of mu(t) = d(ln_od)/dt.

    Analytically differentiates the cubic polynomial through rescaled space.
    The intercept column drops out (derivative of a constant is zero).

    mu(t) = y_std * sum_j [beta_j / X_std[j]] * dX_j/dt

    where X = [1, t, t^2, t^3] and dX/dt = [0, 1, 2*t, 3*t^2]

    Args:
        model: fitted BayesLinearRegress with rescale_data=True
        t_eval: (n_t,) time in original coordinates

    Returns:
        (n_mcmc, n_t) posterior samples of growth rate (h^-1)
    """
    X_std = model.X_std  # [X_std_intercept, X_std_t, X_std_t2, X_std_t3]

    # beta[:, 0] is intercept (drops out), beta[:, 1:4] are polynomial coeffs
    b1 = model.beta[:, 1]  # (n_mcmc,)
    b2 = model.beta[:, 2]
    b3 = model.beta[:, 3]

    mu = model.y_std * (
        b1[:, None] / X_std[1]
        + 2 * b2[:, None] * t_eval[None, :] / X_std[2]
        + 3 * b3[:, None] * t_eval[None, :]**2 / X_std[3]
    )
    return mu


def compute_r_squared(model, t_data, y_data):
    """R^2 using posterior-mean predictions."""
    y_hat = predict_curve(model, t_data).mean(axis=0)
    ss_res = np.sum((y_data - y_hat) ** 2)
    ss_tot = np.sum((y_data - y_data.mean()) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def compute_rep_curves(model, t_dense):
    """Compute fit curves and growth rate curves for one rep.

    Returns dict with all curve arrays needed for plotting and export.
    """
    # Cubic fit curves
    y_curve = predict_curve(model, t_dense)   # (n_mcmc, 200)
    y_mean = y_curve.mean(axis=0)
    y_lo = np.percentile(y_curve, 100 * CI_LO, axis=0)
    y_hi = np.percentile(y_curve, 100 * CI_HI, axis=0)

    # Posterior predictive (with noise)
    sigma_orig = model.sigma * model.y_std
    y_pred = y_curve + np.random.randn(*y_curve.shape) * sigma_orig[:, None]
    pp_lo = np.percentile(y_pred, 100 * CI_LO, axis=0)
    pp_hi = np.percentile(y_pred, 100 * CI_HI, axis=0)

    # Growth rate curves
    mu = compute_growth_rate(model, t_dense)
    mu_mean = mu.mean(axis=0)
    mu_lo = np.percentile(mu, 100 * CI_LO, axis=0)
    mu_hi = np.percentile(mu, 100 * CI_HI, axis=0)

    return {
        "y_mean": y_mean, "y_lo": y_lo, "y_hi": y_hi,
        "pp_lo": pp_lo, "pp_hi": pp_hi,
        "mu_mean": mu_mean, "mu_lo": mu_lo, "mu_hi": mu_hi,
    }


def plot_species(species, reps_data, reps_curves, t_dense, save_path):
    """Two-panel plot: ln(OD) fits and growth rate, 3 reps overlaid."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    for rep in sorted(reps_data.keys()):
        t_data, y_data = reps_data[rep]
        c = reps_curves[rep]
        color = REP_COLORS[rep]

        # --- Left: cubic fit ---
        ax1.scatter(t_data, y_data, s=40, color=color, zorder=3, alpha=0.8)
        ax1.plot(t_dense, c["y_mean"], color=color, linewidth=1.5, label=f"Rep {rep}")
        ax1.fill_between(t_dense, c["y_lo"], c["y_hi"], alpha=0.2, color=color)
        ax1.fill_between(t_dense, c["pp_lo"], c["pp_hi"], alpha=0.08, color=color)

        # --- Right: growth rate ---
        ax2.plot(t_dense, c["mu_mean"], color=color, linewidth=1.5, label=f"Rep {rep}")
        ax2.fill_between(t_dense, c["mu_lo"], c["mu_hi"], alpha=0.15, color=color)

    ax1.set_xlabel("Time (h)")
    ax1.set_ylabel("ln(OD_ind)")
    ax1.set_title(f"{species}: Cubic Fit")
    ax1.legend(fontsize=7)

    ax2.axhline(0, color="gray", linestyle="--", linewidth=0.5)
    ax2.set_xlabel("Time (h)")
    ax2.set_ylabel("d(ln OD)/dt  [1/h]")
    ax2.set_title(f"{species}: Growth Rate")
    ax2.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def main():
    args = parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    run_dir = os.path.join(script_dir, "results", args.name)
    models_dir = os.path.join(run_dir, "models")
    plots_dir = os.path.join(run_dir, "plots")
    plot_data_dir = os.path.join(run_dir, "plot_data")
    summary_path = os.path.join(run_dir, "summary.csv")

    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(plot_data_dir, exist_ok=True)
    data = load_data(args.data)

    species_list = sorted(set(k[0] for k in data.keys()),
                          key=lambda s: int(s[1:]))

    print(f"Analyzing {len(species_list)} species  (CI level = {CI_LEVEL*100:.0f}%)")

    summary_rows = []
    scatter_rows = []
    fit_curves_rows = []
    gr_curves_rows = []

    t_dense = np.linspace(3.5, 10.5, 200)

    for species in species_list:
        reps_data = {}
        reps_curves = {}
        npz_dict = {"t_dense": t_dense}

        for rep in [1, 2, 3]:
            key = (species, rep)
            t_raw, y = data[key]
            order = np.argsort(t_raw)
            t_raw, y = t_raw[order], y[order]
            reps_data[rep] = (t_raw, y)

            model_path = os.path.join(models_dir, f"{species}_rep{rep}.pkl")
            model = BayesLinearRegress.load(model_path)

            r2 = compute_r_squared(model, t_raw, y)

            # Compute curves for plotting and export
            curves = compute_rep_curves(model, t_dense)
            reps_curves[rep] = curves

            # Store in NPZ dict
            npz_dict[f"t_data_rep{rep}"] = t_raw
            npz_dict[f"y_data_rep{rep}"] = y
            for k, v in curves.items():
                npz_dict[f"{k}_rep{rep}"] = v

            # Accumulate scatter CSV rows
            for ti, yi in zip(t_raw, y):
                scatter_rows.append({
                    "species": species, "rep": rep,
                    "t": f"{ti:.6f}", "ln_od": f"{yi:.6f}",
                })

            # Accumulate fit curves CSV rows
            for i in range(len(t_dense)):
                fit_curves_rows.append({
                    "species": species, "rep": rep,
                    "t": f"{t_dense[i]:.6f}",
                    "y_mean": f"{curves['y_mean'][i]:.6f}",
                    "y_lo": f"{curves['y_lo'][i]:.6f}",
                    "y_hi": f"{curves['y_hi'][i]:.6f}",
                    "pp_lo": f"{curves['pp_lo'][i]:.6f}",
                    "pp_hi": f"{curves['pp_hi'][i]:.6f}",
                })

            # Accumulate growth rate curves CSV rows
            for i in range(len(t_dense)):
                gr_curves_rows.append({
                    "species": species, "rep": rep,
                    "t": f"{t_dense[i]:.6f}",
                    "mu_mean": f"{curves['mu_mean'][i]:.6f}",
                    "mu_lo": f"{curves['mu_lo'][i]:.6f}",
                    "mu_hi": f"{curves['mu_hi'][i]:.6f}",
                })

            # Growth rate at each observed time point (for summary)
            mu = compute_growth_rate(model, TIME_POINTS)
            for i, t in enumerate(TIME_POINTS):
                mu_lo = np.percentile(mu[:, i], 100 * CI_LO)
                mu_hi = np.percentile(mu[:, i], 100 * CI_HI)
                excludes_zero = (mu_lo > 0) or (mu_hi < 0)
                summary_rows.append({
                    "species": species,
                    "rep": rep,
                    "time": f"{t:.0f}",
                    "mu_mean": f"{mu[:, i].mean():.6f}",
                    "mu_ci_lo": f"{mu_lo:.6f}",
                    "mu_ci_hi": f"{mu_hi:.6f}",
                    "mu_excludes_zero": excludes_zero,
                    "r_squared": f"{r2:.6f}",
                    "sigma_mean": f"{model.sigma.mean():.6f}",
                    "model_file": f"models/{species}_rep{rep}.pkl",
                })

            def s(flag):
                return "*" if flag else " "

            mu_at_5 = mu[:, 1]  # t=5
            mu_lo5 = np.percentile(mu_at_5, 100 * CI_LO)
            mu_hi5 = np.percentile(mu_at_5, 100 * CI_HI)
            ex5 = (mu_lo5 > 0) or (mu_hi5 < 0)
            print(f"  {species:>4s} rep{rep}  R2={r2:.3f}  "
                  f"mu(t=5)={mu_at_5.mean():.4f}[{mu_lo5:.4f},{mu_hi5:.4f}]{s(ex5)}  "
                  f"sigma={model.sigma.mean():.4f}")

        # Save NPZ for this species (all reps)
        np.savez(os.path.join(plot_data_dir, f"{species}.npz"), **npz_dict)

        plot_path = os.path.join(plots_dir, f"{species}.pdf")
        plot_species(species, reps_data, reps_curves, t_dense, plot_path)

    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "species", "rep", "time",
            "mu_mean", "mu_ci_lo", "mu_ci_hi", "mu_excludes_zero",
            "r_squared", "sigma_mean", "model_file",
        ])
        writer.writeheader()
        writer.writerows(summary_rows)

    # Write plot-data CSVs and tab-separated TXT copies
    def write_csv_and_txt(basename, fieldnames, rows):
        csv_path = os.path.join(plot_data_dir, basename + ".csv")
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        txt_path = os.path.join(plot_data_dir, basename + ".txt")
        with open(txt_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    write_csv_and_txt("scatter", ["species", "rep", "t", "ln_od"], scatter_rows)
    write_csv_and_txt("fit_curves", [
        "species", "rep", "t", "y_mean", "y_lo", "y_hi", "pp_lo", "pp_hi",
    ], fit_curves_rows)
    write_csv_and_txt("growth_rate_curves", [
        "species", "rep", "t", "mu_mean", "mu_lo", "mu_hi",
    ], gr_curves_rows)

    print(f"\nSummary   -> {summary_path}")
    print(f"Plots     -> {plots_dir}/")
    print(f"Plot data -> {plot_data_dir}/")


if __name__ == "__main__":
    main()
