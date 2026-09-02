"""
Load fitted per-strain models, compute summary statistics (R^2, 95%
credible intervals for slope & intercept), save a summary CSV, and
produce a scatter + regression-line plot for each strain.
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

DEFAULT_DATA = os.path.join(project_root, "data", "first_data-lnPTR-growth_rate-plot_57strains.csv")


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze fitted per-strain linear models")
    parser.add_argument("--data", default=DEFAULT_DATA, help="Path to CSV data file")
    parser.add_argument("--name", default="default", help="Run name (results under results/<name>/)")
    return parser.parse_args()

CI_LEVEL = 0.95  # 95 % credible interval
CI_LO = (1 - CI_LEVEL) / 2        # 0.025
CI_HI = 1 - (1 - CI_LEVEL) / 2    # 0.975


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


def original_space_params(model):
    """Convert every MCMC sample of (beta0, beta1) to original-space
    slope and intercept arrays."""
    slopes = (model.y_std * model.beta[:, 1]) / model.X_std[1]
    intercepts = (model.y_mean
                  + model.y_std * model.beta[:, 0]
                  - model.y_std * model.beta[:, 1] * model.X_mean[1] / model.X_std[1])
    return slopes, intercepts


def credible_interval(samples):
    """Return (lo, hi, excludes_zero) for the 99.9 % equal-tailed CI."""
    lo = np.percentile(samples, 100 * CI_LO)
    hi = np.percentile(samples, 100 * CI_HI)
    excludes_zero = (lo > 0) or (hi < 0)
    return lo, hi, excludes_zero


def compute_r_squared(model, y):
    """Bayesian R^2 using the posterior-mean predictions."""
    beta_mean = model.beta.mean(axis=0)
    y_hat_rescaled = model.X_rescaled @ beta_mean
    y_hat = y_hat_rescaled * model.y_std + model.y_mean
    ss_res = np.sum((y - y_hat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def predict_grid(model, x_grid_raw):
    """Compute posterior mean-line predictions on a 1-D grid in original space.

    Goes through the model's own rescaled-space path so the CI band
    is not inflated by the unidentifiable intercept column.

    Returns:
        y_samples: (n_mcmc, n_grid) predictions in original y space
    """
    # Build design matrix the same way the model does
    if model.use_intercept:
        X_grid = np.column_stack((np.ones(len(x_grid_raw)), x_grid_raw))
    else:
        X_grid = x_grid_raw
    if model.rescale_data:
        X_grid = (X_grid - model.X_mean) / model.X_std

    # (n_grid, n_features) @ (n_features, n_mcmc) -> (n_grid, n_mcmc)
    y_rescaled = X_grid @ model.beta.T
    if model.rescale_data:
        y_samples = y_rescaled * model.y_std + model.y_mean
    else:
        y_samples = y_rescaled
    return y_samples.T  # (n_mcmc, n_grid)


def compute_plot_data(model, X, y, slopes):
    """Compute all grid quantities needed for plotting and data export.

    Returns dict with scatter data, dense grid, CI bands, and scalars.
    """
    x_raw = X[:, 0]
    x_grid = np.linspace(x_raw.min() - 0.1, x_raw.max() + 0.1, 200).reshape(-1, 1)

    # Posterior predictions via rescaled-space path (stable CI band)
    mu = predict_grid(model, x_grid)  # (n_mcmc, n_grid) — mean line only
    y_mean_line = mu.mean(axis=0)
    mu_lo = np.percentile(mu, 100 * CI_LO, axis=0)
    mu_hi = np.percentile(mu, 100 * CI_HI, axis=0)

    # Posterior predictive (add observation noise)
    sigma_orig = model.sigma * model.y_std if model.rescale_data else model.sigma
    y_pred = mu + np.random.randn(*mu.shape) * sigma_orig[:, np.newaxis]
    pp_lo = np.percentile(y_pred, 100 * CI_LO, axis=0)
    pp_hi = np.percentile(y_pred, 100 * CI_HI, axis=0)

    return {
        "x_scatter": x_raw,
        "y_scatter": y,
        "x_grid": x_grid[:, 0],
        "y_mean": y_mean_line,
        "mu_lo": mu_lo,
        "mu_hi": mu_hi,
        "pp_lo": pp_lo,
        "pp_hi": pp_hi,
    }


def plot_strain(strain, pd, slopes, r2, save_path):
    """Scatter plot with posterior mean line, mean CI, and predictive CI."""
    slope_ci = credible_interval(slopes)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(pd["x_scatter"], pd["y_scatter"], s=20, alpha=0.6, color="black", zorder=3)
    ax.plot(pd["x_grid"], pd["y_mean"], color="blue", linewidth=1.5, label="Posterior mean")
    ax.fill_between(pd["x_grid"], pd["mu_lo"], pd["mu_hi"], alpha=0.3, color="blue", label="95% CI (mean)")
    ax.fill_between(pd["x_grid"], pd["pp_lo"], pd["pp_hi"], alpha=0.12, color="blue", label="95% CI (pred)")
    ax.set_xlabel("ln(PTR)")
    ax.set_ylabel("Growth rate")
    ax.set_title(f"{strain}   R\u00b2={r2:.3f}")

    sig_s = "*" if slope_ci[2] else ""
    stats_text = (
        f"slope = {slopes.mean():.3f} [{slope_ci[0]:.3f}, {slope_ci[1]:.3f}]{sig_s}\n"
        f"n = {len(pd['y_scatter'])}"
    )
    ax.text(0.03, 0.97, stats_text, transform=ax.transAxes,
            fontsize=7, verticalalignment="top", fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    ax.legend(fontsize=6.5, loc="lower right")
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
    strain_names = sorted(data.keys(), key=lambda s: int(s[1:]))

    print(f"Analyzing {len(strain_names)} strains  (CI level = {CI_LEVEL*100:.1f}%)")

    summary_rows = []
    scatter_rows = []
    curves_rows = []

    for strain in strain_names:
        X, y = data[strain]
        n = len(y)

        model_path = os.path.join(models_dir, f"{strain}.pkl")
        model = BayesLinearRegress.load(model_path)
        # Restore training data for R^2 computation
        if model.use_intercept:
            model.X = np.column_stack((np.ones(len(X)), X))
        else:
            model.X = X
        model.y = y
        model.X_rescaled = (model.X - model.X_mean) / model.X_std if model.rescale_data else model.X
        model.y_rescaled = (model.y - model.y_mean) / model.y_std if model.rescale_data else model.y

        slopes, intercepts = original_space_params(model)
        r2 = compute_r_squared(model, y)
        slope_lo, slope_hi, slope_sig = credible_interval(slopes)
        intercept_lo, intercept_hi, _ = credible_interval(intercepts)


        sig_label = "*" if slope_sig else ""
        print(f"  {strain:>4s}  n={n:>2d}  R²={r2:.3f}  "
              f"slope={slopes.mean():.3f} [{slope_lo:.3f},{slope_hi:.3f}]{sig_label}")

        pd = compute_plot_data(model, X, y, slopes)

        plot_path = os.path.join(plots_dir, f"{strain}.pdf")
        plot_strain(strain, pd, slopes, r2, plot_path)

        # Save NPZ
        np.savez(
            os.path.join(plot_data_dir, f"{strain}.npz"),
            x_scatter=pd["x_scatter"], y_scatter=pd["y_scatter"],
            x_grid=pd["x_grid"], y_mean=pd["y_mean"],
            mu_lo=pd["mu_lo"], mu_hi=pd["mu_hi"],
            pp_lo=pd["pp_lo"], pp_hi=pd["pp_hi"],
            slope_mean=slopes.mean(), slope_ci_lo=slope_lo, slope_ci_hi=slope_hi,
            intercept_mean = intercepts.mean(), intercept_ci_lo=intercept_lo, intercept_ci_hi=intercept_hi,
            r_squared=r2, n_obs=n,
        )

        # Accumulate CSV rows
        for xi, yi in zip(pd["x_scatter"], pd["y_scatter"]):
            scatter_rows.append({"strain": strain, "x": f"{xi:.6f}", "y": f"{yi:.6f}"})
        for i in range(len(pd["x_grid"])):
            curves_rows.append({
                "strain": strain,
                "x": f"{pd['x_grid'][i]:.6f}",
                "y_mean": f"{pd['y_mean'][i]:.6f}",
                "mu_lo": f"{pd['mu_lo'][i]:.6f}",
                "mu_hi": f"{pd['mu_hi'][i]:.6f}",
                "pp_lo": f"{pd['pp_lo'][i]:.6f}",
                "pp_hi": f"{pd['pp_hi'][i]:.6f}",
            })

        summary_rows.append({
            "strain": strain,
            "n_obs": n,
            "r_squared": f"{r2:.6f}",
            "slope_mean": f"{slopes.mean():.6f}",
            "slope_ci_lo": f"{slope_lo:.6f}",
            "slope_ci_hi": f"{slope_hi:.6f}",
            "slope_excludes_zero": slope_sig,
            "intercept_mean": f"{intercepts.mean():.6f}",
            "intercept_ci_lo": f"{intercept_lo:.6f}",
            "intercept_ci_hi": f"{intercept_hi:.6f}",
            "model_file": f"models/{strain}.pkl",
        })

    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "strain", "n_obs", "r_squared",
            "slope_mean", "slope_ci_lo", "slope_ci_hi", "slope_excludes_zero",
            "intercept_mean", "intercept_ci_lo", "intercept_ci_hi",
            "model_file",
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

    write_csv_and_txt("scatter", ["strain", "x", "y"], scatter_rows)
    write_csv_and_txt("curves", [
        "strain", "x", "y_mean", "mu_lo", "mu_hi", "pp_lo", "pp_hi",
    ], curves_rows)

    print(f"\nSummary   -> {summary_path}")
    print(f"Plots     -> {plots_dir}/")
    print(f"Plot data -> {plot_data_dir}/")


if __name__ == "__main__":
    main()
