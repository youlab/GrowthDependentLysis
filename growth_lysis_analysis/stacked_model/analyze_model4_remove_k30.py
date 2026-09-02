"""
Analyze Model 4 (Stacked Regression) results from composition sampling.

Produces per-rep scatter plots of G_hat vs L with:
  - Error bars on both axes (upstream posterior uncertainty)
  - Regression line with mean CI and predictive CI
  - Summary CSV
"""
import os
import csv
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze Model 4 stacked regression results")
    parser.add_argument("--name", default="default", help="Run name (results under results/<name>/)")
    return parser.parse_args()

# 68% CI for data-point error bars (1-sigma)
EB_LEVEL = 0.68
EB_LO = (1 - EB_LEVEL) / 2
EB_HI = 1 - EB_LO

# 95% CI for regression bands and significance tests
CI_LEVEL = 0.95
CI_LO = (1 - CI_LEVEL) / 2
CI_HI = 1 - CI_LO


def main():
    args = parse_args()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    run_dir = os.path.join(script_dir, "results", args.name)
    results_dir = run_dir
    plots_dir = os.path.join(run_dir, "plots")
    plot_data_dir = os.path.join(run_dir, "plot_data")
    summary_path = os.path.join(run_dir, "summary.csv")


    # --- START: species exclusion rule ---
    # If input/run name contains "resistant", exclude k30; otherwise exclude none
    if "resistant" in args.name.lower():
        EXCLUDE_SPECIES = {"k30"}
    else:
        EXCLUDE_SPECIES = set()
    # --- END: species exclusion rule ---


    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(plot_data_dir, exist_ok=True)

    summary_rows = []
    scatter_rows = []
    reg_curves_rows = []

    for rep in [1, 2, 3]:
        path = os.path.join(results_dir, f"model4_rep{rep}.pkl")
        with open(path, "rb") as f:
            r = pickle.load(f)

        G = r["G_hat"]         # (S, n_obs)
        L = r["L"]             # (S, n_obs)
        g0 = r["gamma_0"]     # (S,)
        g1 = r["gamma_1"]     # (S,)
        sig = r["sigma_eta"]  # (S,)
        obs_info = r["obs_info"]
        
        
        # filter out species in EXCLUDE_SPECIES
        mask = np.array([o.get("species") not in EXCLUDE_SPECIES for o in obs_info], dtype=bool)
        if not mask.any():
            print(f"Rep {rep}: all observations excluded by EXCLUDE_SPECIES={EXCLUDE_SPECIES}; skipping.")
            continue
        G = G[:, mask]
        L = L[:, mask]
        obs_info = [o for o, m in zip(obs_info, mask) if m]

        # Data point means and 68% CI error bars
        G_mean = G.mean(axis=0)
        G_lo = np.percentile(G, 100 * EB_LO, axis=0)
        G_hi = np.percentile(G, 100 * EB_HI, axis=0)
        L_mean = L.mean(axis=0)
        L_lo = np.percentile(L, 100 * EB_LO, axis=0)
        L_hi = np.percentile(L, 100 * EB_HI, axis=0)

        # R-squared using posterior means
        L_hat = g0.mean() + g1.mean() * G_mean
        ss_res = np.sum((L_mean - L_hat) ** 2)
        ss_tot = np.sum((L_mean - L_mean.mean()) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

        # Regression line on a dense grid
        pad = 0.05 * (G_mean.max() - G_mean.min() + 0.1)
        x_grid = np.linspace(G_mean.min() - pad, G_mean.max() + pad, 200)

        # Mean-only CI (epistemic uncertainty in gamma)
        lines = g0[:, None] + g1[:, None] * x_grid[None, :]  # (S, 200)
        line_mean = lines.mean(axis=0)
        line_lo = np.percentile(lines, 100 * CI_LO, axis=0)
        line_hi = np.percentile(lines, 100 * CI_HI, axis=0)

        # Predictive CI (add observation noise sigma_eta)
        pred = lines + np.random.randn(*lines.shape) * sig[:, None]
        pred_lo = np.percentile(pred, 100 * CI_LO, axis=0)
        pred_hi = np.percentile(pred, 100 * CI_HI, axis=0)

        # --- Plot ---
        fig, ax = plt.subplots(figsize=(8, 5))

        # Assign colors per strain
        strains = [o["species"] for o in obs_info]
        unique_strains = sorted(set(strains), key=lambda s: int(s[1:]))
        strain_to_idx = {s: i for i, s in enumerate(unique_strains)}
        cmap = plt.colormaps.get_cmap("turbo").resampled(len(unique_strains))
        colors = [cmap(strain_to_idx[s]) for s in strains]

        # Error bars on data points, colored by strain
        xerr = np.array([G_mean - G_lo, G_hi - G_mean])
        yerr = np.array([L_mean - L_lo, L_hi - L_mean])
        for i in range(len(G_mean)):
            ax.errorbar(
                G_mean[i], L_mean[i], xerr=[[xerr[0, i]], [xerr[1, i]]],
                yerr=[[yerr[0, i]], [yerr[1, i]]],
                fmt="o", ms=3, color=colors[i], alpha=0.6,
                elinewidth=0.5, capsize=1, zorder=3,
            )

        # Colorbar for strains
        sm = plt.cm.ScalarMappable(
            cmap=cmap,
            norm=plt.Normalize(0, len(unique_strains) - 1),
        )
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.02, aspect=30)
        n_ticks = min(10, len(unique_strains))
        tick_idx = np.linspace(0, len(unique_strains) - 1, n_ticks, dtype=int)
        cbar.set_ticks(tick_idx)
        cbar.set_ticklabels([unique_strains[i] for i in tick_idx])
        cbar.ax.tick_params(labelsize=6)

        # Regression line + bands
        ax.plot(x_grid, line_mean, "b-", lw=1.5, label="Posterior mean")
        ax.fill_between(
            x_grid, line_lo, line_hi, alpha=0.3, color="blue",
            label="95% CI (mean)",
        )
        ax.fill_between(
            x_grid, pred_lo, pred_hi, alpha=0.12, color="blue",
            label="95% CI (pred)",
        )

        ax.set_xlabel("G_hat (predicted growth rate)")
        ax.set_ylabel("L = G_hat - dy/dt")
        ax.set_title(f"Model 4: Stacked Regression (Rep {rep})")

        # Stats annotation
        s_lo = np.percentile(g1, 100 * CI_LO)
        s_hi = np.percentile(g1, 100 * CI_HI)
        
        # intercept CI
        g0_lo = np.percentile(g0, 100 * CI_LO)
        g0_hi = np.percentile(g0, 100 * CI_HI)

        # sigma CI
        sig_lo = np.percentile(sig, 100 * CI_LO)
        sig_hi = np.percentile(sig, 100 * CI_HI)

        ex = (s_lo > 0) or (s_hi < 0)
        sig_s = "*" if ex else ""

        stats = (
            f"slope = {g1.mean():.4f} [{s_lo:.4f}, {s_hi:.4f}]{sig_s}\n"
            f"intercept = {g0.mean():.4f}\n"
            f"sigma = {sig.mean():.4f}\n"
            f"R² = {r2:.4f}\n"
            f"n = {G.shape[1]}"
        )
        ax.text(
            0.03, 0.97, stats, transform=ax.transAxes,
            fontsize=7, va="top", fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
        )

        ax.legend(fontsize=7, loc="lower right")
        plt.tight_layout()

        plot_path = os.path.join(plots_dir, f"model4_rep{rep}.pdf")
        fig.savefig(plot_path)
        plt.close(fig)

        # Save NPZ
        strain_labels = np.array(strains)
        np.savez(
            os.path.join(plot_data_dir, f"model4_rep{rep}.npz"),
            G_mean=G_mean, L_mean=L_mean,
            G_lo=G_lo, G_hi=G_hi, L_lo=L_lo, L_hi=L_hi,
            strain_labels=strain_labels,
            x_grid=x_grid, line_mean=line_mean,
            line_lo=line_lo, line_hi=line_hi,
            pred_lo=pred_lo, pred_hi=pred_hi,
            gamma_1_mean=g1.mean(), gamma_1_ci_lo=s_lo, gamma_1_ci_hi=s_hi,
            gamma_0_mean=g0.mean(), sigma_mean=sig.mean(), r_squared=r2,
        )

        # Accumulate scatter CSV rows
        for i in range(len(G_mean)):
            scatter_rows.append({
                "rep": rep, "species": strains[i],
                "G_mean": f"{G_mean[i]:.6f}", "L_mean": f"{L_mean[i]:.6f}",
                "G_lo": f"{G_lo[i]:.6f}", "G_hi": f"{G_hi[i]:.6f}",
                "L_lo": f"{L_lo[i]:.6f}", "L_hi": f"{L_hi[i]:.6f}",
            })

        # Accumulate regression curves CSV rows
        for i in range(len(x_grid)):
            # === MODIFIED: include regression stats in each row ===
            reg_curves_rows.append({
                "rep": rep,
                "x": f"{x_grid[i]:.6f}",
                "line_mean": f"{line_mean[i]:.6f}",
                "line_lo": f"{line_lo[i]:.6f}",
                "line_hi": f"{line_hi[i]:.6f}",
                "pred_lo": f"{pred_lo[i]:.6f}",
                "pred_hi": f"{pred_hi[i]:.6f}",
                "gamma_1_mean": f"{g1.mean():.6f}",
                "gamma_1_ci_lo": f"{s_lo:.6f}",
                "gamma_1_ci_hi": f"{s_hi:.6f}",
                "gamma_0_mean": f"{g0.mean():.6f}",
                "gamma_0_ci_lo": f"{g0_lo:.6f}",
                "gamma_0_ci_hi": f"{g0_hi:.6f}",
                "sigma_mean": f"{sig.mean():.6f}",
                "sigma_ci_lo": f"{sig_lo:.6f}",
                "sigma_ci_hi": f"{sig_hi:.6f}",
                "r_squared": f"{r2:.6f}",
            })
            # === END MODIFIED ===

        print(
            f"Rep {rep}: slope={g1.mean():.4f} "
            f"[{s_lo:.4f}, {s_hi:.4f}]{sig_s}  R²={r2:.4f}"
        )
        print(f"  Plot -> {plot_path}")

        summary_rows.append({
            "rep": rep,
            "n_obs": G.shape[1],
            "r_squared": f"{r2:.6f}",
            "gamma_1_mean": f"{g1.mean():.6f}",
            "gamma_1_ci_lo": f"{s_lo:.6f}",
            "gamma_1_ci_hi": f"{s_hi:.6f}",
            "gamma_1_excludes_zero": ex,
            "gamma_0_mean": f"{g0.mean():.6f}",
            "sigma_eta_mean": f"{sig.mean():.6f}",
        })

    # Save summary CSV
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "rep", "n_obs", "r_squared",
            "gamma_1_mean", "gamma_1_ci_lo", "gamma_1_ci_hi",
            "gamma_1_excludes_zero",
            "gamma_0_mean", "sigma_eta_mean",
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

    write_csv_and_txt("scatter", [
        "rep", "species", "G_mean", "L_mean", "G_lo", "G_hi", "L_lo", "L_hi",
    ], scatter_rows)
    # === MODIFIED: include regression stats in the written fieldnames ===
    write_csv_and_txt("regression_curves", [
        "rep", "x", "line_mean", "line_lo", "line_hi", "pred_lo", "pred_hi",
        "gamma_1_mean", "gamma_1_ci_lo", "gamma_1_ci_hi",
        "gamma_0_mean", "gamma_0_ci_lo", "gamma_0_ci_hi",
        "sigma_mean", "sigma_ci_lo", "sigma_ci_hi",
        "r_squared",
    ], reg_curves_rows)

    print(f"\nSummary   -> {summary_path}")
    print(f"Plots     -> {plots_dir}/")
    print(f"Plot data -> {plot_data_dir}/")


if __name__ == "__main__":
    main()
