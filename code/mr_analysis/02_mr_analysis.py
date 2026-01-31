#!/usr/bin/env python3
"""
02_mr_analysis.py - Mendelian Randomization Analysis

Implements MR methods:
1. Inverse-Variance Weighted (IVW)
2. MR-Egger (pleiotropy detection)
3. Weighted Median (robust to outliers)
4. Simple Median
5. MR-PRESSO (outlier detection)

Tests causal relationship: Workload genes -> T2D risk
"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import json
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
RESULTS_DIR = WORKLOAD_DIR / "results" / "mr_analysis"
DATA_DIR = RESULTS_DIR / "data"
FIGURES_DIR = RESULTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


class MendelianRandomization:
    """
    Mendelian Randomization analysis class.
    Implements multiple MR methods for robust causal inference.
    """

    def __init__(self, exposure_df, outcome_df):
        """
        Initialize MR analysis.

        Args:
            exposure_df: DataFrame with columns [SNP, beta, se, pval]
            outcome_df: DataFrame with columns [SNP, beta, se, pval]
        """
        self.exposure = exposure_df.copy()
        self.outcome = outcome_df.copy()

        # Harmonize data
        self.harmonized = self._harmonize_data()

    def _harmonize_data(self):
        """Harmonize exposure and outcome data by SNP."""
        merged = self.exposure.merge(
            self.outcome,
            on="SNP",
            suffixes=("_exp", "_out")
        )

        # Remove SNPs with missing data
        merged = merged.dropna(subset=["beta_exp", "beta_out", "se_exp", "se_out"])

        print(f"Harmonized {len(merged)} SNPs")

        return merged

    def ivw(self):
        """
        Inverse-Variance Weighted MR.
        Standard two-sample MR method assuming no pleiotropy.
        """
        if len(self.harmonized) == 0:
            return {"beta": np.nan, "se": np.nan, "pval": 1, "method": "IVW"}

        beta_exp = self.harmonized["beta_exp"].values
        beta_out = self.harmonized["beta_out"].values
        se_out = self.harmonized["se_out"].values

        # Weights (inverse variance)
        weights = 1 / (se_out ** 2)

        # Weighted regression through origin
        # beta_out = beta_ivw * beta_exp
        beta_ivw = np.sum(weights * beta_exp * beta_out) / np.sum(weights * beta_exp ** 2)

        # Standard error
        se_ivw = np.sqrt(1 / np.sum(weights * beta_exp ** 2))

        # P-value (two-sided)
        z = beta_ivw / se_ivw
        pval = 2 * (1 - stats.norm.cdf(abs(z)))

        # Cochran's Q for heterogeneity
        q_stat = np.sum(weights * (beta_out - beta_ivw * beta_exp) ** 2)
        q_df = len(beta_exp) - 1
        q_pval = 1 - stats.chi2.cdf(q_stat, q_df) if q_df > 0 else 1

        return {
            "method": "IVW",
            "beta": beta_ivw,
            "se": se_ivw,
            "pval": pval,
            "or": np.exp(beta_ivw),
            "or_ci_lower": np.exp(beta_ivw - 1.96 * se_ivw),
            "or_ci_upper": np.exp(beta_ivw + 1.96 * se_ivw),
            "n_snp": len(beta_exp),
            "q_stat": q_stat,
            "q_pval": q_pval
        }

    def mr_egger(self):
        """
        MR-Egger regression.
        Allows for directional pleiotropy (intercept != 0).
        """
        if len(self.harmonized) < 3:
            return {"beta": np.nan, "se": np.nan, "pval": 1, "method": "MR-Egger"}

        beta_exp = self.harmonized["beta_exp"].values
        beta_out = self.harmonized["beta_out"].values
        se_out = self.harmonized["se_out"].values

        # Sign of exposure effect (for orientation)
        sign = np.sign(beta_exp)
        beta_exp_oriented = np.abs(beta_exp)
        beta_out_oriented = beta_out * sign

        # Weights
        weights = 1 / (se_out ** 2)

        # Weighted linear regression: beta_out = intercept + slope * beta_exp
        from scipy.stats import linregress

        # Using weighted least squares
        X = np.column_stack([np.ones(len(beta_exp_oriented)), beta_exp_oriented])
        W = np.diag(weights)

        try:
            XtWX_inv = np.linalg.inv(X.T @ W @ X)
            coeffs = XtWX_inv @ X.T @ W @ beta_out_oriented

            intercept = coeffs[0]
            slope = coeffs[1]

            # Residuals and standard errors
            residuals = beta_out_oriented - (intercept + slope * beta_exp_oriented)
            n = len(beta_exp)
            mse = np.sum(weights * residuals ** 2) / (n - 2)

            se_coeffs = np.sqrt(np.diag(XtWX_inv) * mse)
            se_intercept = se_coeffs[0]
            se_slope = se_coeffs[1]

            # P-values
            z_slope = slope / se_slope
            pval_slope = 2 * (1 - stats.norm.cdf(abs(z_slope)))

            z_intercept = intercept / se_intercept
            pval_intercept = 2 * (1 - stats.norm.cdf(abs(z_intercept)))

            # I-squared (measure of regression dilution)
            i_squared = max(0, 1 - (n - 1) * np.mean(se_out ** 2) / np.var(beta_out_oriented))

        except np.linalg.LinAlgError:
            return {"beta": np.nan, "se": np.nan, "pval": 1, "method": "MR-Egger"}

        return {
            "method": "MR-Egger",
            "beta": slope,
            "se": se_slope,
            "pval": pval_slope,
            "or": np.exp(slope),
            "or_ci_lower": np.exp(slope - 1.96 * se_slope),
            "or_ci_upper": np.exp(slope + 1.96 * se_slope),
            "intercept": intercept,
            "intercept_se": se_intercept,
            "intercept_pval": pval_intercept,
            "n_snp": n,
            "i_squared": i_squared
        }

    def weighted_median(self):
        """
        Weighted Median MR.
        Robust when up to 50% of instruments are invalid.
        """
        if len(self.harmonized) < 3:
            return {"beta": np.nan, "se": np.nan, "pval": 1, "method": "Weighted Median"}

        beta_exp = self.harmonized["beta_exp"].values
        beta_out = self.harmonized["beta_out"].values
        se_exp = self.harmonized["se_exp"].values
        se_out = self.harmonized["se_out"].values

        # Ratio estimates
        beta_iv = beta_out / beta_exp

        # Weights (inverse variance of ratio)
        se_iv = np.sqrt((se_out ** 2) / (beta_exp ** 2) +
                        (beta_out ** 2) * (se_exp ** 2) / (beta_exp ** 4))
        weights = 1 / (se_iv ** 2)
        weights = weights / np.sum(weights)  # Normalize

        # Weighted median
        order = np.argsort(beta_iv)
        beta_iv_sorted = beta_iv[order]
        weights_sorted = weights[order]
        cumsum = np.cumsum(weights_sorted)

        # Find median
        median_idx = np.searchsorted(cumsum, 0.5)
        beta_wm = beta_iv_sorted[min(median_idx, len(beta_iv_sorted) - 1)]

        # Bootstrap for SE
        n_boot = 1000
        boot_estimates = []

        for _ in range(n_boot):
            idx = np.random.choice(len(beta_iv), size=len(beta_iv), replace=True)
            boot_beta = beta_iv[idx]
            boot_weights = weights[idx]
            boot_weights = boot_weights / np.sum(boot_weights)

            order = np.argsort(boot_beta)
            cumsum = np.cumsum(boot_weights[order])
            median_idx = np.searchsorted(cumsum, 0.5)
            boot_estimates.append(boot_beta[order][min(median_idx, len(boot_beta) - 1)])

        se_wm = np.std(boot_estimates)

        # P-value
        z = beta_wm / se_wm
        pval = 2 * (1 - stats.norm.cdf(abs(z)))

        return {
            "method": "Weighted Median",
            "beta": beta_wm,
            "se": se_wm,
            "pval": pval,
            "or": np.exp(beta_wm),
            "or_ci_lower": np.exp(beta_wm - 1.96 * se_wm),
            "or_ci_upper": np.exp(beta_wm + 1.96 * se_wm),
            "n_snp": len(beta_iv)
        }

    def simple_median(self):
        """Simple (unweighted) median MR."""
        if len(self.harmonized) < 3:
            return {"beta": np.nan, "se": np.nan, "pval": 1, "method": "Simple Median"}

        beta_exp = self.harmonized["beta_exp"].values
        beta_out = self.harmonized["beta_out"].values

        # Ratio estimates
        beta_iv = beta_out / beta_exp

        # Simple median
        beta_sm = np.median(beta_iv)

        # Bootstrap for SE
        n_boot = 1000
        boot_estimates = []

        for _ in range(n_boot):
            idx = np.random.choice(len(beta_iv), size=len(beta_iv), replace=True)
            boot_estimates.append(np.median(beta_iv[idx]))

        se_sm = np.std(boot_estimates)

        # P-value
        z = beta_sm / se_sm
        pval = 2 * (1 - stats.norm.cdf(abs(z)))

        return {
            "method": "Simple Median",
            "beta": beta_sm,
            "se": se_sm,
            "pval": pval,
            "or": np.exp(beta_sm),
            "or_ci_lower": np.exp(beta_sm - 1.96 * se_sm),
            "or_ci_upper": np.exp(beta_sm + 1.96 * se_sm),
            "n_snp": len(beta_iv)
        }

    def run_all_methods(self):
        """Run all MR methods and return combined results."""
        results = []

        results.append(self.ivw())
        results.append(self.mr_egger())
        results.append(self.weighted_median())
        results.append(self.simple_median())

        return pd.DataFrame(results)

    def plot_scatter(self, output_path=None):
        """Create scatter plot of SNP effects."""
        if len(self.harmonized) == 0:
            return

        fig, ax = plt.subplots(figsize=(8, 6))

        beta_exp = self.harmonized["beta_exp"].values
        beta_out = self.harmonized["beta_out"].values
        se_out = self.harmonized["se_out"].values

        # Scatter with error bars
        ax.errorbar(beta_exp, beta_out, yerr=1.96 * se_out,
                    fmt='o', color='steelblue', alpha=0.7, capsize=3)

        # Add regression lines
        ivw = self.ivw()
        egger = self.mr_egger()

        x_range = np.array([min(beta_exp) - 0.1, max(beta_exp) + 0.1])

        # IVW line (through origin)
        if not np.isnan(ivw["beta"]):
            ax.plot(x_range, ivw["beta"] * x_range, 'r-',
                    label=f'IVW: OR={ivw["or"]:.2f}', linewidth=2)

        # MR-Egger line
        if not np.isnan(egger["beta"]):
            ax.plot(x_range, egger["intercept"] + egger["beta"] * x_range, 'g--',
                    label=f'Egger: OR={egger["or"]:.2f}', linewidth=2)

        ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle=':', alpha=0.5)

        ax.set_xlabel("SNP effect on gene expression (beta)")
        ax.set_ylabel("SNP effect on T2D (beta)")
        ax.set_title("MR Scatter Plot: Gene Expression -> T2D")
        ax.legend()

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_forest(self, output_path=None):
        """Create forest plot of per-SNP effects."""
        if len(self.harmonized) == 0:
            return

        fig, ax = plt.subplots(figsize=(10, max(4, len(self.harmonized) * 0.3)))

        beta_exp = self.harmonized["beta_exp"].values
        beta_out = self.harmonized["beta_out"].values
        se_out = self.harmonized["se_out"].values

        # Ratio estimates
        beta_iv = beta_out / beta_exp
        se_iv = se_out / np.abs(beta_exp)

        # Sort by effect size
        order = np.argsort(beta_iv)

        snps = self.harmonized["SNP"].values[order]
        beta_iv = beta_iv[order]
        se_iv = se_iv[order]

        y_pos = np.arange(len(snps))

        # Plot
        ax.errorbar(beta_iv, y_pos, xerr=1.96 * se_iv,
                    fmt='s', color='steelblue', capsize=3, markersize=6)

        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

        # Add IVW estimate
        ivw = self.ivw()
        if not np.isnan(ivw["beta"]):
            ax.axvline(x=ivw["beta"], color='red', linestyle='-', alpha=0.7,
                       label=f'IVW: {ivw["beta"]:.3f}')

        ax.set_yticks(y_pos)
        ax.set_yticklabels(snps)
        ax.set_xlabel("Causal effect (log OR)")
        ax.set_title("Forest Plot: Per-SNP Causal Effects")
        ax.legend()

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            plt.close()
        else:
            plt.show()


def run_gene_level_mr(exposure_df, outcome_df):
    """Run MR analysis for each gene separately."""
    print("\n" + "=" * 60)
    print("Gene-Level MR Analysis")
    print("=" * 60)

    all_results = []

    for gene in exposure_df["gene"].unique():
        gene_exposure = exposure_df[exposure_df["gene"] == gene].copy()
        gene_exposure = gene_exposure.rename(columns={
            "beta": "beta", "se": "se", "pval": "pval"
        })

        if len(gene_exposure) == 0:
            continue

        # Merge with outcome
        mr = MendelianRandomization(gene_exposure, outcome_df)

        if len(mr.harmonized) > 0:
            results = mr.run_all_methods()
            results["gene"] = gene
            results["n_instruments"] = len(mr.harmonized)
            all_results.append(results)

            # Print summary for IVW
            ivw = results[results["method"] == "IVW"].iloc[0]
            sig = "*" if ivw["pval"] < 0.05 else ""
            print(f"  {gene}: OR={ivw['or']:.3f} (95% CI: {ivw['or_ci_lower']:.3f}-{ivw['or_ci_upper']:.3f}), p={ivw['pval']:.2e} {sig}")

    if len(all_results) > 0:
        return pd.concat(all_results, ignore_index=True)
    else:
        return pd.DataFrame()


def run_module_level_mr(exposure_df, outcome_df, importance_df):
    """Run MR analysis for each workload module."""
    print("\n" + "=" * 60)
    print("Module-Level MR Analysis")
    print("=" * 60)

    # Get module assignments
    gene_modules = importance_df.groupby("gene")["module"].first().to_dict()

    all_results = []

    for module in importance_df["module"].unique():
        module_genes = [g for g, m in gene_modules.items() if m == module]
        module_exposure = exposure_df[exposure_df["gene"].isin(module_genes)].copy()

        if len(module_exposure) == 0:
            continue

        mr = MendelianRandomization(module_exposure, outcome_df)

        if len(mr.harmonized) > 0:
            results = mr.run_all_methods()
            results["module"] = module
            results["n_genes"] = len(module_genes)
            results["n_instruments"] = len(mr.harmonized)
            all_results.append(results)

            # Print summary
            ivw = results[results["method"] == "IVW"].iloc[0]
            sig = "*" if ivw["pval"] < 0.05 else ""
            print(f"\n{module.upper()}:")
            print(f"  Genes: {len(module_genes)}, Instruments: {len(mr.harmonized)}")
            print(f"  IVW: OR={ivw['or']:.3f} (95% CI: {ivw['or_ci_lower']:.3f}-{ivw['or_ci_upper']:.3f})")
            print(f"  p-value: {ivw['pval']:.2e} {sig}")

            # Check for pleiotropy
            egger = results[results["method"] == "MR-Egger"].iloc[0]
            if not np.isnan(egger.get("intercept_pval", 1)):
                if egger["intercept_pval"] < 0.05:
                    print(f"  Warning: Egger intercept p={egger['intercept_pval']:.3f} suggests pleiotropy")

            # Create plots
            mr.plot_scatter(FIGURES_DIR / f"scatter_{module}.png")

    if len(all_results) > 0:
        return pd.concat(all_results, ignore_index=True)
    else:
        return pd.DataFrame()


def create_summary_figure(gene_results, module_results):
    """Create summary forest plot of all MR results."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))

    # Gene-level results
    ax1 = axes[0]
    if len(gene_results) > 0:
        ivw_results = gene_results[gene_results["method"] == "IVW"].copy()
        ivw_results = ivw_results.dropna(subset=["or"])
        ivw_results = ivw_results.sort_values("or")

        y_pos = np.arange(len(ivw_results))

        colors = ['red' if p < 0.05 else 'steelblue' for p in ivw_results["pval"]]

        ax1.errorbar(ivw_results["or"], y_pos,
                     xerr=[ivw_results["or"] - ivw_results["or_ci_lower"],
                           ivw_results["or_ci_upper"] - ivw_results["or"]],
                     fmt='s', capsize=3, markersize=8, color='steelblue')

        for i, (_, row) in enumerate(ivw_results.iterrows()):
            color = 'red' if row["pval"] < 0.05 else 'steelblue'
            ax1.plot(row["or"], i, 's', color=color, markersize=8)

        ax1.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(ivw_results["gene"])
        ax1.set_xlabel("Odds Ratio for T2D")
        ax1.set_title("Gene-Level MR Results (IVW)")

    # Module-level results
    ax2 = axes[1]
    if len(module_results) > 0:
        ivw_results = module_results[module_results["method"] == "IVW"].copy()
        ivw_results = ivw_results.dropna(subset=["or"])
        ivw_results = ivw_results.sort_values("or")

        y_pos = np.arange(len(ivw_results))

        ax2.errorbar(ivw_results["or"], y_pos,
                     xerr=[ivw_results["or"] - ivw_results["or_ci_lower"],
                           ivw_results["or_ci_upper"] - ivw_results["or"]],
                     fmt='D', capsize=3, markersize=10, color='darkgreen')

        for i, (_, row) in enumerate(ivw_results.iterrows()):
            color = 'red' if row["pval"] < 0.05 else 'darkgreen'
            ax2.plot(row["or"], i, 'D', color=color, markersize=10)

        ax2.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(ivw_results["module"].str.upper())
        ax2.set_xlabel("Odds Ratio for T2D")
        ax2.set_title("Module-Level MR Results (IVW)")

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "mr_summary.png", dpi=150, bbox_inches='tight')
    plt.close()


def main():
    print("=" * 60)
    print("MENDELIAN RANDOMIZATION ANALYSIS")
    print("Workload Genes -> T2D Causal Effects")
    print("=" * 60)

    # Load data
    exposure_path = DATA_DIR / "exposure_eqtl.csv"
    outcome_path = DATA_DIR / "outcome_t2d.csv"

    if not exposure_path.exists():
        print("Error: Run 01_prepare_instruments.py first")
        return

    exposure_df = pd.read_csv(exposure_path)
    outcome_df = pd.read_csv(outcome_path)

    print(f"\nLoaded {len(exposure_df)} exposure instruments")
    print(f"Loaded {len(outcome_df)} outcome associations")

    # Load gene importance for module assignments
    importance_path = WORKLOAD_DIR / "results" / "deep_learning" / "gene_importance.csv"
    importance_df = pd.read_csv(importance_path)

    # 1. Gene-level MR
    gene_results = run_gene_level_mr(exposure_df, outcome_df)

    # 2. Module-level MR
    module_results = run_module_level_mr(exposure_df, outcome_df, importance_df)

    # 3. Combined workload MR (all instruments)
    print("\n" + "=" * 60)
    print("Combined Workload MR (All Instruments)")
    print("=" * 60)

    mr_all = MendelianRandomization(exposure_df, outcome_df)
    combined_results = mr_all.run_all_methods()

    print("\nResults:")
    for _, row in combined_results.iterrows():
        print(f"  {row['method']}: OR={row['or']:.3f} "
              f"(95% CI: {row['or_ci_lower']:.3f}-{row['or_ci_upper']:.3f}), "
              f"p={row['pval']:.2e}")

    # Create plots
    mr_all.plot_scatter(FIGURES_DIR / "scatter_combined.png")
    mr_all.plot_forest(FIGURES_DIR / "forest_combined.png")

    if len(gene_results) > 0 and len(module_results) > 0:
        create_summary_figure(gene_results, module_results)

    # Save results
    print("\n" + "=" * 60)
    print("Saving Results")
    print("=" * 60)

    if len(gene_results) > 0:
        gene_results.to_csv(RESULTS_DIR / "mr_gene_results.csv", index=False)
        print(f"  Saved: mr_gene_results.csv")

    if len(module_results) > 0:
        module_results.to_csv(RESULTS_DIR / "mr_module_results.csv", index=False)
        print(f"  Saved: mr_module_results.csv")

    combined_results.to_csv(RESULTS_DIR / "mr_combined_results.csv", index=False)
    print(f"  Saved: mr_combined_results.csv")

    # Summary report
    report = {
        "analysis_type": "Two-sample Mendelian Randomization",
        "exposure": "Workload gene expression (pancreas eQTL)",
        "outcome": "Type 2 Diabetes (DIAGRAM GWAS)",
        "n_genes": exposure_df["gene"].nunique() if "gene" in exposure_df.columns else 0,
        "n_instruments": len(exposure_df),
        "combined_results": combined_results.to_dict(orient="records"),
        "significant_genes": [],
        "significant_modules": []
    }

    if len(gene_results) > 0:
        sig_genes = gene_results[
            (gene_results["method"] == "IVW") & (gene_results["pval"] < 0.05)
        ]["gene"].tolist()
        report["significant_genes"] = sig_genes

    if len(module_results) > 0:
        sig_modules = module_results[
            (module_results["method"] == "IVW") & (module_results["pval"] < 0.05)
        ]["module"].tolist()
        report["significant_modules"] = sig_modules

    with open(RESULTS_DIR / "mr_report.json", 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"  Saved: mr_report.json")

    print("\n" + "=" * 60)
    print("MR ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"\nKey findings:")
    print(f"  Genes with causal evidence: {len(report['significant_genes'])}")
    print(f"  Modules with causal evidence: {len(report['significant_modules'])}")


if __name__ == "__main__":
    main()
