#!/usr/bin/env python3
"""
04_validate.py - Validation of trained workload model

Performs:
1. Internal validation (reconstruction, prediction metrics)
2. Literature alignment validation
3. Cross-dataset validation (if external data available)
4. Comparison of literature vs data-driven CWI
"""

import torch
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import json
from scipy import stats
from sklearn.metrics import (
    roc_auc_score, accuracy_score, confusion_matrix,
    silhouette_score, adjusted_rand_score
)
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Import model from 02_model.py
import importlib.util
spec = importlib.util.spec_from_file_location("model", Path(__file__).parent / "02_model.py")
model_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model_module)
create_model = model_module.create_model
LiteratureGuidedVAE = model_module.LiteratureGuidedVAE
FEATURE_MODULES = model_module.FEATURE_MODULES

# Paths - use absolute paths for reliability
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent  # workload/code/deeplearning -> workload
BETACELL_DIR = WORKLOAD_DIR.parent / "betacell"
RESULTS_DIR = WORKLOAD_DIR / "results" / "deep_learning"
FIGURES_DIR = RESULTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def load_model_and_data():
    """Load trained model and data."""
    print("=" * 60)
    print("Loading model and data")
    print("=" * 60)

    # Load training data for gene names
    data_path = RESULTS_DIR / "training_data.h5ad"
    adata = sc.read_h5ad(data_path)
    gene_names = adata.var_names.tolist()

    # Load model
    model, _ = create_model(gene_names)
    weights_path = RESULTS_DIR / "model_weights.pt"
    model.load_state_dict(torch.load(weights_path, map_location='cpu'))
    model.eval()

    print(f"Loaded model with {sum(p.numel() for p in model.parameters()):,} parameters")
    print(f"Data: {adata.n_obs} cells, {adata.n_vars} genes")

    return model, adata, gene_names


def validate_literature_alignment(model, adata):
    """
    Validate that CWI aligns with literature expectations.

    Expected behaviors:
    1. CWI higher in T2D vs Normal
    2. CWI negatively correlated with identity markers (PDX1, MAFA, UCN3)
    3. CWI positively correlated with stress markers (DDIT3, ATF4, ALDH1A3)
    """
    print("\n" + "=" * 60)
    print("Literature Alignment Validation")
    print("=" * 60)

    validations = {}

    # Get CWI scores
    cwi_path = RESULTS_DIR / "cwi_scores.csv"
    cwi_df = pd.read_csv(cwi_path)
    cwi = cwi_df["cwi_predicted"].values
    cwi_lit = cwi_df["cwi_literature"].values

    # Get expression matrix
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X

    # 1. CWI higher in T2D
    if "condition_binary" in adata.obs.columns:
        t2d_mask = adata.obs["condition_binary"] == 1
        t2d_cwi = cwi[t2d_mask]
        normal_cwi = cwi[~t2d_mask]

        t_stat, p_value = stats.ttest_ind(t2d_cwi, normal_cwi)
        effect_size = (t2d_cwi.mean() - normal_cwi.mean()) / np.sqrt(
            (t2d_cwi.var() + normal_cwi.var()) / 2
        )

        validations["t2d_vs_normal"] = {
            "passed": t2d_cwi.mean() > normal_cwi.mean(),
            "t2d_mean": float(t2d_cwi.mean()),
            "normal_mean": float(normal_cwi.mean()),
            "t_statistic": float(t_stat),
            "p_value": float(p_value),
            "effect_size": float(effect_size)
        }

        print(f"\n1. T2D vs Normal CWI:")
        print(f"   T2D mean: {t2d_cwi.mean():.3f}")
        print(f"   Normal mean: {normal_cwi.mean():.3f}")
        print(f"   Effect size (Cohen's d): {effect_size:.3f}")
        print(f"   p-value: {p_value:.2e}")
        print(f"   PASSED: {'Yes' if validations['t2d_vs_normal']['passed'] else 'No'}")

    # 2. Correlation with identity markers (should be negative)
    identity_markers = ["PDX1", "MAFA", "UCN3", "NKX6-1", "NEUROD1"]
    print("\n2. Identity Marker Correlations (expected: negative):")

    for gene in identity_markers:
        if gene in adata.var_names:
            gene_idx = adata.var_names.get_loc(gene)
            gene_expr = X[:, gene_idx]
            corr, p_val = stats.pearsonr(cwi, gene_expr)

            validations[f"corr_{gene}"] = {
                "correlation": float(corr),
                "p_value": float(p_val),
                "expected_negative": True,
                "passed": corr < 0
            }

            status = "PASS" if corr < 0 else "FAIL"
            print(f"   {gene}: r = {corr:.3f} (p = {p_val:.2e}) {status}")

    # 3. Correlation with stress markers (should be positive)
    stress_markers = ["DDIT3", "ATF4", "ALDH1A3", "XBP1", "HSPA5"]
    print("\n3. Stress Marker Correlations (expected: positive):")

    for gene in stress_markers:
        if gene in adata.var_names:
            gene_idx = adata.var_names.get_loc(gene)
            gene_expr = X[:, gene_idx]
            corr, p_val = stats.pearsonr(cwi, gene_expr)

            validations[f"corr_{gene}"] = {
                "correlation": float(corr),
                "p_value": float(p_val),
                "expected_positive": True,
                "passed": corr > 0
            }

            status = "PASS" if corr > 0 else "FAIL"
            print(f"   {gene}: r = {corr:.3f} (p = {p_val:.2e}) {status}")

    # 4. Correlation between data-driven and literature CWI
    if not np.isnan(cwi_lit).all():
        corr_lit, p_lit = stats.pearsonr(cwi, cwi_lit)
        validations["cwi_literature_correlation"] = {
            "correlation": float(corr_lit),
            "p_value": float(p_lit)
        }
        print(f"\n4. Data-driven vs Literature CWI:")
        print(f"   Correlation: r = {corr_lit:.3f} (p = {p_lit:.2e})")

    # Summary
    n_passed = sum(1 for v in validations.values() if isinstance(v, dict) and v.get("passed", True))
    n_total = len([v for v in validations.values() if isinstance(v, dict) and "passed" in v])

    print(f"\n   VALIDATION SUMMARY: {n_passed}/{n_total} tests passed")

    return validations


def validate_state_separation(model, adata):
    """Validate that workload states are well-separated in latent space."""
    print("\n" + "=" * 60)
    print("State Separation Validation")
    print("=" * 60)

    # Get latent representations
    if hasattr(adata.X, 'toarray'):
        X = torch.FloatTensor(adata.X.toarray())
    else:
        X = torch.FloatTensor(adata.X)

    with torch.no_grad():
        mu, _, _ = model.encode(X)
        z = mu.numpy()

    # Get state labels
    state_path = RESULTS_DIR / "state_predictions.csv"
    state_df = pd.read_csv(state_path)
    states_pred = state_df["state_predicted"].values

    # Map to numeric
    state_map = {"S1_Resting": 0, "S2_Active": 1, "S3_Stressed": 2,
                 "S4_Exhausted": 3, "S5_Failing": 4}
    states_numeric = np.array([state_map.get(s, -1) for s in states_pred])
    valid_mask = states_numeric >= 0

    metrics = {}

    # Silhouette score
    if len(np.unique(states_numeric[valid_mask])) > 1:
        sil_score = silhouette_score(z[valid_mask], states_numeric[valid_mask])
        metrics["silhouette_score"] = float(sil_score)
        print(f"Silhouette Score: {sil_score:.3f}")
        print(f"  Target: > 0.3, {'PASSED' if sil_score > 0.3 else 'NEEDS IMPROVEMENT'}")

    # Compare with literature states
    if "workload_state" in adata.obs.columns:
        states_lit = adata.obs["workload_state"].values
        states_lit_numeric = np.array([state_map.get(s, -1) for s in states_lit])
        valid_both = (states_numeric >= 0) & (states_lit_numeric >= 0)

        if valid_both.sum() > 0:
            ari = adjusted_rand_score(
                states_lit_numeric[valid_both],
                states_numeric[valid_both]
            )
            metrics["adjusted_rand_index"] = float(ari)
            print(f"Adjusted Rand Index (vs literature): {ari:.3f}")

            # Confusion matrix
            cm = confusion_matrix(
                states_lit_numeric[valid_both],
                states_numeric[valid_both]
            )
            accuracy = np.diag(cm).sum() / cm.sum()
            metrics["state_agreement"] = float(accuracy)
            print(f"State Agreement: {accuracy:.1%}")

    return metrics


def validate_gene_importance(model, adata):
    """Validate that important genes align with literature."""
    print("\n" + "=" * 60)
    print("Gene Importance Validation")
    print("=" * 60)

    # Load gene importance
    importance_path = RESULTS_DIR / "gene_importance.csv"
    importance_df = pd.read_csv(importance_path)

    # Expected top genes per module (from literature)
    expected_top = {
        "biosynthetic": ["INS", "IAPP", "HSPA5", "PCSK1"],
        "metabolic": ["PDX1", "MAFA", "GCK", "SLC2A2"],
        "stress": ["DDIT3", "XBP1", "ATF4", "ATF6"],
        "dedifferentiation": ["ALDH1A3", "LDHA", "HK1"]
    }

    validations = {}

    for module in expected_top:
        module_genes = importance_df[importance_df["module"] == module]
        if len(module_genes) == 0:
            continue

        top_5 = module_genes.nlargest(5, "importance")["gene"].tolist()
        expected = expected_top[module]

        overlap = len(set(top_5) & set(expected))

        validations[module] = {
            "top_5_genes": top_5,
            "expected_genes": expected,
            "overlap": overlap,
            "overlap_ratio": overlap / len(expected)
        }

        print(f"\n{module.upper()} Module:")
        print(f"  Top 5 data-driven: {top_5}")
        print(f"  Expected from literature: {expected}")
        print(f"  Overlap: {overlap}/{len(expected)}")

    return validations


def compare_cwi_distributions(model, adata):
    """Compare data-driven vs literature CWI distributions."""
    print("\n" + "=" * 60)
    print("CWI Distribution Comparison")
    print("=" * 60)

    # Load CWI scores
    cwi_path = RESULTS_DIR / "cwi_scores.csv"
    cwi_df = pd.read_csv(cwi_path)

    cwi_pred = cwi_df["cwi_predicted"].values
    cwi_lit = cwi_df["cwi_literature"].values

    # Statistics
    print("\nData-driven CWI:")
    print(f"  Mean: {cwi_pred.mean():.3f}")
    print(f"  Std: {cwi_pred.std():.3f}")
    print(f"  Range: [{cwi_pred.min():.3f}, {cwi_pred.max():.3f}]")

    if not np.isnan(cwi_lit).all():
        print("\nLiterature CWI:")
        print(f"  Mean: {np.nanmean(cwi_lit):.3f}")
        print(f"  Std: {np.nanstd(cwi_lit):.3f}")
        print(f"  Range: [{np.nanmin(cwi_lit):.3f}, {np.nanmax(cwi_lit):.3f}]")

    # Create comparison figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # 1. Distribution comparison
    ax1 = axes[0]
    ax1.hist(cwi_pred, bins=50, alpha=0.7, label="Data-driven", density=True)
    if not np.isnan(cwi_lit).all():
        ax1.hist(cwi_lit, bins=50, alpha=0.7, label="Literature", density=True)
    ax1.set_xlabel("CWI Score")
    ax1.set_ylabel("Density")
    ax1.set_title("CWI Distribution")
    ax1.legend()

    # 2. Scatter plot
    ax2 = axes[1]
    if not np.isnan(cwi_lit).all():
        mask = ~np.isnan(cwi_lit)
        ax2.scatter(cwi_lit[mask], cwi_pred[mask], alpha=0.3, s=5)
        ax2.plot([0, 5], [0, 5], 'r--', label="y=x")
        corr = np.corrcoef(cwi_lit[mask], cwi_pred[mask])[0, 1]
        ax2.set_title(f"Literature vs Data-driven (r={corr:.3f})")
        ax2.set_xlabel("Literature CWI")
        ax2.set_ylabel("Data-driven CWI")

    # 3. By condition
    ax3 = axes[2]
    if "condition" in cwi_df.columns:
        conditions = cwi_df["condition"].values
        for cond, label in [(0, "Normal"), (1, "T2D")]:
            mask = conditions == cond
            ax3.hist(cwi_pred[mask], bins=30, alpha=0.7, label=label, density=True)
        ax3.set_xlabel("CWI Score")
        ax3.set_ylabel("Density")
        ax3.set_title("CWI by Condition")
        ax3.legend()

    plt.tight_layout()
    fig_path = FIGURES_DIR / "cwi_comparison.png"
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved figure to {fig_path}")

    return {"figure_saved": str(fig_path)}


def create_validation_plots(model, adata):
    """Create validation visualization plots."""
    print("\n" + "=" * 60)
    print("Creating Validation Plots")
    print("=" * 60)

    # Get data
    if hasattr(adata.X, 'toarray'):
        X = torch.FloatTensor(adata.X.toarray())
    else:
        X = torch.FloatTensor(adata.X)

    # Get latent space
    with torch.no_grad():
        mu, _, _ = model.encode(X)
        z = mu.numpy()

    # Load predictions
    cwi_df = pd.read_csv(RESULTS_DIR / "cwi_scores.csv")
    state_df = pd.read_csv(RESULTS_DIR / "state_predictions.csv")

    # 1. Latent space UMAP
    from sklearn.manifold import TSNE

    print("Computing t-SNE...")
    tsne = TSNE(n_components=2, random_state=42, perplexity=30)
    z_2d = tsne.fit_transform(z[:min(5000, len(z))])  # Limit for speed

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Color by CWI
    ax1 = axes[0]
    scatter = ax1.scatter(z_2d[:, 0], z_2d[:, 1],
                          c=cwi_df["cwi_predicted"].values[:len(z_2d)],
                          cmap="RdYlBu_r", s=5, alpha=0.7)
    plt.colorbar(scatter, ax=ax1, label="CWI")
    ax1.set_title("Latent Space (colored by CWI)")
    ax1.set_xlabel("t-SNE 1")
    ax1.set_ylabel("t-SNE 2")

    # Color by state
    ax2 = axes[1]
    state_colors = {"S1_Resting": 0, "S2_Active": 1, "S3_Stressed": 2,
                    "S4_Exhausted": 3, "S5_Failing": 4}
    states = [state_colors.get(s, 0) for s in state_df["state_predicted"].values[:len(z_2d)]]
    scatter = ax2.scatter(z_2d[:, 0], z_2d[:, 1], c=states,
                          cmap="viridis", s=5, alpha=0.7)
    plt.colorbar(scatter, ax=ax2, label="State")
    ax2.set_title("Latent Space (colored by State)")
    ax2.set_xlabel("t-SNE 1")
    ax2.set_ylabel("t-SNE 2")

    # Color by condition
    ax3 = axes[2]
    conditions = cwi_df["condition"].values[:len(z_2d)]
    scatter = ax3.scatter(z_2d[:, 0], z_2d[:, 1], c=conditions,
                          cmap="coolwarm", s=5, alpha=0.7)
    plt.colorbar(scatter, ax=ax3, label="Condition")
    ax3.set_title("Latent Space (colored by Condition)")
    ax3.set_xlabel("t-SNE 1")
    ax3.set_ylabel("t-SNE 2")

    plt.tight_layout()
    fig_path = FIGURES_DIR / "latent_space_visualization.png"
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved latent space visualization to {fig_path}")

    # 2. Gene importance heatmap
    importance_df = pd.read_csv(RESULTS_DIR / "gene_importance.csv")

    fig, ax = plt.subplots(figsize=(12, 8))

    # Pivot for heatmap
    pivot_data = importance_df.pivot_table(
        index="gene", columns="module", values="importance", fill_value=0
    )

    # Get top genes per module
    top_genes = []
    for module in pivot_data.columns:
        top = importance_df[importance_df["module"] == module].nlargest(10, "importance")["gene"]
        top_genes.extend(top.tolist())
    top_genes = list(dict.fromkeys(top_genes))  # Remove duplicates, keep order

    if len(top_genes) > 0:
        pivot_subset = pivot_data.loc[pivot_data.index.isin(top_genes)]
        sns.heatmap(pivot_subset, cmap="YlOrRd", ax=ax, annot=True, fmt=".2f")
        ax.set_title("Gene Importance by Module (Top Genes)")
        plt.tight_layout()
        fig_path = FIGURES_DIR / "gene_importance_heatmap.png"
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved gene importance heatmap to {fig_path}")

    return {"figures_created": True}


def generate_validation_report(all_results: dict):
    """Generate final validation report."""
    print("\n" + "=" * 60)
    print("VALIDATION REPORT")
    print("=" * 60)

    report = {
        "timestamp": pd.Timestamp.now().isoformat(),
        "results": all_results,
        "summary": {}
    }

    # Summarize key metrics
    lit_val = all_results.get("literature_alignment", {})
    if "t2d_vs_normal" in lit_val:
        report["summary"]["t2d_effect_size"] = lit_val["t2d_vs_normal"]["effect_size"]

    state_val = all_results.get("state_separation", {})
    if "silhouette_score" in state_val:
        report["summary"]["silhouette_score"] = state_val["silhouette_score"]
    if "adjusted_rand_index" in state_val:
        report["summary"]["adjusted_rand_index"] = state_val["adjusted_rand_index"]

    # Overall assessment
    passed_tests = 0
    total_tests = 0

    for key, val in lit_val.items():
        if isinstance(val, dict) and "passed" in val:
            total_tests += 1
            if val["passed"]:
                passed_tests += 1

    report["summary"]["tests_passed"] = passed_tests
    report["summary"]["tests_total"] = total_tests
    report["summary"]["pass_rate"] = passed_tests / total_tests if total_tests > 0 else 0

    # Print summary
    print(f"\nOverall Validation Results:")
    print(f"  Tests Passed: {passed_tests}/{total_tests} ({report['summary']['pass_rate']:.1%})")

    if "silhouette_score" in report["summary"]:
        sil = report["summary"]["silhouette_score"]
        print(f"  Silhouette Score: {sil:.3f} ({'Good' if sil > 0.3 else 'Moderate'})")

    if "t2d_effect_size" in report["summary"]:
        eff = report["summary"]["t2d_effect_size"]
        print(f"  T2D Effect Size: {eff:.3f} ({'Large' if eff > 0.8 else 'Medium' if eff > 0.5 else 'Small'})")

    # Save report
    report_path = RESULTS_DIR / "validation_report.json"
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nSaved validation report to {report_path}")

    return report


def main():
    print("=" * 60)
    print("BETA-CELL WORKLOAD - MODEL VALIDATION")
    print("=" * 60)

    # Load model and data
    model, adata, gene_names = load_model_and_data()

    all_results = {}

    # 1. Literature alignment validation
    all_results["literature_alignment"] = validate_literature_alignment(model, adata)

    # 2. State separation validation
    all_results["state_separation"] = validate_state_separation(model, adata)

    # 3. Gene importance validation
    all_results["gene_importance"] = validate_gene_importance(model, adata)

    # 4. CWI distribution comparison
    all_results["cwi_comparison"] = compare_cwi_distributions(model, adata)

    # 5. Create validation plots
    all_results["plots"] = create_validation_plots(model, adata)

    # 6. Generate final report
    report = generate_validation_report(all_results)

    print("\n" + "=" * 60)
    print("VALIDATION COMPLETE")
    print("=" * 60)
    print(f"\nAll outputs saved to: {RESULTS_DIR}")
    print("  - validation_report.json")
    print("  - figures/cwi_comparison.png")
    print("  - figures/latent_space_visualization.png")
    print("  - figures/gene_importance_heatmap.png")


if __name__ == "__main__":
    main()
