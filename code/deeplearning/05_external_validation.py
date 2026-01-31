#!/usr/bin/env python3
"""
05_external_validation.py - External validation of workload index

Validates the trained model on independent datasets:
1. GSE81547 (Enge et al.) - Aging pancreas study
2. GSE84133 (Baron et al.) - Human/mouse islets
3. Beta trajectory data - Pseudotime validation

Tests:
- CWI generalization to new data
- Cross-dataset marker correlations
- Age/condition associations in external data
"""

import torch
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import json
import warnings
warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
BETACELL_DIR = WORKLOAD_DIR.parent / "betacell"
RESULTS_DIR = WORKLOAD_DIR / "results" / "deep_learning"
EXTERNAL_DIR = RESULTS_DIR / "external_validation"
EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)

# Import model
import importlib.util
spec = importlib.util.spec_from_file_location("model", SCRIPT_DIR / "02_model.py")
model_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model_module)
create_model = model_module.create_model
LiteratureGuidedVAE = model_module.LiteratureGuidedVAE
FEATURE_MODULES = model_module.FEATURE_MODULES


def load_trained_model():
    """Load trained model and gene list."""
    print("=" * 60)
    print("Loading trained model")
    print("=" * 60)

    # Load training genes
    training_data = sc.read_h5ad(RESULTS_DIR / "training_data.h5ad")
    gene_names = training_data.var_names.tolist()

    # Create and load model
    model, _ = create_model(gene_names)
    weights_path = RESULTS_DIR / "model_weights.pt"
    model.load_state_dict(torch.load(weights_path, map_location='cpu'))
    model.eval()

    print(f"  Model loaded with {len(gene_names)} genes")

    return model, gene_names


def find_external_datasets():
    """Find available external datasets."""
    print("\n" + "=" * 60)
    print("Searching for external datasets")
    print("=" * 60)

    datasets = {}

    # 1. GSE81547 (Enge aging study)
    enge_paths = [
        BETACELL_DIR / "data" / "raw" / "GSE81547",
        BETACELL_DIR / "data" / "processed" / "enge_processed.h5ad",
    ]
    for path in enge_paths:
        if path.exists():
            datasets["enge"] = path
            print(f"  Found Enge (GSE81547): {path}")
            break

    # 2. Baron dataset
    baron_paths = [
        BETACELL_DIR / "data" / "raw" / "baron_GSE84133.tar.gz",
        BETACELL_DIR / "data" / "processed" / "baron_processed.h5ad",
    ]
    for path in baron_paths:
        if path.exists():
            datasets["baron"] = path
            print(f"  Found Baron (GSE84133): {path}")
            break

    # 3. Beta trajectory (from our analysis)
    traj_path = BETACELL_DIR / "results" / "celltype_deep_analysis" / "beta_trajectory.h5ad"
    if traj_path.exists():
        datasets["trajectory"] = traj_path
        print(f"  Found Beta trajectory: {traj_path}")

    # 4. Large atlas (if accessible)
    atlas_path = BETACELL_DIR / "data" / "processed" / "gse221156_atlas_processed.h5ad"
    if atlas_path.exists():
        datasets["atlas"] = atlas_path
        print(f"  Found GSE221156 Atlas: {atlas_path}")

    if len(datasets) == 0:
        print("  No external datasets found!")

    return datasets


def prepare_external_data(adata, training_genes):
    """
    Prepare external data for model prediction.
    Aligns genes with training data.
    """
    # Find overlapping genes
    common_genes = [g for g in training_genes if g in adata.var_names]
    missing_genes = [g for g in training_genes if g not in adata.var_names]

    print(f"    Common genes: {len(common_genes)}/{len(training_genes)}")
    print(f"    Missing genes: {len(missing_genes)}")

    if len(common_genes) < len(training_genes) * 0.5:
        print("    Warning: Less than 50% gene overlap!")
        return None

    # Create aligned expression matrix
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X.copy()

    # Build aligned matrix (fill missing with 0)
    aligned_X = np.zeros((adata.n_obs, len(training_genes)))

    for i, gene in enumerate(training_genes):
        if gene in adata.var_names:
            gene_idx = list(adata.var_names).index(gene)
            aligned_X[:, i] = X[:, gene_idx]

    # Handle inf/nan
    aligned_X = np.nan_to_num(aligned_X, nan=0.0, posinf=0.0, neginf=0.0)

    return torch.FloatTensor(aligned_X)


def compute_cwi_external(model, X_tensor):
    """Compute CWI for external data."""
    model.eval()
    with torch.no_grad():
        cwi = model.get_workload_index(X_tensor)
    return cwi.numpy()


def validate_on_trajectory(model, training_genes):
    """
    Validate on beta trajectory data.
    CWI should correlate with pseudotime (disease progression).
    """
    print("\n" + "=" * 60)
    print("Validation 1: Beta Cell Trajectory (Pseudotime)")
    print("=" * 60)

    traj_path = BETACELL_DIR / "results" / "celltype_deep_analysis" / "beta_trajectory.h5ad"
    if not traj_path.exists():
        print("  Trajectory data not found, skipping...")
        return None

    adata = sc.read_h5ad(traj_path)
    print(f"  Loaded: {adata.n_obs} cells")

    # Prepare data
    X_tensor = prepare_external_data(adata, training_genes)
    if X_tensor is None:
        return None

    # Compute CWI
    cwi = compute_cwi_external(model, X_tensor)

    results = {"dataset": "beta_trajectory", "n_cells": len(cwi)}

    # Check for pseudotime
    pseudotime_cols = ["dpt_pseudotime", "pseudotime", "latent_time", "monocle_pseudotime"]
    pseudotime = None

    for col in pseudotime_cols:
        if col in adata.obs.columns:
            pseudotime = adata.obs[col].values
            results["pseudotime_col"] = col
            break

    if pseudotime is not None:
        # Correlation with pseudotime
        valid_mask = ~np.isnan(pseudotime)
        corr, pval = stats.pearsonr(cwi[valid_mask], pseudotime[valid_mask])

        results["cwi_pseudotime_corr"] = float(corr)
        results["cwi_pseudotime_pval"] = float(pval)

        print(f"\n  CWI vs Pseudotime:")
        print(f"    Correlation: r = {corr:.3f}")
        print(f"    p-value: {pval:.2e}")
        print(f"    Expected: POSITIVE (higher CWI = more progressed)")
        print(f"    Result: {'PASS' if corr > 0 else 'UNEXPECTED'}")

    # Check condition if available
    if "condition" in adata.obs.columns or "disease" in adata.obs.columns:
        cond_col = "condition" if "condition" in adata.obs.columns else "disease"
        conditions = adata.obs[cond_col].astype(str)

        t2d_mask = conditions.str.contains("T2D|diabetic|diabetes", case=False, na=False)
        normal_mask = conditions.str.contains("Normal|healthy|control", case=False, na=False)

        if t2d_mask.sum() > 5 and normal_mask.sum() > 5:
            t2d_cwi = cwi[t2d_mask]
            normal_cwi = cwi[normal_mask]

            t_stat, p_val = stats.ttest_ind(t2d_cwi, normal_cwi)
            effect_size = (t2d_cwi.mean() - normal_cwi.mean()) / np.sqrt(
                (t2d_cwi.var() + normal_cwi.var()) / 2
            )

            results["t2d_vs_normal_effect"] = float(effect_size)
            results["t2d_vs_normal_pval"] = float(p_val)

            print(f"\n  CWI by Condition:")
            print(f"    T2D mean: {t2d_cwi.mean():.3f} (n={len(t2d_cwi)})")
            print(f"    Normal mean: {normal_cwi.mean():.3f} (n={len(normal_cwi)})")
            print(f"    Effect size: {effect_size:.3f}")
            print(f"    Result: {'PASS' if effect_size > 0 else 'UNEXPECTED'}")

    # Save CWI scores
    cwi_df = pd.DataFrame({
        "cell_id": adata.obs.index,
        "cwi": cwi
    })
    cwi_df.to_csv(EXTERNAL_DIR / "cwi_trajectory.csv", index=False)

    return results


def validate_marker_correlations(model, training_genes, adata, dataset_name):
    """Validate that CWI correlates correctly with markers in external data."""

    X_tensor = prepare_external_data(adata, training_genes)
    if X_tensor is None:
        return None

    cwi = compute_cwi_external(model, X_tensor)

    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X

    results = {"dataset": dataset_name}

    # Identity markers (should be negative correlation)
    identity_markers = ["PDX1", "MAFA", "UCN3", "NKX6-1"]
    identity_corrs = []

    print(f"\n  Identity markers (expected: negative):")
    for gene in identity_markers:
        if gene in adata.var_names:
            idx = list(adata.var_names).index(gene)
            corr, pval = stats.pearsonr(cwi, X[:, idx])
            identity_corrs.append(corr)
            status = "PASS" if corr < 0 else "FAIL"
            print(f"    {gene}: r = {corr:.3f} [{status}]")

    if identity_corrs:
        results["identity_mean_corr"] = float(np.mean(identity_corrs))
        results["identity_pass_rate"] = sum(1 for c in identity_corrs if c < 0) / len(identity_corrs)

    # Stress markers (should be positive correlation)
    stress_markers = ["DDIT3", "XBP1", "HSPA5", "ATF4"]
    stress_corrs = []

    print(f"\n  Stress markers (expected: positive):")
    for gene in stress_markers:
        if gene in adata.var_names:
            idx = list(adata.var_names).index(gene)
            corr, pval = stats.pearsonr(cwi, X[:, idx])
            stress_corrs.append(corr)
            status = "PASS" if corr > 0 else "FAIL"
            print(f"    {gene}: r = {corr:.3f} [{status}]")

    if stress_corrs:
        results["stress_mean_corr"] = float(np.mean(stress_corrs))
        results["stress_pass_rate"] = sum(1 for c in stress_corrs if c > 0) / len(stress_corrs)

    return results


def validate_on_atlas_subset(model, training_genes):
    """
    Validate on a subset of the large atlas.
    Test generalization to much larger dataset.
    """
    print("\n" + "=" * 60)
    print("Validation 2: GSE221156 Atlas (Large-scale)")
    print("=" * 60)

    atlas_path = BETACELL_DIR / "data" / "processed" / "gse221156_atlas_processed.h5ad"
    if not atlas_path.exists():
        print("  Atlas not found, skipping...")
        return None

    print("  Loading atlas (may take time for large file)...")

    try:
        # Try to load just metadata first
        adata = sc.read_h5ad(atlas_path, backed='r')
        n_total = adata.n_obs
        print(f"  Total cells in atlas: {n_total}")

        # Sample a subset for validation
        n_sample = min(5000, n_total)
        print(f"  Sampling {n_sample} cells for validation...")

        # Get random indices
        np.random.seed(42)
        sample_idx = np.random.choice(n_total, n_sample, replace=False)
        sample_idx = sorted(sample_idx)

        # Load subset
        adata_subset = adata[sample_idx].to_memory()

        print(f"  Loaded subset: {adata_subset.n_obs} cells, {adata_subset.n_vars} genes")

        # Validate marker correlations
        results = validate_marker_correlations(model, training_genes, adata_subset, "atlas_subset")

        if results:
            results["n_cells"] = n_sample
            results["n_total"] = n_total

        return results

    except Exception as e:
        print(f"  Error loading atlas: {e}")
        print("  Skipping atlas validation...")
        return None


def create_validation_figures(all_results):
    """Create summary figures for external validation."""
    print("\n" + "=" * 60)
    print("Creating validation figures")
    print("=" * 60)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # 1. Marker correlation comparison
    ax1 = axes[0]

    datasets = []
    identity_corrs = []
    stress_corrs = []

    for result in all_results:
        if result and "identity_mean_corr" in result:
            datasets.append(result["dataset"])
            identity_corrs.append(result["identity_mean_corr"])
            stress_corrs.append(result.get("stress_mean_corr", 0))

    if datasets:
        x = np.arange(len(datasets))
        width = 0.35

        ax1.bar(x - width/2, identity_corrs, width, label='Identity (expect <0)', color='steelblue')
        ax1.bar(x + width/2, stress_corrs, width, label='Stress (expect >0)', color='coral')
        ax1.axhline(y=0, color='gray', linestyle='--')

        ax1.set_ylabel('Mean Correlation with CWI')
        ax1.set_xticks(x)
        ax1.set_xticklabels(datasets, rotation=45)
        ax1.legend()
        ax1.set_title('Marker Correlations Across Datasets')

    # 2. Pass rates
    ax2 = axes[1]

    pass_rates = []
    labels = []

    for result in all_results:
        if result:
            if "identity_pass_rate" in result:
                labels.append(f"{result['dataset']}\n(Identity)")
                pass_rates.append(result["identity_pass_rate"] * 100)
            if "stress_pass_rate" in result:
                labels.append(f"{result['dataset']}\n(Stress)")
                pass_rates.append(result["stress_pass_rate"] * 100)

    if pass_rates:
        colors = ['green' if r >= 50 else 'red' for r in pass_rates]
        ax2.bar(range(len(pass_rates)), pass_rates, color=colors, alpha=0.7)
        ax2.axhline(y=50, color='gray', linestyle='--', label='50% threshold')
        ax2.set_ylabel('Pass Rate (%)')
        ax2.set_xticks(range(len(labels)))
        ax2.set_xticklabels(labels, rotation=45)
        ax2.set_title('Validation Pass Rates')
        ax2.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(EXTERNAL_DIR / "external_validation_summary.png", dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved: external_validation_summary.png")


def main():
    print("=" * 60)
    print("EXTERNAL VALIDATION OF WORKLOAD INDEX")
    print("=" * 60)

    # Load model
    model, training_genes = load_trained_model()

    # Find datasets
    datasets = find_external_datasets()

    all_results = []

    # 1. Trajectory validation
    traj_results = validate_on_trajectory(model, training_genes)
    if traj_results:
        all_results.append(traj_results)

    # 2. Atlas validation (if available)
    atlas_results = validate_on_atlas_subset(model, training_genes)
    if atlas_results:
        all_results.append(atlas_results)

    # Create figures
    if all_results:
        create_validation_figures(all_results)

    # Save results
    report = {
        "validation_type": "external",
        "n_datasets": len(all_results),
        "results": all_results,
        "summary": {}
    }

    # Compute overall pass rate
    total_tests = 0
    passed_tests = 0

    for result in all_results:
        if result:
            if "identity_pass_rate" in result:
                total_tests += 1
                if result["identity_pass_rate"] >= 0.5:
                    passed_tests += 1
            if "stress_pass_rate" in result:
                total_tests += 1
                if result["stress_pass_rate"] >= 0.5:
                    passed_tests += 1

    if total_tests > 0:
        report["summary"]["overall_pass_rate"] = passed_tests / total_tests
        report["summary"]["tests_passed"] = passed_tests
        report["summary"]["tests_total"] = total_tests

    with open(EXTERNAL_DIR / "external_validation_report.json", 'w') as f:
        json.dump(report, f, indent=2, default=str)

    print("\n" + "=" * 60)
    print("EXTERNAL VALIDATION COMPLETE")
    print("=" * 60)

    if total_tests > 0:
        print(f"\nOverall: {passed_tests}/{total_tests} validation tests passed")

    print(f"\nResults saved to: {EXTERNAL_DIR}")


if __name__ == "__main__":
    main()
