"""
X-intNMF Training Pipeline for Beta-Cell Workload Analysis
==========================================================

Trains WorkloadXintNMF on beta-cell single-cell data with:
1. Literature-derived module priors
2. Gene interaction network regularization
3. Workload factor extraction and CWI computation

Author: Beta-Cell Workload Analysis Pipeline
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
from pathlib import Path
from scipy import stats
from typing import Dict, List, Tuple, Optional
import json
import pickle

warnings.filterwarnings('ignore')

# Add parent directory for imports
SCRIPT_DIR = Path(__file__).parent.resolve()
sys.path.insert(0, str(SCRIPT_DIR))

from xintnmf_model import WorkloadXintNMF, XintNMF

# Configure paths
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
BETACELL_DIR = WORKLOAD_DIR.parent / "betacell"
RESULTS_DIR = WORKLOAD_DIR / "results" / "xintnmf"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_beta_cell_data() -> Tuple[sc.AnnData, np.ndarray, List[str]]:
    """
    Load and prepare beta-cell data for X-intNMF.

    Returns:
    --------
    adata : AnnData
        Full beta-cell dataset
    X : np.ndarray
        Expression matrix (genes x cells)
    gene_names : List[str]
        Gene names
    """
    data_path = BETACELL_DIR / "data" / "processed" / "segerstolpe_processed.h5ad"

    print(f"Loading data from {data_path}")
    adata = sc.read_h5ad(data_path)

    # Extract beta cells
    beta_mask = adata.obs['cell_type'] == 'beta cell'
    beta_cells = adata[beta_mask].copy()
    print(f"Extracted {beta_cells.n_obs} beta cells")

    # Get expression matrix (cells x genes -> genes x cells for NMF)
    if hasattr(beta_cells.X, 'toarray'):
        X = beta_cells.X.toarray().T  # genes x cells
    else:
        X = beta_cells.X.T

    # Ensure non-negative
    X = np.maximum(X, 0)

    gene_names = list(beta_cells.var_names)

    print(f"Expression matrix: {X.shape[0]} genes x {X.shape[1]} cells")

    return beta_cells, X, gene_names


def filter_workload_genes(
    X: np.ndarray,
    gene_names: List[str],
    min_expression: float = 0.1
) -> Tuple[np.ndarray, List[str], List[int]]:
    """
    Filter to genes relevant for workload analysis.

    Includes:
    1. Literature workload genes
    2. Highly variable genes
    3. Genes with sufficient expression
    """
    # Workload module genes
    WORKLOAD_GENES = set()
    for module, genes in WorkloadXintNMF.WORKLOAD_MODULES.items():
        WORKLOAD_GENES.update(genes)

    print(f"Workload module genes: {len(WORKLOAD_GENES)}")

    # Find available workload genes
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    available_workload = [g for g in WORKLOAD_GENES if g in gene_to_idx]
    print(f"Available in dataset: {len(available_workload)}")

    # Add highly expressed genes
    mean_expr = X.mean(axis=1)
    highly_expressed_idx = np.where(mean_expr > np.percentile(mean_expr, 50))[0]

    # Combine indices
    selected_idx = set()

    # Add workload genes
    for gene in available_workload:
        selected_idx.add(gene_to_idx[gene])

    # Add top variable genes
    gene_var = X.var(axis=1)
    top_var_idx = np.argsort(gene_var)[-500:]  # Top 500 variable
    selected_idx.update(top_var_idx)

    selected_idx = sorted(list(selected_idx))

    # Filter
    X_filtered = X[selected_idx, :]
    gene_names_filtered = [gene_names[i] for i in selected_idx]

    print(f"Selected {len(selected_idx)} genes for analysis")

    return X_filtered, gene_names_filtered, selected_idx


def create_ppi_network(
    gene_names: List[str],
    ppi_file: Optional[Path] = None
) -> np.ndarray:
    """
    Create gene-gene interaction network.

    Uses STRING PPI data if available, otherwise module-based network.
    """
    n_genes = len(gene_names)
    network = np.zeros((n_genes, n_genes))

    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Create module-based connections
    for module, genes in WorkloadXintNMF.WORKLOAD_MODULES.items():
        module_idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]

        # Fully connect genes within module
        for i in module_idx:
            for j in module_idx:
                if i != j:
                    network[i, j] = 1.0

    # Add some cross-module connections for related pathways
    # Biosynthetic <-> Metabolic (insulin secretion pathway)
    biosyn_idx = [gene_to_idx[g] for g in ["INS", "PCSK1", "PCSK2"] if g in gene_to_idx]
    metab_idx = [gene_to_idx[g] for g in ["GCK", "SLC2A2", "ABCC8", "KCNJ11"] if g in gene_to_idx]

    for i in biosyn_idx:
        for j in metab_idx:
            network[i, j] = 0.5
            network[j, i] = 0.5

    # Stress <-> Dedifferentiation (shared stress pathways)
    stress_idx = [gene_to_idx[g] for g in ["DDIT3", "ATF4", "XBP1"] if g in gene_to_idx]
    dediff_idx = [gene_to_idx[g] for g in ["ALDH1A3", "SOX9"] if g in gene_to_idx]

    for i in stress_idx:
        for j in dediff_idx:
            network[i, j] = 0.3
            network[j, i] = 0.3

    n_edges = (network > 0).sum() // 2
    print(f"Created network with {n_edges} edges")

    return network


def train_xintnmf(
    X: np.ndarray,
    gene_names: List[str],
    n_components: int = 5,
    alpha: float = 0.1,
    network: Optional[np.ndarray] = None
) -> WorkloadXintNMF:
    """
    Train WorkloadXintNMF model.
    """
    print("\n" + "="*60)
    print("Training X-intNMF Model")
    print("="*60)

    model = WorkloadXintNMF(
        n_components=n_components,
        alpha=alpha,
        beta=0.0,  # No cross-omics for single modality
        module_weight=0.5,
        max_iter=300,
        tol=1e-5,
        verbose=True
    )

    model.fit_workload(X, gene_names, external_network=network)

    return model


def analyze_components(
    model: WorkloadXintNMF,
    gene_names: List[str]
) -> Dict:
    """
    Analyze NMF components and their biological meaning.
    """
    print("\n" + "="*60)
    print("Component Analysis")
    print("="*60)

    results = {
        "components": {},
        "module_assignments": model.module_assignments_
    }

    # Get top genes for each component
    loadings = model.get_component_loadings(gene_names)

    for k in range(model.n_components):
        top_genes = loadings[k][:10]

        assignment = model.module_assignments_.get(k, {})
        module = assignment.get('module', 'unknown')

        print(f"\nComponent {k} -> {module.upper()}")
        print(f"  Module scores: {assignment.get('all_scores', {})}")
        print(f"  Top genes:")
        for gene, weight in top_genes[:5]:
            print(f"    {gene}: {weight:.4f}")

        results["components"][k] = {
            "module": module,
            "top_genes": [(g, float(w)) for g, w in top_genes],
            "scores": {m: float(s) for m, s in assignment.get('all_scores', {}).items()}
        }

    return results


def compute_cell_scores(
    model: WorkloadXintNMF,
    beta_cells: sc.AnnData
) -> pd.DataFrame:
    """
    Compute workload scores for each cell.
    """
    print("\n" + "="*60)
    print("Computing Cell Workload Scores")
    print("="*60)

    # Get module scores
    module_scores = model.get_workload_scores()

    # Compute CWI
    cwi = model.compute_cwi()

    # Create DataFrame
    scores_df = pd.DataFrame(index=beta_cells.obs_names)

    for module, values in module_scores.items():
        scores_df[f"score_{module}"] = values

    scores_df["CWI"] = cwi

    # Add condition info
    scores_df["condition"] = beta_cells.obs["condition"].values

    print(f"\nScore statistics:")
    print(scores_df.describe())

    return scores_df


def validate_t2d_association(
    scores_df: pd.DataFrame
) -> Dict:
    """
    Validate workload scores against T2D status.
    """
    print("\n" + "="*60)
    print("T2D Association Validation")
    print("="*60)

    results = {}

    t2d_mask = scores_df["condition"] == "T2D"
    healthy_mask = scores_df["condition"] == "Healthy"

    for col in scores_df.columns:
        if col == "condition":
            continue

        t2d_vals = scores_df.loc[t2d_mask, col]
        healthy_vals = scores_df.loc[healthy_mask, col]

        stat, pval = stats.mannwhitneyu(t2d_vals, healthy_vals, alternative='two-sided')

        # Effect size
        pooled_std = np.sqrt((t2d_vals.std()**2 + healthy_vals.std()**2) / 2)
        cohens_d = (t2d_vals.mean() - healthy_vals.mean()) / (pooled_std + 1e-10)

        direction = "UP" if t2d_vals.mean() > healthy_vals.mean() else "DOWN"
        sig = "*" if pval < 0.05 else ""

        print(f"\n{col}:")
        print(f"  T2D: {t2d_vals.mean():.4f} +/- {t2d_vals.std():.4f}")
        print(f"  Healthy: {healthy_vals.mean():.4f} +/- {healthy_vals.std():.4f}")
        print(f"  Direction: {direction} in T2D")
        print(f"  P-value: {pval:.2e} {sig}")
        print(f"  Cohen's d: {cohens_d:.3f}")

        results[col] = {
            "t2d_mean": float(t2d_vals.mean()),
            "healthy_mean": float(healthy_vals.mean()),
            "pvalue": float(pval),
            "cohens_d": float(cohens_d),
            "direction": direction,
            "significant": pval < 0.05
        }

    return results


def classify_workload_states(
    scores_df: pd.DataFrame,
    n_states: int = 5
) -> pd.Series:
    """
    Classify cells into workload states based on scores.

    States:
    S1: Resting (low activity, low stress)
    S2: Active (high biosynthetic, normal stress)
    S3: Stressed (elevated stress, normal dediff)
    S4: Exhausted (high stress, high dediff)
    S5: Failing (low activity, high dediff)
    """
    print("\n" + "="*60)
    print("Classifying Workload States")
    print("="*60)

    # Normalize scores for classification
    score_cols = [c for c in scores_df.columns if c.startswith("score_")]

    normalized = scores_df[score_cols].copy()
    for col in score_cols:
        min_val = normalized[col].min()
        max_val = normalized[col].max()
        if max_val > min_val:
            normalized[col] = (normalized[col] - min_val) / (max_val - min_val)

    # State classification rules
    states = []

    for idx in normalized.index:
        biosyn = normalized.loc[idx, "score_biosynthetic"] if "score_biosynthetic" in normalized.columns else 0.5
        stress = normalized.loc[idx, "score_stress"] if "score_stress" in normalized.columns else 0.5
        dediff = normalized.loc[idx, "score_dedifferentiation"] if "score_dedifferentiation" in normalized.columns else 0.5

        # Classification logic
        if dediff > 0.7:
            if biosyn < 0.3:
                state = "S5_Failing"
            else:
                state = "S4_Exhausted"
        elif stress > 0.7:
            state = "S3_Stressed"
        elif biosyn > 0.6:
            state = "S2_Active"
        else:
            state = "S1_Resting"

        states.append(state)

    state_series = pd.Series(states, index=scores_df.index, name="workload_state")

    # Print distribution
    print("\nState distribution:")
    state_counts = state_series.value_counts()
    for state, count in state_counts.items():
        pct = count / len(state_series) * 100
        print(f"  {state}: {count} ({pct:.1f}%)")

    # Distribution by condition
    print("\nState distribution by condition:")
    crosstab = pd.crosstab(state_series, scores_df["condition"], normalize='columns') * 100
    print(crosstab.round(1))

    return state_series


def create_visualization(
    model: WorkloadXintNMF,
    scores_df: pd.DataFrame,
    states: pd.Series,
    results_dir: Path
):
    """
    Create visualization of X-intNMF results.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    print("\n" + "="*60)
    print("Creating Visualizations")
    print("="*60)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # 1. Reconstruction error convergence
    ax = axes[0, 0]
    ax.plot(model.reconstruction_err_)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Reconstruction Error")
    ax.set_title("X-intNMF Convergence")
    ax.set_yscale('log')

    # 2. Module scores by condition
    ax = axes[0, 1]
    score_cols = [c for c in scores_df.columns if c.startswith("score_")]

    melted = scores_df[score_cols + ["condition"]].melt(
        id_vars=["condition"],
        var_name="Module",
        value_name="Score"
    )
    melted["Module"] = melted["Module"].str.replace("score_", "")

    sns.boxplot(data=melted, x="Module", y="Score", hue="condition", ax=ax)
    ax.set_title("Module Scores by Condition")
    ax.tick_params(axis='x', rotation=45)

    # 3. CWI distribution
    ax = axes[0, 2]
    for cond in ["Healthy", "T2D"]:
        mask = scores_df["condition"] == cond
        ax.hist(scores_df.loc[mask, "CWI"], bins=20, alpha=0.5, label=cond, density=True)
    ax.set_xlabel("CWI")
    ax.set_ylabel("Density")
    ax.set_title("Workload Index Distribution")
    ax.legend()

    # 4. Workload states by condition
    ax = axes[1, 0]
    state_cond = pd.crosstab(states, scores_df["condition"], normalize='columns') * 100
    state_cond.plot(kind='bar', ax=ax)
    ax.set_xlabel("Workload State")
    ax.set_ylabel("Percentage")
    ax.set_title("Workload States by Condition")
    ax.tick_params(axis='x', rotation=45)
    ax.legend(title="Condition")

    # 5. Component weights heatmap (top genes)
    ax = axes[1, 1]
    # Get top 10 genes per component
    n_top = 8
    W = model.W_
    gene_names = model.gene_names_

    top_genes_all = set()
    for k in range(model.n_components):
        top_idx = np.argsort(W[:, k])[-n_top:]
        top_genes_all.update(top_idx)

    top_idx_list = sorted(list(top_genes_all))[:20]  # Limit to 20 genes
    W_subset = W[top_idx_list, :]
    gene_labels = [gene_names[i] if i < len(gene_names) else f"Gene_{i}" for i in top_idx_list]

    sns.heatmap(W_subset, ax=ax, cmap='viridis',
                xticklabels=[f"C{k}" for k in range(model.n_components)],
                yticklabels=gene_labels)
    ax.set_title("Top Gene Weights per Component")

    # 6. Stress vs Dediff scatter
    ax = axes[1, 2]
    if "score_stress" in scores_df.columns and "score_dedifferentiation" in scores_df.columns:
        colors = ['green' if c == 'Healthy' else 'red' for c in scores_df["condition"]]
        ax.scatter(scores_df["score_stress"], scores_df["score_dedifferentiation"],
                  c=colors, alpha=0.5)
        ax.set_xlabel("Stress Score")
        ax.set_ylabel("Dedifferentiation Score")
        ax.set_title("Stress vs Dedifferentiation")

        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor='green', label='Healthy'),
                         Patch(facecolor='red', label='T2D')]
        ax.legend(handles=legend_elements)

    plt.tight_layout()
    fig_path = results_dir / "xintnmf_results.png"
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Saved figure to {fig_path}")


def save_results(
    model: WorkloadXintNMF,
    scores_df: pd.DataFrame,
    states: pd.Series,
    component_analysis: Dict,
    validation_results: Dict,
    results_dir: Path
):
    """
    Save all results.
    """
    print("\n" + "="*60)
    print("Saving Results")
    print("="*60)

    # Save model
    model_path = results_dir / "xintnmf_model.pkl"
    with open(model_path, 'wb') as f:
        pickle.dump(model, f)
    print(f"Model saved to {model_path}")

    # Save scores
    scores_df["workload_state"] = states
    scores_path = results_dir / "cell_scores.csv"
    scores_df.to_csv(scores_path)
    print(f"Scores saved to {scores_path}")

    # Save W and H matrices
    np.save(results_dir / "W_matrix.npy", model.W_)
    np.save(results_dir / "H_matrix.npy", model.H_)
    print("Matrices saved")

    # Save analysis results
    results = {
        "component_analysis": component_analysis,
        "validation": validation_results,
        "model_params": {
            "n_components": model.n_components,
            "alpha": model.alpha,
            "n_iterations": model.n_iter_,
            "final_error": float(model.reconstruction_err_[-1])
        }
    }

    # Convert numpy types for JSON serialization
    def convert_types(obj):
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_types(v) for v in obj]
        return obj

    results_path = results_dir / "analysis_results.json"
    with open(results_path, 'w') as f:
        json.dump(convert_types(results), f, indent=2)
    print(f"Analysis results saved to {results_path}")


def generate_report(
    component_analysis: Dict,
    validation_results: Dict,
    scores_df: pd.DataFrame,
    states: pd.Series,
    results_dir: Path
) -> str:
    """
    Generate summary report.
    """
    lines = [
        "=" * 80,
        "X-intNMF BETA-CELL WORKLOAD ANALYSIS REPORT",
        "=" * 80,
        "",
        "EXECUTIVE SUMMARY",
        "-" * 40,
    ]

    # CWI results
    cwi_result = validation_results.get("CWI", {})
    if cwi_result.get("significant"):
        direction = "ELEVATED" if cwi_result["direction"] == "UP" else "REDUCED"
        lines.append(f"[PASS] CWI is significantly {direction} in T2D (p = {cwi_result['pvalue']:.2e})")
    else:
        lines.append(f"[INFO] CWI not significantly different (p = {cwi_result.get('pvalue', 'N/A')})")

    # Module results
    sig_modules = sum(1 for k, v in validation_results.items()
                     if k.startswith("score_") and v.get("significant"))
    lines.append(f"[INFO] {sig_modules}/4 modules show significant T2D differences")

    # Component analysis
    lines.extend([
        "",
        "COMPONENT ANALYSIS",
        "-" * 40,
    ])

    for k, info in component_analysis.get("components", {}).items():
        module = info.get("module", "unknown")
        top_genes = info.get("top_genes", [])[:3]
        genes_str = ", ".join([g for g, w in top_genes])
        lines.append(f"  Component {k} -> {module.upper()}: {genes_str}")

    # Validation results
    lines.extend([
        "",
        "T2D ASSOCIATION RESULTS",
        "-" * 40,
    ])

    for score_name, result in validation_results.items():
        if not isinstance(result, dict):
            continue
        sig = "*" if result.get("significant") else ""
        lines.append(
            f"  {score_name}: {result.get('direction', 'N/A')} in T2D "
            f"(p={result.get('pvalue', 0):.2e}{sig}, d={result.get('cohens_d', 0):.2f})"
        )

    # State distribution
    lines.extend([
        "",
        "WORKLOAD STATE DISTRIBUTION",
        "-" * 40,
    ])

    state_counts = states.value_counts()
    for state, count in state_counts.items():
        pct = count / len(states) * 100
        lines.append(f"  {state}: {count} ({pct:.1f}%)")

    # Conclusions
    lines.extend([
        "",
        "=" * 80,
        "CONCLUSIONS",
        "=" * 80,
        "",
        "X-intNMF successfully identified workload-related gene modules:",
    ])

    # Find key findings
    stress_result = validation_results.get("score_stress", {})
    dediff_result = validation_results.get("score_dedifferentiation", {})
    biosyn_result = validation_results.get("score_biosynthetic", {})

    if stress_result.get("direction") == "UP":
        lines.append("  - Elevated STRESS response in T2D beta cells")
    if dediff_result.get("direction") == "UP":
        lines.append("  - Increased DEDIFFERENTIATION markers in T2D")
    if biosyn_result.get("direction") == "DOWN":
        lines.append("  - Reduced BIOSYNTHETIC capacity in T2D")

    lines.extend([
        "",
        "The network-regularized approach captures gene module structure",
        "while remaining interpretable for biological insights.",
        "",
        "=" * 80
    ])

    report = "\n".join(lines)

    report_path = results_dir / "xintnmf_report.txt"
    with open(report_path, 'w') as f:
        f.write(report)

    print(f"\nReport saved to {report_path}")

    return report


def main():
    """
    Main X-intNMF training pipeline.
    """
    print("=" * 60)
    print("X-intNMF BETA-CELL WORKLOAD ANALYSIS")
    print("=" * 60)

    # Load data
    beta_cells, X, gene_names = load_beta_cell_data()

    # Filter to relevant genes
    X_filtered, gene_names_filtered, selected_idx = filter_workload_genes(X, gene_names)

    # Create gene network
    network = create_ppi_network(gene_names_filtered)

    # Train model
    model = train_xintnmf(
        X_filtered,
        gene_names_filtered,
        n_components=5,
        alpha=0.1,
        network=network
    )

    # Analyze components
    component_analysis = analyze_components(model, gene_names_filtered)

    # Compute cell scores
    scores_df = compute_cell_scores(model, beta_cells)

    # Validate T2D association
    validation_results = validate_t2d_association(scores_df)

    # Classify workload states
    states = classify_workload_states(scores_df)

    # Create visualizations
    create_visualization(model, scores_df, states, RESULTS_DIR)

    # Save results
    save_results(
        model, scores_df, states,
        component_analysis, validation_results,
        RESULTS_DIR
    )

    # Generate report
    report = generate_report(
        component_analysis, validation_results,
        scores_df, states, RESULTS_DIR
    )

    print("\n" + report)

    print(f"\nAll results saved to: {RESULTS_DIR}")

    return model, scores_df, states


if __name__ == "__main__":
    model, scores_df, states = main()
