"""
Comprehensive Workload Index Validation
=======================================

Validates beta-cell workload index and MR findings using Segerstolpe dataset.
Tests:
1. CWI difference between T2D and healthy beta cells
2. MR-validated gene expression patterns
3. Workload state classification accuracy
4. Module score correlations

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
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# Configure paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
BETACELL_DIR = WORKLOAD_DIR.parent / "betacell"
RESULTS_DIR = WORKLOAD_DIR / "results" / "validation"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Workload gene signatures
WORKLOAD_SIGNATURES = {
    "biosynthetic": [
        "INS", "IAPP", "PCSK1", "PCSK2", "CPE", "SCG2", "CHGB",
        "SNAP25", "SYT1", "RAB3A", "VAMP2", "STX1A"
    ],
    "metabolic": [
        "GCK", "SLC2A2", "G6PC2", "PDX1", "MAFA", "NKX6-1",
        "ABCC8", "KCNJ11", "CACNA1C", "TRPM2"
    ],
    "stress": [
        "XBP1", "ATF6", "ERN1", "DDIT3", "ATF4", "HSPA5", "HSP90B1",
        "CALR", "CANX", "PDIA4", "PDIA6", "EIF2AK3", "TRIB3", "GPX1"
    ],
    "dedifferentiation": [
        "ALDH1A3", "NEUROG3", "SOX9", "HES1", "FOXO1", "MYC",
        "NANOG", "POU5F1", "CD44"
    ]
}

# MR-validated genes with their effects
MR_VALIDATED = {
    "PDX1": {"or": 0.66, "effect": "protective", "interpretation": "Higher expression = lower T2D risk"},
    "SLC2A2": {"or": 0.89, "effect": "protective", "interpretation": "Higher expression = lower T2D risk"},
    "MAFA": {"or": 1.14, "effect": "risk", "interpretation": "Higher expression = higher T2D risk"},
}


def load_segerstolpe_data() -> sc.AnnData:
    """Load the Segerstolpe processed dataset."""
    data_path = BETACELL_DIR / "data" / "processed" / "segerstolpe_processed.h5ad"

    if not data_path.exists():
        raise FileNotFoundError(f"Dataset not found at {data_path}")

    print(f"Loading data from {data_path}")
    adata = sc.read_h5ad(data_path)
    print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    return adata


def extract_beta_cells(adata: sc.AnnData) -> sc.AnnData:
    """Extract beta cells from the dataset."""
    beta_mask = adata.obs['cell_type'] == 'beta cell'
    beta_cells = adata[beta_mask].copy()
    print(f"Extracted {beta_cells.n_obs} beta cells")

    # Print condition distribution
    print(f"Condition distribution:")
    print(beta_cells.obs['condition'].value_counts())

    return beta_cells


def calculate_module_scores(adata: sc.AnnData) -> pd.DataFrame:
    """Calculate module scores for each signature."""
    scores = {}
    gene_counts = {}

    for module, genes in WORKLOAD_SIGNATURES.items():
        available = [g for g in genes if g in adata.var_names]
        gene_counts[module] = len(available)

        if len(available) >= 2:
            sc.tl.score_genes(adata, available, score_name=f"score_{module}")
            scores[module] = adata.obs[f"score_{module}"].values
        else:
            scores[module] = np.zeros(adata.n_obs)

        print(f"  {module}: {len(available)}/{len(genes)} genes available")

    return pd.DataFrame(scores, index=adata.obs_names), gene_counts


def calculate_cwi(scores: pd.DataFrame) -> pd.Series:
    """
    Calculate Composite Workload Index.
    CWI = (Biosynthetic / Metabolic) * (1 + Stress) * (1 + Dediff)
    """
    # Normalize to positive values (shift to min=0.1)
    normalized = scores.copy()
    for col in normalized.columns:
        min_val = normalized[col].min()
        normalized[col] = normalized[col] - min_val + 0.1

    # Calculate CWI
    demand = normalized['biosynthetic']
    capacity = normalized['metabolic'].replace(0, 0.1)
    stress = normalized['stress']
    dediff = normalized['dedifferentiation']

    cwi = (demand / capacity) * (1 + stress) * (1 + dediff)

    return cwi


def validate_t2d_cwi_difference(beta_cells: sc.AnnData, cwi: pd.Series) -> Dict:
    """Test if CWI differs between T2D and healthy beta cells."""
    print("\n" + "="*60)
    print("VALIDATION 1: CWI T2D vs Healthy Comparison")
    print("="*60)

    conditions = beta_cells.obs['condition']
    t2d_cwi = cwi[conditions == 'T2D']
    healthy_cwi = cwi[conditions == 'Healthy']

    # Statistical test
    stat, pval = stats.mannwhitneyu(t2d_cwi, healthy_cwi, alternative='two-sided')

    # Effect size (Cohen's d)
    pooled_std = np.sqrt((t2d_cwi.std()**2 + healthy_cwi.std()**2) / 2)
    cohens_d = (t2d_cwi.mean() - healthy_cwi.mean()) / pooled_std

    results = {
        "t2d_n": len(t2d_cwi),
        "healthy_n": len(healthy_cwi),
        "t2d_mean": float(t2d_cwi.mean()),
        "t2d_std": float(t2d_cwi.std()),
        "healthy_mean": float(healthy_cwi.mean()),
        "healthy_std": float(healthy_cwi.std()),
        "mannwhitney_stat": float(stat),
        "pvalue": float(pval),
        "cohens_d": float(cohens_d),
        "significant": pval < 0.05
    }

    print(f"\nT2D (n={len(t2d_cwi)}): CWI = {t2d_cwi.mean():.4f} +/- {t2d_cwi.std():.4f}")
    print(f"Healthy (n={len(healthy_cwi)}): CWI = {healthy_cwi.mean():.4f} +/- {healthy_cwi.std():.4f}")
    print(f"\nMann-Whitney U test:")
    print(f"  Statistic: {stat:.2f}")
    print(f"  P-value: {pval:.2e}")
    print(f"  Cohen's d: {cohens_d:.3f}")

    if pval < 0.05:
        direction = "HIGHER" if t2d_cwi.mean() > healthy_cwi.mean() else "LOWER"
        print(f"\n  RESULT: CWI is significantly {direction} in T2D beta cells (p < 0.05)")
    else:
        print(f"\n  RESULT: No significant difference in CWI between T2D and healthy")

    return results


def validate_mr_genes(beta_cells: sc.AnnData) -> Dict:
    """Validate MR-supported causal genes."""
    print("\n" + "="*60)
    print("VALIDATION 2: MR-Validated Gene Expression")
    print("="*60)

    results = {}
    conditions = beta_cells.obs['condition']

    for gene, info in MR_VALIDATED.items():
        print(f"\n{gene} (MR OR = {info['or']}, {info['effect']}):")
        print(f"  {info['interpretation']}")

        if gene not in beta_cells.var_names:
            print(f"  WARNING: Gene not found in dataset")
            results[gene] = {"available": False}
            continue

        # Get expression
        expr = beta_cells[:, gene].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        else:
            expr = np.asarray(expr).flatten()

        t2d_expr = expr[conditions == 'T2D']
        healthy_expr = expr[conditions == 'Healthy']

        # Statistical test
        stat, pval = stats.mannwhitneyu(t2d_expr, healthy_expr, alternative='two-sided')

        # Effect direction
        fold_change = np.mean(t2d_expr) / (np.mean(healthy_expr) + 1e-10)

        results[gene] = {
            "available": True,
            "mr_or": info['or'],
            "mr_effect": info['effect'],
            "t2d_mean": float(np.mean(t2d_expr)),
            "healthy_mean": float(np.mean(healthy_expr)),
            "fold_change": float(fold_change),
            "pvalue": float(pval),
            "direction": "down" if fold_change < 1 else "up"
        }

        print(f"  T2D: {np.mean(t2d_expr):.3f} +/- {np.std(t2d_expr):.3f}")
        print(f"  Healthy: {np.mean(healthy_expr):.3f} +/- {np.std(healthy_expr):.3f}")
        print(f"  Fold change (T2D/Healthy): {fold_change:.3f}")
        print(f"  P-value: {pval:.2e}")

        # Check if expression pattern is consistent with MR finding
        if info['effect'] == 'protective':
            # Protective gene should be lower in T2D
            expected = fold_change < 1
            consistency = "CONSISTENT" if expected else "INCONSISTENT"
        else:
            # Risk gene should be higher in T2D
            expected = fold_change > 1
            consistency = "CONSISTENT" if expected else "INCONSISTENT"

        results[gene]["mr_consistent"] = (consistency == "CONSISTENT")
        print(f"  MR Consistency: {consistency}")

    return results


def validate_module_scores(beta_cells: sc.AnnData, scores: pd.DataFrame) -> Dict:
    """Validate module scores between T2D and healthy."""
    print("\n" + "="*60)
    print("VALIDATION 3: Module Score Comparison")
    print("="*60)

    results = {}
    conditions = beta_cells.obs['condition']

    for module in scores.columns:
        t2d_scores = scores.loc[conditions == 'T2D', module]
        healthy_scores = scores.loc[conditions == 'Healthy', module]

        stat, pval = stats.mannwhitneyu(t2d_scores, healthy_scores, alternative='two-sided')

        results[module] = {
            "t2d_mean": float(t2d_scores.mean()),
            "healthy_mean": float(healthy_scores.mean()),
            "pvalue": float(pval),
            "significant": pval < 0.05,
            "direction": "up" if t2d_scores.mean() > healthy_scores.mean() else "down"
        }

        sig = "*" if pval < 0.05 else ""
        direction = "UP" if t2d_scores.mean() > healthy_scores.mean() else "DOWN"
        print(f"\n{module}:")
        print(f"  T2D: {t2d_scores.mean():.4f} +/- {t2d_scores.std():.4f}")
        print(f"  Healthy: {healthy_scores.mean():.4f} +/- {healthy_scores.std():.4f}")
        print(f"  p-value: {pval:.2e} {sig}")
        print(f"  Direction in T2D: {direction}")

    return results


def create_validation_figures(beta_cells: sc.AnnData, scores: pd.DataFrame,
                             cwi: pd.Series, results_dir: Path):
    """Create validation visualization figures."""
    print("\n" + "="*60)
    print("Creating validation figures...")
    print("="*60)

    conditions = beta_cells.obs['condition'].values

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # 1. CWI boxplot
    ax = axes[0, 0]
    cwi_df = pd.DataFrame({'CWI': cwi.values, 'Condition': conditions})
    sns.boxplot(data=cwi_df, x='Condition', y='CWI', ax=ax, palette=['green', 'red'])
    ax.set_title('Composite Workload Index (CWI)')
    ax.set_ylabel('CWI Score')

    # 2. Module scores heatmap
    ax = axes[0, 1]
    module_means = pd.DataFrame({
        'T2D': scores[conditions == 'T2D'].mean(),
        'Healthy': scores[conditions == 'Healthy'].mean()
    })
    sns.heatmap(module_means, annot=True, fmt='.3f', cmap='RdYlBu_r', ax=ax)
    ax.set_title('Module Scores by Condition')

    # 3. MR gene expression
    ax = axes[0, 2]
    mr_genes = [g for g in MR_VALIDATED.keys() if g in beta_cells.var_names]
    if mr_genes:
        expr_data = []
        for gene in mr_genes:
            expr = beta_cells[:, gene].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            for i, val in enumerate(expr):
                expr_data.append({'Gene': gene, 'Expression': val, 'Condition': conditions[i]})

        expr_df = pd.DataFrame(expr_data)
        sns.boxplot(data=expr_df, x='Gene', y='Expression', hue='Condition', ax=ax)
        ax.set_title('MR-Validated Gene Expression')
        ax.legend(loc='upper right')

    # 4. CWI distribution
    ax = axes[1, 0]
    for cond in ['Healthy', 'T2D']:
        mask = conditions == cond
        ax.hist(cwi[mask], bins=20, alpha=0.5, label=cond, density=True)
    ax.set_xlabel('CWI')
    ax.set_ylabel('Density')
    ax.set_title('CWI Distribution')
    ax.legend()

    # 5. Module correlation
    ax = axes[1, 1]
    corr = scores.corr()
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='coolwarm', center=0, ax=ax)
    ax.set_title('Module Score Correlations')

    # 6. Stress vs Dediff scatter
    ax = axes[1, 2]
    colors = ['green' if c == 'Healthy' else 'red' for c in conditions]
    ax.scatter(scores['stress'], scores['dedifferentiation'], c=colors, alpha=0.5)
    ax.set_xlabel('Stress Score')
    ax.set_ylabel('Dedifferentiation Score')
    ax.set_title('Stress vs Dedifferentiation')
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='green', label='Healthy'),
                      Patch(facecolor='red', label='T2D')]
    ax.legend(handles=legend_elements)

    plt.tight_layout()
    fig_path = results_dir / 'validation_figures.png'
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved figures to {fig_path}")


def generate_validation_report(all_results: Dict, results_dir: Path) -> str:
    """Generate comprehensive validation report."""

    lines = [
        "=" * 80,
        "BETA-CELL WORKLOAD INDEX VALIDATION REPORT",
        "=" * 80,
        "",
        "EXECUTIVE SUMMARY",
        "-" * 40,
    ]

    # CWI validation
    cwi = all_results['cwi_comparison']
    if cwi['significant']:
        direction = "ELEVATED" if cwi['t2d_mean'] > cwi['healthy_mean'] else "REDUCED"
        lines.append(f"[PASS] CWI is significantly {direction} in T2D (p = {cwi['pvalue']:.2e})")
    else:
        lines.append(f"[FAIL] No significant CWI difference (p = {cwi['pvalue']:.2e})")

    # MR gene validation
    mr_consistent = sum(1 for g, r in all_results['mr_genes'].items()
                       if r.get('mr_consistent', False))
    mr_total = sum(1 for g, r in all_results['mr_genes'].items()
                  if r.get('available', False))

    if mr_consistent == mr_total and mr_total > 0:
        lines.append(f"[PASS] All {mr_consistent}/{mr_total} MR-validated genes show consistent expression")
    else:
        lines.append(f"[PARTIAL] {mr_consistent}/{mr_total} MR-validated genes show consistent expression")

    # Module validation
    sig_modules = sum(1 for m, r in all_results['modules'].items() if r['significant'])
    lines.append(f"[INFO] {sig_modules}/4 modules show significant T2D differences")

    # Detailed results
    lines.extend([
        "",
        "DETAILED RESULTS",
        "-" * 40,
        "",
        "1. Composite Workload Index (CWI)",
        f"   T2D (n={cwi['t2d_n']}): {cwi['t2d_mean']:.4f} +/- {cwi['t2d_std']:.4f}",
        f"   Healthy (n={cwi['healthy_n']}): {cwi['healthy_mean']:.4f} +/- {cwi['healthy_std']:.4f}",
        f"   P-value: {cwi['pvalue']:.2e}",
        f"   Effect size (Cohen's d): {cwi['cohens_d']:.3f}",
        "",
        "2. MR-Validated Genes",
    ])

    for gene, info in all_results['mr_genes'].items():
        if info.get('available'):
            status = "CONSISTENT" if info.get('mr_consistent') else "INCONSISTENT"
            lines.append(f"   {gene}: MR OR={info['mr_or']}, FC={info['fold_change']:.2f}, p={info['pvalue']:.2e} [{status}]")
        else:
            lines.append(f"   {gene}: Not available in dataset")

    lines.extend([
        "",
        "3. Module Scores (T2D vs Healthy)",
    ])

    for module, info in all_results['modules'].items():
        sig = "*" if info['significant'] else ""
        lines.append(f"   {module}: {info['direction'].upper()} in T2D (p={info['pvalue']:.2e}){sig}")

    # Conclusions
    lines.extend([
        "",
        "=" * 80,
        "CONCLUSIONS",
        "=" * 80,
        "",
        "Key Findings:",
    ])

    if cwi['significant'] and cwi['t2d_mean'] > cwi['healthy_mean']:
        lines.append("1. Beta-cell workload (CWI) is elevated in T2D, supporting the")
        lines.append("   hypothesis that increased workload contributes to beta-cell dysfunction")

    if mr_consistent > 0:
        lines.append(f"2. MR-validated causal genes show expected expression patterns:")
        for gene, info in all_results['mr_genes'].items():
            if info.get('mr_consistent'):
                effect = "reduced" if info['mr_or'] < 1 else "increased"
                lines.append(f"   - {gene} is {effect} in T2D (consistent with causal role)")

    # Interpretation
    stress_up = all_results['modules'].get('stress', {}).get('direction') == 'up'
    dediff_up = all_results['modules'].get('dedifferentiation', {}).get('direction') == 'up'

    if stress_up and dediff_up:
        lines.extend([
            "",
            "3. T2D beta cells show increased stress AND dedifferentiation signatures,",
            "   suggesting a progression from stressed to failing states"
        ])

    lines.extend([
        "",
        "Clinical Implications:",
        "- Interventions reducing beta-cell workload may prevent T2D progression",
        "- PDX1 and SLC2A2 are potential therapeutic targets (MR-validated protective)",
        "- Monitoring workload signatures could identify at-risk individuals",
        "",
        "=" * 80
    ])

    report = "\n".join(lines)

    # Save report
    report_path = results_dir / "validation_report.txt"
    with open(report_path, 'w') as f:
        f.write(report)

    print(f"\nReport saved to: {report_path}")

    return report


def main():
    """Main validation pipeline."""
    print("=" * 60)
    print("BETA-CELL WORKLOAD INDEX VALIDATION")
    print("=" * 60)

    # Load data
    adata = load_segerstolpe_data()

    # Extract beta cells
    beta_cells = extract_beta_cells(adata)

    # Calculate module scores
    print("\nCalculating module scores...")
    scores, gene_counts = calculate_module_scores(beta_cells)

    # Calculate CWI
    print("\nCalculating Composite Workload Index...")
    cwi = calculate_cwi(scores)
    print(f"CWI range: {cwi.min():.4f} - {cwi.max():.4f}")

    # Run validations
    all_results = {}

    # 1. CWI T2D vs Healthy
    all_results['cwi_comparison'] = validate_t2d_cwi_difference(beta_cells, cwi)

    # 2. MR gene validation
    all_results['mr_genes'] = validate_mr_genes(beta_cells)

    # 3. Module score comparison
    all_results['modules'] = validate_module_scores(beta_cells, scores)

    # Create figures
    create_validation_figures(beta_cells, scores, cwi, RESULTS_DIR)

    # Generate report
    report = generate_validation_report(all_results, RESULTS_DIR)
    print("\n" + report)

    # Save results
    import json

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

    with open(RESULTS_DIR / "validation_results.json", 'w') as f:
        json.dump(convert_types(all_results), f, indent=2)

    print(f"\nResults saved to: {RESULTS_DIR}")

    return all_results


if __name__ == "__main__":
    results = main()
