"""
Multi-Dataset External Validation for Beta-Cell Workload Index
==============================================================

Validates our workload findings across multiple public single-cell datasets:
1. GSE81608 (Xin et al.) - 1,492 islet cells from ND/T2D donors
2. GSE86469 (Lawlor et al.) - 617 islet cells from ND/T2D donors
3. GSE84133 (Baron et al.) - Human pancreatic islet cells
4. GSE85241 (Muraro et al.) - Human pancreatic cells
5. PANC-DB/HPAP - 222,077 cells from 67 donors

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
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# Configure paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
RESULTS_DIR = WORKLOAD_DIR / "results" / "external_validation"
DATA_DIR = WORKLOAD_DIR / "data" / "external"

# Create directories
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Workload gene signatures from our analysis
WORKLOAD_SIGNATURES = {
    "biosynthetic_capacity": [
        "INS", "IAPP", "PCSK1", "PCSK2", "CPE", "SCG2", "CHGB",
        "SNAP25", "SYT1", "RAB3A", "VAMP2", "STX1A"
    ],
    "metabolic_sensing": [
        "GCK", "SLC2A2", "G6PC2", "PDX1", "MAFA", "NKX6-1",
        "ABCC8", "KCNJ11", "CACNA1C", "TRPM2"
    ],
    "stress_response": [
        "XBP1", "ATF6", "ERN1", "DDIT3", "ATF4", "HSPA5", "HSP90B1",
        "CALR", "CANX", "PDIA4", "PDIA6", "EIF2AK3", "TRIB3", "GPX1"
    ],
    "dedifferentiation": [
        "ALDH1A3", "NEUROG3", "SOX9", "HES1", "FOXO1", "MYC",
        "NANOG", "POU5F1", "CD44", "GLUT1"
    ],
    "inflammation": [
        "IL1B", "TNF", "IL6", "CCL2", "CXCL8", "NFKB1", "STAT3"
    ]
}

# MR-validated causal genes (from our analysis)
MR_VALIDATED_GENES = {
    "protective": ["PDX1", "SLC2A2"],  # OR < 1
    "risk": ["MAFA"],  # OR > 1
}

# Dataset information
DATASETS = {
    "GSE81608": {
        "name": "Xin et al. 2016",
        "description": "1,492 islet cells from ND/T2D donors",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81608",
        "cell_type_col": "cell_type",
        "condition_col": "disease",
        "beta_label": "beta"
    },
    "GSE86469": {
        "name": "Lawlor et al. 2017",
        "description": "617 islet cells from ND/T2D donors",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86469",
        "cell_type_col": "cell_type",
        "condition_col": "status",
        "beta_label": "Beta"
    },
    "GSE84133": {
        "name": "Baron et al. 2016",
        "description": "Human pancreatic islet cells",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133",
        "cell_type_col": "assigned_cluster",
        "condition_col": None,
        "beta_label": "beta"
    },
    "GSE85241": {
        "name": "Muraro et al. 2016",
        "description": "Human pancreatic cells CEL-Seq2",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241",
        "cell_type_col": "cell_type1",
        "condition_col": None,
        "beta_label": "beta"
    },
    "E-MTAB-5061": {
        "name": "Segerstolpe et al. 2016",
        "description": "Human islet scRNA-seq T2D vs healthy",
        "url": "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/",
        "cell_type_col": "cell_type",
        "condition_col": "disease",
        "beta_label": "beta"
    }
}


def download_geo_dataset(geo_id: str, output_dir: Path) -> Optional[Path]:
    """
    Download and preprocess a GEO dataset.
    Returns path to h5ad file or None if download fails.
    """
    h5ad_path = output_dir / f"{geo_id}.h5ad"

    if h5ad_path.exists():
        print(f"  Dataset {geo_id} already exists at {h5ad_path}")
        return h5ad_path

    print(f"  Attempting to download {geo_id}...")

    try:
        # Try using scanpy's built-in datasets first for common ones
        if geo_id == "GSE84133":
            # Baron dataset is available in scanpy
            adata = sc.datasets.baron_human()
            adata.write(h5ad_path)
            return h5ad_path
    except Exception as e:
        print(f"    Built-in download failed: {e}")

    # For other datasets, provide instructions
    print(f"""
    Manual download required for {geo_id}:
    1. Visit: {DATASETS[geo_id]['url']}
    2. Download the processed count matrix
    3. Save as {h5ad_path}

    Or use GEOparse:
    ```python
    import GEOparse
    gse = GEOparse.get_GEO(geo='{geo_id}', destdir='{output_dir}')
    ```
    """)

    return None


def calculate_module_scores(adata: sc.AnnData,
                           signatures: Dict[str, List[str]]) -> pd.DataFrame:
    """
    Calculate module scores for workload gene signatures.
    """
    scores = {}

    for module_name, genes in signatures.items():
        # Find available genes
        available = [g for g in genes if g in adata.var_names]

        if len(available) < 2:
            print(f"    Warning: Only {len(available)} genes available for {module_name}")
            scores[module_name] = np.zeros(adata.n_obs)
            continue

        # Calculate score using scanpy
        sc.tl.score_genes(adata, available, score_name=f"score_{module_name}")
        scores[module_name] = adata.obs[f"score_{module_name}"].values

        print(f"    {module_name}: {len(available)}/{len(genes)} genes available")

    return pd.DataFrame(scores, index=adata.obs_names)


def calculate_workload_index(scores: pd.DataFrame) -> pd.Series:
    """
    Calculate Composite Workload Index (CWI) from module scores.
    CWI = (Demand / Capacity) * Stress * (1 + Dediff)

    Where:
    - Demand = biosynthetic_capacity (normalized)
    - Capacity = metabolic_sensing (normalized)
    - Stress = stress_response (normalized)
    - Dediff = dedifferentiation (normalized)
    """
    # Normalize scores to 0-1 range
    normalized = scores.copy()
    for col in normalized.columns:
        min_val = normalized[col].min()
        max_val = normalized[col].max()
        if max_val > min_val:
            normalized[col] = (normalized[col] - min_val) / (max_val - min_val)
        else:
            normalized[col] = 0.5

    # Calculate CWI components
    demand = normalized.get('biosynthetic_capacity', pd.Series(0.5, index=scores.index))
    capacity = normalized.get('metabolic_sensing', pd.Series(0.5, index=scores.index))
    stress = normalized.get('stress_response', pd.Series(0.5, index=scores.index))
    dediff = normalized.get('dedifferentiation', pd.Series(0.5, index=scores.index))

    # Avoid division by zero
    capacity = capacity.replace(0, 0.01)

    # CWI formula
    cwi = (demand / capacity) * (1 + stress) * (1 + dediff)

    return cwi


def validate_t2d_association(adata: sc.AnnData,
                             dataset_info: dict,
                             results_dir: Path) -> Dict:
    """
    Validate workload index association with T2D status.
    """
    results = {
        "dataset": dataset_info["name"],
        "n_cells": adata.n_obs,
        "n_beta_cells": 0,
        "module_scores": {},
        "cwi_stats": {},
        "t2d_comparison": None,
        "gene_validation": {}
    }

    # Extract beta cells
    cell_type_col = dataset_info.get("cell_type_col")
    beta_label = dataset_info.get("beta_label", "beta")

    if cell_type_col and cell_type_col in adata.obs.columns:
        # Handle various beta cell labels
        beta_mask = adata.obs[cell_type_col].str.lower().str.contains('beta', na=False)
        beta_cells = adata[beta_mask].copy()
        results["n_beta_cells"] = beta_cells.n_obs
        print(f"  Found {beta_cells.n_obs} beta cells")
    else:
        beta_cells = adata.copy()
        results["n_beta_cells"] = adata.n_obs
        print(f"  No cell type column found, using all {adata.n_obs} cells")

    if beta_cells.n_obs < 10:
        print(f"  Warning: Too few beta cells ({beta_cells.n_obs}), skipping")
        return results

    # Calculate module scores
    print("  Calculating module scores...")
    scores = calculate_module_scores(beta_cells, WORKLOAD_SIGNATURES)
    results["module_scores"] = {
        col: {"mean": scores[col].mean(), "std": scores[col].std()}
        for col in scores.columns
    }

    # Calculate CWI
    cwi = calculate_workload_index(scores)
    results["cwi_stats"] = {
        "mean": float(cwi.mean()),
        "std": float(cwi.std()),
        "median": float(cwi.median()),
        "min": float(cwi.min()),
        "max": float(cwi.max())
    }

    # Compare T2D vs ND if condition info available
    condition_col = dataset_info.get("condition_col")
    if condition_col and condition_col in beta_cells.obs.columns:
        conditions = beta_cells.obs[condition_col].str.lower()

        # Identify T2D and ND cells
        t2d_mask = conditions.str.contains('t2d|type.*2|diabetic', na=False, regex=True)
        nd_mask = conditions.str.contains('normal|healthy|nd|control', na=False, regex=True)

        t2d_cwi = cwi[t2d_mask]
        nd_cwi = cwi[nd_mask]

        if len(t2d_cwi) > 5 and len(nd_cwi) > 5:
            # Statistical comparison
            stat, pval = stats.mannwhitneyu(t2d_cwi, nd_cwi, alternative='two-sided')
            effect_size = (t2d_cwi.mean() - nd_cwi.mean()) / np.sqrt(
                (t2d_cwi.std()**2 + nd_cwi.std()**2) / 2
            )

            results["t2d_comparison"] = {
                "n_t2d": len(t2d_cwi),
                "n_nd": len(nd_cwi),
                "t2d_cwi_mean": float(t2d_cwi.mean()),
                "nd_cwi_mean": float(nd_cwi.mean()),
                "mannwhitney_stat": float(stat),
                "pvalue": float(pval),
                "cohens_d": float(effect_size)
            }

            print(f"  T2D vs ND comparison:")
            print(f"    T2D (n={len(t2d_cwi)}): CWI = {t2d_cwi.mean():.3f} +/- {t2d_cwi.std():.3f}")
            print(f"    ND (n={len(nd_cwi)}): CWI = {nd_cwi.mean():.3f} +/- {nd_cwi.std():.3f}")
            print(f"    p-value = {pval:.2e}, Cohen's d = {effect_size:.3f}")

    # Validate MR-supported genes
    print("  Validating MR-supported genes...")
    for gene_type, genes in MR_VALIDATED_GENES.items():
        for gene in genes:
            if gene in beta_cells.var_names:
                expr = beta_cells[:, gene].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.asarray(expr).flatten()

                results["gene_validation"][gene] = {
                    "type": gene_type,
                    "mean_expr": float(np.mean(expr)),
                    "pct_expressed": float(np.sum(expr > 0) / len(expr) * 100),
                    "available": True
                }

                # Check correlation with CWI
                valid_mask = ~np.isnan(expr) & ~np.isnan(cwi.values)
                if np.sum(valid_mask) > 10:
                    corr, corr_p = stats.spearmanr(expr[valid_mask], cwi.values[valid_mask])
                    results["gene_validation"][gene]["cwi_correlation"] = float(corr)
                    results["gene_validation"][gene]["cwi_corr_pval"] = float(corr_p)
            else:
                results["gene_validation"][gene] = {
                    "type": gene_type,
                    "available": False
                }

    return results


def create_validation_report(all_results: List[Dict], output_path: Path):
    """
    Create a comprehensive validation report.
    """
    report_lines = [
        "=" * 80,
        "MULTI-DATASET EXTERNAL VALIDATION REPORT",
        "Beta-Cell Workload Index Validation",
        "=" * 80,
        "",
        "SUMMARY",
        "-" * 40,
    ]

    # Summary statistics
    total_cells = sum(r.get("n_beta_cells", 0) for r in all_results)
    datasets_with_t2d = sum(1 for r in all_results if r.get("t2d_comparison"))

    report_lines.extend([
        f"Total datasets analyzed: {len(all_results)}",
        f"Total beta cells: {total_cells:,}",
        f"Datasets with T2D comparison: {datasets_with_t2d}",
        "",
    ])

    # Per-dataset results
    report_lines.extend([
        "DATASET-SPECIFIC RESULTS",
        "-" * 40,
    ])

    for result in all_results:
        report_lines.extend([
            "",
            f"Dataset: {result['dataset']}",
            f"  Total cells: {result['n_cells']:,}",
            f"  Beta cells: {result['n_beta_cells']:,}",
        ])

        if result.get("cwi_stats"):
            cwi = result["cwi_stats"]
            report_lines.append(
                f"  CWI: {cwi['mean']:.3f} +/- {cwi['std']:.3f} "
                f"(range: {cwi['min']:.3f} - {cwi['max']:.3f})"
            )

        if result.get("t2d_comparison"):
            comp = result["t2d_comparison"]
            sig = "***" if comp["pvalue"] < 0.001 else "**" if comp["pvalue"] < 0.01 else "*" if comp["pvalue"] < 0.05 else ""
            report_lines.extend([
                f"  T2D vs ND Comparison:",
                f"    T2D (n={comp['n_t2d']}): CWI = {comp['t2d_cwi_mean']:.3f}",
                f"    ND (n={comp['n_nd']}): CWI = {comp['nd_cwi_mean']:.3f}",
                f"    p-value: {comp['pvalue']:.2e} {sig}",
                f"    Effect size (Cohen's d): {comp['cohens_d']:.3f}",
            ])

    # Gene validation summary
    report_lines.extend([
        "",
        "MR-VALIDATED GENE ANALYSIS",
        "-" * 40,
    ])

    gene_summary = {}
    for result in all_results:
        for gene, info in result.get("gene_validation", {}).items():
            if gene not in gene_summary:
                gene_summary[gene] = {"type": info.get("type"), "datasets": []}
            if info.get("available"):
                gene_summary[gene]["datasets"].append({
                    "dataset": result["dataset"],
                    "expr": info.get("mean_expr", 0),
                    "pct": info.get("pct_expressed", 0),
                    "cwi_corr": info.get("cwi_correlation", np.nan)
                })

    for gene, info in gene_summary.items():
        report_lines.append(f"\n{gene} ({info['type']} gene):")
        for ds in info["datasets"]:
            corr_str = f", CWI r={ds['cwi_corr']:.3f}" if not np.isnan(ds.get('cwi_corr', np.nan)) else ""
            report_lines.append(
                f"  {ds['dataset']}: mean={ds['expr']:.3f}, "
                f"% expressed={ds['pct']:.1f}%{corr_str}"
            )

    # Conclusions
    report_lines.extend([
        "",
        "=" * 80,
        "CONCLUSIONS",
        "=" * 80,
    ])

    # Check if CWI is elevated in T2D across datasets
    t2d_elevated = sum(
        1 for r in all_results
        if r.get("t2d_comparison") and r["t2d_comparison"]["t2d_cwi_mean"] > r["t2d_comparison"]["nd_cwi_mean"]
    )

    if t2d_elevated > 0 and datasets_with_t2d > 0:
        report_lines.append(
            f"- CWI elevated in T2D in {t2d_elevated}/{datasets_with_t2d} datasets with condition info"
        )

    # Write report
    report_text = "\n".join(report_lines)
    with open(output_path, 'w') as f:
        f.write(report_text)

    print(f"\nReport saved to: {output_path}")
    return report_text


def try_load_cellxgene_census():
    """
    Attempt to load data from CellxGene Census (large-scale human pancreas atlas).
    """
    try:
        import cellxgene_census
        print("\nAttempting to access CellxGene Census...")

        with cellxgene_census.open_soma() as census:
            # Query for pancreatic beta cells
            adata = cellxgene_census.get_anndata(
                census,
                organism="Homo sapiens",
                obs_value_filter="tissue_general == 'pancreas' and cell_type == 'type B pancreatic cell'",
                var_value_filter="feature_name in " + str(list(WORKLOAD_SIGNATURES["biosynthetic_capacity"][:5])),
                obs_column_names=["cell_type", "tissue", "disease", "donor_id"]
            )

            print(f"  Retrieved {adata.n_obs} beta cells from CellxGene Census")
            return adata

    except ImportError:
        print("  CellxGene Census not installed. Install with: pip install cellxgene-census")
    except Exception as e:
        print(f"  CellxGene Census access failed: {e}")

    return None


def main():
    """
    Main validation pipeline.
    """
    print("=" * 60)
    print("MULTI-DATASET EXTERNAL VALIDATION")
    print("Beta-Cell Workload Index")
    print("=" * 60)

    all_results = []

    # Try to load Baron dataset (available in scanpy)
    print("\n" + "-" * 40)
    print("1. Loading Baron et al. (GSE84133) dataset...")
    print("-" * 40)

    try:
        baron_path = DATA_DIR / "GSE84133.h5ad"
        if baron_path.exists():
            adata_baron = sc.read_h5ad(baron_path)
        else:
            print("  Downloading Baron dataset from scanpy...")
            adata_baron = sc.datasets.baron_human()
            adata_baron.write(baron_path)

        print(f"  Loaded {adata_baron.n_obs} cells, {adata_baron.n_vars} genes")

        # Run validation
        baron_results = validate_t2d_association(
            adata_baron,
            DATASETS["GSE84133"],
            RESULTS_DIR
        )
        all_results.append(baron_results)

    except Exception as e:
        print(f"  Error loading Baron dataset: {e}")

    # Try CellxGene Census
    print("\n" + "-" * 40)
    print("2. Attempting CellxGene Census access...")
    print("-" * 40)

    census_adata = try_load_cellxgene_census()
    if census_adata is not None:
        census_results = validate_t2d_association(
            census_adata,
            {
                "name": "CellxGene Census",
                "cell_type_col": "cell_type",
                "condition_col": "disease",
                "beta_label": "type B pancreatic cell"
            },
            RESULTS_DIR
        )
        all_results.append(census_results)

    # Check for manually downloaded datasets
    print("\n" + "-" * 40)
    print("3. Checking for manually downloaded datasets...")
    print("-" * 40)

    for geo_id, info in DATASETS.items():
        if geo_id == "GSE84133":  # Already processed
            continue

        dataset_path = DATA_DIR / f"{geo_id}.h5ad"
        if dataset_path.exists():
            print(f"\n  Found {geo_id}:")
            try:
                adata = sc.read_h5ad(dataset_path)
                print(f"    Loaded {adata.n_obs} cells, {adata.n_vars} genes")

                results = validate_t2d_association(adata, info, RESULTS_DIR)
                all_results.append(results)

            except Exception as e:
                print(f"    Error processing: {e}")
        else:
            print(f"  {geo_id} not found at {dataset_path}")
            print(f"    Download from: {info['url']}")

    # Generate report
    if all_results:
        print("\n" + "=" * 60)
        print("GENERATING VALIDATION REPORT")
        print("=" * 60)

        report_path = RESULTS_DIR / "multi_dataset_validation_report.txt"
        report = create_validation_report(all_results, report_path)
        print("\n" + report)

        # Save detailed results as JSON
        import json

        # Convert numpy types to Python types for JSON serialization
        def convert_types(obj):
            if isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert_types(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_types(v) for v in obj]
            return obj

        results_json = convert_types(all_results)

        with open(RESULTS_DIR / "validation_results.json", 'w') as f:
            json.dump(results_json, f, indent=2)

        print(f"\nDetailed results saved to: {RESULTS_DIR / 'validation_results.json'}")

    # Print download instructions for missing datasets
    print("\n" + "=" * 60)
    print("DOWNLOAD INSTRUCTIONS FOR ADDITIONAL DATASETS")
    print("=" * 60)
    print("""
To validate on additional datasets, download them to:
{data_dir}

Available datasets:
1. GSE81608 (Xin et al.) - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81608
2. GSE86469 (Lawlor et al.) - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86469
3. E-MTAB-5061 (Segerstolpe) - https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/
4. PANC-DB - https://hpap.pmacs.upenn.edu/

For easy download of processed h5ad files:
- CellxGene: https://cellxgene.cziscience.com/collections/a238e9fa-2bdf-41df-8522-69046f99baff
- Human Cell Atlas: https://explore.data.humancellatlas.org/projects/ae71be1d-ddd8-4feb-9bed-24c3ddb6e1ad

Install cellxgene-census for direct API access:
  pip install cellxgene-census
""".format(data_dir=DATA_DIR))

    return all_results


if __name__ == "__main__":
    results = main()
