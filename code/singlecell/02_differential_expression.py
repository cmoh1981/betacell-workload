#!/usr/bin/env python3
"""
02_differential_expression.py - Differential expression analysis
Beta-Cell T2D Workload Project
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import yaml

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

RESULTS = Path(config["paths"]["results"])

def run_de_analysis(adata, groupby: str, reference: str, target: str) -> pd.DataFrame:
    """Run differential expression between two groups."""
    print(f"Running DE: {target} vs {reference}")

    # Subset to groups of interest
    mask = adata.obs[groupby].isin([reference, target])
    adata_sub = adata[mask].copy()

    # Run rank_genes_groups
    sc.tl.rank_genes_groups(
        adata_sub,
        groupby=groupby,
        groups=[target],
        reference=reference,
        method="wilcoxon"
    )

    # Extract results
    result = sc.get.rank_genes_groups_df(adata_sub, group=target)
    result["comparison"] = f"{target}_vs_{reference}"

    return result

def analyze_key_genes(de_results: pd.DataFrame, gene_sets: dict) -> pd.DataFrame:
    """Analyze expression of key gene sets."""
    rows = []
    for set_name, genes in gene_sets.items():
        for gene in genes:
            gene_data = de_results[de_results["names"] == gene]
            if not gene_data.empty:
                rows.append({
                    "gene_set": set_name,
                    "gene": gene,
                    "logfoldchange": gene_data["logfoldchanges"].values[0],
                    "pval": gene_data["pvals"].values[0],
                    "pval_adj": gene_data["pvals_adj"].values[0],
                })
    return pd.DataFrame(rows)

def main():
    print("=" * 60)
    print("Differential Expression Analysis")
    print("=" * 60)

    # Load data (example with segerstolpe)
    adata = sc.read_h5ad(config["data_sources"]["segerstolpe"]["path"])

    # Run DE if condition column exists
    if "condition" in adata.obs.columns:
        de_results = run_de_analysis(adata, "condition", "Normal", "T2D")

        # Save full results
        output_path = RESULTS / "tables" / "de_results_T2D_vs_Normal.csv"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        de_results.to_csv(output_path, index=False)
        print(f"Saved DE results to {output_path}")

        # Analyze key genes
        key_genes = analyze_key_genes(de_results, config["gene_sets"])
        key_path = RESULTS / "tables" / "key_genes_de.csv"
        key_genes.to_csv(key_path, index=False)
        print(f"Saved key gene analysis to {key_path}")

if __name__ == "__main__":
    main()
