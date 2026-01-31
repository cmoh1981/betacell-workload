#!/usr/bin/env python3
"""
helpers.py - Utility functions for beta-cell analysis
"""

import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path


def filter_to_beta_cells(adata, cell_type_col: str = "cell_type"):
    """Filter AnnData to beta cells only."""
    beta_patterns = ["beta", "Beta", "BETA", "b cell", "B cell"]
    mask = adata.obs[cell_type_col].str.contains("|".join(beta_patterns), case=False, na=False)
    return adata[mask].copy()


def calculate_gene_scores(adata, gene_list: list, score_name: str):
    """Calculate module score for a gene list."""
    genes_found = [g for g in gene_list if g in adata.var_names]
    if genes_found:
        sc.tl.score_genes(adata, genes_found, score_name=score_name)
    return adata


def get_top_de_genes(de_df: pd.DataFrame, n: int = 50, direction: str = "both"):
    """Get top differentially expressed genes."""
    if direction == "up":
        return de_df[de_df["logfoldchanges"] > 0].nsmallest(n, "pvals_adj")
    elif direction == "down":
        return de_df[de_df["logfoldchanges"] < 0].nsmallest(n, "pvals_adj")
    else:
        return de_df.nsmallest(n, "pvals_adj")


def save_figure(fig, name: str, results_dir: Path, formats=["png", "pdf"]):
    """Save figure in multiple formats."""
    fig_dir = results_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(fig_dir / f"{name}.{fmt}", dpi=300, bbox_inches="tight")
