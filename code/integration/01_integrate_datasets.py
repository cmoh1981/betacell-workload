#!/usr/bin/env python3
"""
01_integrate_datasets.py - Multi-dataset integration
Beta-Cell T2D Workload Project
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
import yaml

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

DATA_PROCESSED = Path(config["paths"]["data_processed"])
RESULTS = Path(config["paths"]["results"])

def harmonize_genes(adatas: list, names: list) -> list:
    """Find common genes across datasets."""
    gene_sets = [set(adata.var_names) for adata in adatas]
    common_genes = gene_sets[0].intersection(*gene_sets[1:])
    print(f"Common genes across {len(adatas)} datasets: {len(common_genes)}")

    harmonized = []
    for adata, name in zip(adatas, names):
        adata_sub = adata[:, list(common_genes)].copy()
        adata_sub.obs["dataset"] = name
        harmonized.append(adata_sub)

    return harmonized

def integrate_with_harmony(adata: ad.AnnData, batch_key: str = "dataset") -> ad.AnnData:
    """Integrate datasets using Harmony."""
    import harmonypy

    print("Running Harmony integration...")

    # Ensure PCA is computed
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata, n_comps=config["analysis"]["n_pcs"])

    # Run Harmony
    ho = harmonypy.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        batch_key,
        max_iter_harmony=20
    )
    adata.obsm["X_pca_harmony"] = ho.Z_corr.T

    # Recompute neighbors and UMAP on corrected PCs
    sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=config["analysis"]["resolution"])

    return adata

def main():
    print("=" * 60)
    print("Multi-Dataset Integration")
    print("=" * 60)

    adatas = []
    names = []

    for name, source in config["data_sources"].items():
        try:
            path = Path(source["path"])
            adata = sc.read_h5ad(path)
            adatas.append(adata)
            names.append(name)
            print(f"Loaded {name}: {adata.n_obs} cells")
        except FileNotFoundError:
            print(f"[SKIP] {name}: not found")

    if len(adatas) < 2:
        print("Need at least 2 datasets for integration")
        return

    # Harmonize genes
    harmonized = harmonize_genes(adatas, names)

    # Concatenate
    adata_combined = ad.concat(harmonized, join="inner")
    print(f"Combined: {adata_combined.n_obs} cells, {adata_combined.n_vars} genes")

    # Preprocess
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    sc.pp.log1p(adata_combined)
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=config["analysis"]["n_top_genes"])
    sc.pp.pca(adata_combined)

    # Integrate
    adata_integrated = integrate_with_harmony(adata_combined)

    # Save
    output_path = DATA_PROCESSED / "integrated_harmony.h5ad"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata_integrated.write_h5ad(output_path)
    print(f"Saved integrated data to {output_path}")

if __name__ == "__main__":
    main()
