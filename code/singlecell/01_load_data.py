#!/usr/bin/env python3
"""
01_load_data.py - Load and inspect single-cell datasets
Beta-Cell T2D Workload Project
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
import yaml

# Load configuration
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Setup paths
DATA_RAW = Path(config["paths"]["data_raw"])
DATA_PROCESSED = Path(config["paths"]["data_processed"])
RESULTS = Path(config["paths"]["results"])

def load_processed_data(dataset_name: str) -> ad.AnnData:
    """Load a processed h5ad file from betacell project."""
    source = config["data_sources"].get(dataset_name)
    if source is None:
        raise ValueError(f"Unknown dataset: {dataset_name}")

    path = Path(source["path"])
    print(f"Loading {dataset_name} from {path}")

    adata = sc.read_h5ad(path)
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")

    return adata

def summarize_dataset(adata: ad.AnnData, name: str) -> dict:
    """Generate summary statistics for a dataset."""
    summary = {
        "name": name,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
    }

    # Check for condition/disease status
    for col in ["condition", "disease", "status", "cell_type"]:
        if col in adata.obs.columns:
            summary[col] = adata.obs[col].value_counts().to_dict()

    return summary

def main():
    """Main function to load and inspect all datasets."""
    print("=" * 60)
    print("Beta-Cell T2D Workload - Data Loading")
    print("=" * 60)

    summaries = []

    for dataset_name in config["data_sources"]:
        try:
            adata = load_processed_data(dataset_name)
            summary = summarize_dataset(adata, dataset_name)
            summaries.append(summary)
            print(f"\n{dataset_name}: {summary['n_cells']:,} cells, {summary['n_genes']:,} genes")
        except FileNotFoundError as e:
            print(f"\n[SKIP] {dataset_name}: File not found")
        except Exception as e:
            print(f"\n[ERROR] {dataset_name}: {e}")

    # Save summary
    summary_df = pd.DataFrame(summaries)
    summary_path = RESULTS / "tables" / "dataset_summary.csv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSummary saved to {summary_path}")

if __name__ == "__main__":
    main()
