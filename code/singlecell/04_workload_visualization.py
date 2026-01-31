#!/usr/bin/env python3
"""
04_workload_visualization.py - Visualize Beta-Cell Workload States
"""

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
import yaml

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

FIGURES = Path(config["paths"]["figures"])
FIGURES.mkdir(parents=True, exist_ok=True)

# Color palette for workload states
STATE_COLORS = {
    "S1_Resting": "#2ecc71",    # Green
    "S2_Active": "#3498db",      # Blue
    "S3_Stressed": "#f39c12",    # Orange
    "S4_Exhausted": "#e74c3c",   # Red
    "S5_Failing": "#8e44ad"      # Purple
}


def plot_workload_umap(adata, save_path=None):
    """Plot UMAP colored by workload state."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # UMAP by workload state
    sc.pl.umap(adata, color="Workload_State", palette=STATE_COLORS,
               ax=axes[0], show=False, title="Workload States")

    # UMAP by CWI
    sc.pl.umap(adata, color="CWI", cmap="RdYlBu_r",
               ax=axes[1], show=False, title="Composite Workload Index")

    # UMAP by condition
    if "condition" in adata.obs.columns:
        sc.pl.umap(adata, color="condition",
                   ax=axes[2], show=False, title="Condition")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_cwi_distribution(adata, save_path=None):
    """Plot CWI distribution by condition."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Histogram
    if "condition" in adata.obs.columns:
        for cond in adata.obs["condition"].unique():
            mask = adata.obs["condition"] == cond
            axes[0].hist(adata.obs.loc[mask, "CWI"], bins=50, alpha=0.6, label=cond)
        axes[0].legend()
    else:
        axes[0].hist(adata.obs["CWI"], bins=50, alpha=0.7)
    axes[0].set_xlabel("Composite Workload Index")
    axes[0].set_ylabel("Cell Count")
    axes[0].set_title("CWI Distribution")

    # State bar plot
    state_counts = adata.obs["Workload_State"].value_counts()
    colors = [STATE_COLORS[s] for s in state_counts.index]
    axes[1].bar(state_counts.index, state_counts.values, color=colors)
    axes[1].set_xlabel("Workload State")
    axes[1].set_ylabel("Cell Count")
    axes[1].set_title("Cells per Workload State")
    plt.xticks(rotation=45)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_pillar_comparison(adata, save_path=None):
    """Plot pillar scores comparison between conditions."""
    if "condition" not in adata.obs.columns:
        return

    pillars = ["BIOSYNTHETIC_DEMAND", "METABOLIC_CAPACITY", "STRESS_RESPONSE", "DEDIFFERENTIATION"]
    pillars = [p for p in pillars if p in adata.obs.columns]

    fig, axes = plt.subplots(1, len(pillars), figsize=(4 * len(pillars), 4))
    if len(pillars) == 1:
        axes = [axes]

    for i, pillar in enumerate(pillars):
        sns.boxplot(data=adata.obs, x="condition", y=pillar, ax=axes[i])
        axes[i].set_title(pillar.replace("_", " ").title())
        axes[i].set_xlabel("")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_state_markers_heatmap(adata, save_path=None):
    """Plot heatmap of key markers by workload state."""
    markers = ["INS", "UCN3", "MAFA", "PDX1", "GCK", "HSPA5", "XBP1",
               "ATF4", "DDIT3", "ALDH1A3"]
    markers = [m for m in markers if m in adata.var_names]

    if not markers or "Workload_State" not in adata.obs.columns:
        return

    # Compute mean expression per state
    states = adata.obs["Workload_State"].cat.categories
    mean_expr = pd.DataFrame(index=markers, columns=states)

    for state in states:
        mask = adata.obs["Workload_State"] == state
        if mask.sum() > 0:
            expr = adata[mask, markers].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray()
            mean_expr[state] = expr.mean(axis=0)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(mean_expr.astype(float), cmap="RdBu_r", center=0,
                annot=True, fmt=".2f", ax=ax)
    ax.set_title("Marker Expression by Workload State")
    ax.set_xlabel("Workload State")
    ax.set_ylabel("Gene")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    print("=" * 60)
    print("Beta-Cell Workload Visualization")
    print("=" * 60)

    # Load data with workload scores
    try:
        adata = sc.read_h5ad(config["data_sources"]["segerstolpe"]["path"])
    except FileNotFoundError:
        print("Data not found")
        return

    # Check if workload scores exist
    if "CWI" not in adata.obs.columns:
        print("Run 03_workload_scoring.py first")
        return

    # Generate figures
    print("Generating UMAP plots...")
    plot_workload_umap(adata, FIGURES / "workload_umap.png")

    print("Generating CWI distribution...")
    plot_cwi_distribution(adata, FIGURES / "cwi_distribution.png")

    print("Generating pillar comparison...")
    plot_pillar_comparison(adata, FIGURES / "pillar_comparison.png")

    print("Generating marker heatmap...")
    plot_state_markers_heatmap(adata, FIGURES / "state_markers_heatmap.png")

    print(f"\nFigures saved to {FIGURES}")


if __name__ == "__main__":
    main()
