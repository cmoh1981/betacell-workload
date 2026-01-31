#!/usr/bin/env python3
"""
03_workload_scoring.py - Beta-Cell Workload Quantification
Based on three-pillar framework: Biosynthetic Demand, Metabolic Capacity, Stress Response

References:
- Prentki & Nolan, JCI 2006 (PMC5021190)
- Weir & Bonner-Weir, Diabetes 2004
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import yaml

# Load configuration
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

RESULTS = Path(config["paths"]["results"])

# =============================================================================
# WORKLOAD GENE SIGNATURES
# =============================================================================

WORKLOAD_SIGNATURES = {
    # PILLAR 1: Biosynthetic Demand
    "biosynthetic_demand": {
        "insulin_production": ["INS", "IAPP", "PCSK1", "PCSK2", "CPE", "CHGB", "SCG2"],
        "er_folding": ["HSPA5", "HSP90B1", "PDIA4", "PDIA6", "CALR", "CANX", "DNAJB11"],
        "secretory": ["SLC30A8", "SNAP25", "VAMP2", "STX1A", "SYT7", "RAB3A"]
    },

    # PILLAR 2: Metabolic Capacity
    "metabolic_capacity": {
        "glucose_sensing": ["GCK", "SLC2A2", "G6PC2", "PFKFB2", "PFKL"],
        "mitochondrial": ["TFAM", "NRF1", "PPARGC1A", "MT-ND1", "MT-CO1", "MT-ATP6", "HADH"],
        "beta_identity": ["PDX1", "MAFA", "NKX6-1", "UCN3", "NEUROD1", "PAX6", "NKX2-2"],
        "lipid_metabolism": ["SREBF1", "PPARA", "PPARD", "HNF4A", "ACACA", "FASN"]
    },

    # PILLAR 3: Stress Response
    "stress_response": {
        "upr_adaptive": ["XBP1", "ATF6", "ERN1", "EIF2AK3", "HSPA5"],
        "upr_terminal": ["DDIT3", "ATF4", "TRIB3", "BBC3", "GADD45A"],
        "oxidative": ["NFE2L2", "SOD1", "SOD2", "GPX1", "CAT", "TXN", "HMOX1"],
        "inflammatory": ["NFKB1", "TNFAIP3", "IL1B", "CCL2", "CXCL8"],
        "apoptotic": ["BCL2", "BAX", "CASP3", "CASP9", "BID"]
    },

    # Dedifferentiation markers (workload failure)
    "dedifferentiation": {
        "progenitor": ["ALDH1A3", "SOX9", "HES1", "NEUROG3", "MYC"],
        "stress_induced": ["GASTRIN", "REG1A", "REG3A", "CFTR"]
    }
}

# State-specific marker genes
STATE_MARKERS = {
    "resting": {"up": ["UCN3", "MAFA", "SLC2A2"], "down": ["DDIT3", "ATF4"]},
    "active": {"up": ["INS", "GCK", "PCSK1"], "down": ["ALDH1A3"]},
    "stressed": {"up": ["HSPA5", "XBP1", "ATF6"], "down": []},
    "exhausted": {"up": ["DDIT3", "ATF4", "TRIB3"], "down": ["MAFA", "INS"]},
    "failing": {"up": ["ALDH1A3", "GASTRIN"], "down": ["PDX1", "MAFA", "UCN3"]}
}


def compute_module_score(adata, gene_list: list, score_name: str) -> np.ndarray:
    """Compute module score for a gene list."""
    genes_found = [g for g in gene_list if g in adata.var_names]
    if len(genes_found) == 0:
        print(f"  Warning: No genes found for {score_name}")
        return np.zeros(adata.n_obs)

    sc.tl.score_genes(adata, genes_found, score_name=score_name)
    return adata.obs[score_name].values


def compute_pillar_scores(adata):
    """Compute scores for each pillar of the workload framework."""
    print("Computing Pillar 1: Biosynthetic Demand...")
    for module_name, genes in WORKLOAD_SIGNATURES["biosynthetic_demand"].items():
        score_name = f"biosyn_{module_name}"
        compute_module_score(adata, genes, score_name)

    print("Computing Pillar 2: Metabolic Capacity...")
    for module_name, genes in WORKLOAD_SIGNATURES["metabolic_capacity"].items():
        score_name = f"metab_{module_name}"
        compute_module_score(adata, genes, score_name)

    print("Computing Pillar 3: Stress Response...")
    for module_name, genes in WORKLOAD_SIGNATURES["stress_response"].items():
        score_name = f"stress_{module_name}"
        compute_module_score(adata, genes, score_name)

    print("Computing Dedifferentiation markers...")
    for module_name, genes in WORKLOAD_SIGNATURES["dedifferentiation"].items():
        score_name = f"dediff_{module_name}"
        compute_module_score(adata, genes, score_name)

    return adata


def compute_aggregate_scores(adata):
    """Compute aggregate scores for each pillar."""
    # Biosynthetic Demand (aggregate)
    biosyn_cols = [c for c in adata.obs.columns if c.startswith("biosyn_")]
    if biosyn_cols:
        adata.obs["BIOSYNTHETIC_DEMAND"] = adata.obs[biosyn_cols].mean(axis=1)

    # Metabolic Capacity (aggregate)
    metab_cols = [c for c in adata.obs.columns if c.startswith("metab_")]
    if metab_cols:
        adata.obs["METABOLIC_CAPACITY"] = adata.obs[metab_cols].mean(axis=1)

    # Stress Response (aggregate)
    stress_cols = [c for c in adata.obs.columns if c.startswith("stress_")]
    if stress_cols:
        adata.obs["STRESS_RESPONSE"] = adata.obs[stress_cols].mean(axis=1)

    # Dedifferentiation (aggregate)
    dediff_cols = [c for c in adata.obs.columns if c.startswith("dediff_")]
    if dediff_cols:
        adata.obs["DEDIFFERENTIATION"] = adata.obs[dediff_cols].mean(axis=1)

    return adata


def compute_composite_workload_index(adata):
    """
    Compute the Composite Workload Index (CWI).

    CWI = (Demand / Capacity) × Stress_Modifier × (1 + Dediff_Penalty)

    Higher CWI = Higher workload / stress / dysfunction
    """
    demand = adata.obs["BIOSYNTHETIC_DEMAND"].values
    capacity = adata.obs["METABOLIC_CAPACITY"].values
    stress = adata.obs["STRESS_RESPONSE"].values
    dediff = adata.obs["DEDIFFERENTIATION"].values

    # Normalize to positive range
    demand_norm = (demand - demand.min()) / (demand.max() - demand.min() + 1e-10) + 0.1
    capacity_norm = (capacity - capacity.min()) / (capacity.max() - capacity.min() + 1e-10) + 0.1
    stress_norm = (stress - stress.min()) / (stress.max() - stress.min() + 1e-10)
    dediff_norm = (dediff - dediff.min()) / (dediff.max() - dediff.min() + 1e-10)

    # CWI calculation
    cwi = (demand_norm / capacity_norm) * (1 + stress_norm) * (1 + dediff_norm * 0.5)

    # Normalize to interpretable range (0-5)
    cwi_scaled = (cwi - cwi.min()) / (cwi.max() - cwi.min()) * 5

    adata.obs["CWI"] = cwi_scaled
    return adata


def classify_workload_state(adata):
    """Classify cells into workload states based on CWI."""
    cwi = adata.obs["CWI"].values

    states = []
    for val in cwi:
        if val < 0.5:
            states.append("S1_Resting")
        elif val < 1.0:
            states.append("S2_Active")
        elif val < 1.5:
            states.append("S3_Stressed")
        elif val < 2.5:
            states.append("S4_Exhausted")
        else:
            states.append("S5_Failing")

    adata.obs["Workload_State"] = pd.Categorical(
        states,
        categories=["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"],
        ordered=True
    )
    return adata


def analyze_workload_by_condition(adata, condition_col: str = "condition"):
    """Compare workload metrics between conditions (e.g., Normal vs T2D)."""
    if condition_col not in adata.obs.columns:
        print(f"Column {condition_col} not found")
        return None

    results = []
    for condition in adata.obs[condition_col].unique():
        mask = adata.obs[condition_col] == condition
        subset = adata.obs[mask]

        row = {
            "condition": condition,
            "n_cells": mask.sum(),
            "mean_CWI": subset["CWI"].mean(),
            "median_CWI": subset["CWI"].median(),
            "mean_biosynthetic": subset["BIOSYNTHETIC_DEMAND"].mean(),
            "mean_metabolic": subset["METABOLIC_CAPACITY"].mean(),
            "mean_stress": subset["STRESS_RESPONSE"].mean(),
            "mean_dediff": subset["DEDIFFERENTIATION"].mean(),
        }

        # State distribution
        for state in ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]:
            row[f"pct_{state}"] = (subset["Workload_State"] == state).mean() * 100

        results.append(row)

    return pd.DataFrame(results)


def main():
    print("=" * 60)
    print("Beta-Cell Workload Scoring")
    print("=" * 60)

    # Load beta-cell data
    try:
        adata = sc.read_h5ad(config["data_sources"]["segerstolpe"]["path"])
        print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
    except FileNotFoundError:
        print("Data not found. Please check path in config.yaml")
        return

    # Filter to beta cells if cell_type column exists
    if "cell_type" in adata.obs.columns:
        mask = adata.obs["cell_type"].str.contains("beta|Beta", case=False, na=False)
        adata = adata[mask].copy()
        print(f"Filtered to beta cells: {adata.n_obs} cells")

    # Compute all scores
    adata = compute_pillar_scores(adata)
    adata = compute_aggregate_scores(adata)
    adata = compute_composite_workload_index(adata)
    adata = classify_workload_state(adata)

    # Summary statistics
    print("\n" + "=" * 40)
    print("WORKLOAD STATE DISTRIBUTION")
    print("=" * 40)
    print(adata.obs["Workload_State"].value_counts())

    # Compare by condition
    if "condition" in adata.obs.columns:
        comparison = analyze_workload_by_condition(adata)
        print("\n" + "=" * 40)
        print("WORKLOAD BY CONDITION")
        print("=" * 40)
        print(comparison.to_string(index=False))

        # Save comparison
        output_path = RESULTS / "tables" / "workload_by_condition.csv"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        comparison.to_csv(output_path, index=False)

    # Save all scores
    scores_path = RESULTS / "tables" / "workload_scores.csv"
    score_cols = ["CWI", "Workload_State", "BIOSYNTHETIC_DEMAND",
                  "METABOLIC_CAPACITY", "STRESS_RESPONSE", "DEDIFFERENTIATION"]
    score_cols = [c for c in score_cols if c in adata.obs.columns]
    adata.obs[score_cols].to_csv(scores_path)
    print(f"\nSaved workload scores to {scores_path}")


if __name__ == "__main__":
    main()
