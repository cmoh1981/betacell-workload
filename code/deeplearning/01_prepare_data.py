#!/usr/bin/env python3
"""
01_prepare_data.py - Prepare data for deep learning workload index
Loads existing h5ad files, generates pseudo-labels from literature signatures
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import warnings
warnings.filterwarnings('ignore')

# Paths - use absolute paths for reliability
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent  # workload/code/deeplearning -> workload
BETACELL_DIR = WORKLOAD_DIR.parent / "betacell"  # ../betacell from workload
RESULTS_DIR = WORKLOAD_DIR / "results" / "deep_learning"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Literature-based gene signatures (from architecture document)
LITERATURE_SIGNATURES = {
    "biosynthetic": {
        "insulin_production": ["INS", "IAPP", "PCSK1", "PCSK2", "CPE", "CHGB", "SCG2"],
        "er_folding": ["HSPA5", "HSP90B1", "PDIA4", "PDIA6", "CALR", "CANX"],
        "secretory": ["SLC30A8", "SNAP25"]
    },
    "metabolic": {
        "glucose_sensing": ["GCK", "SLC2A2", "G6PC2", "PFKFB2"],
        "mitochondrial": ["TFAM", "PPARGC1A", "HADH"],  # MT genes often missing
        "beta_identity": ["PDX1", "MAFA", "NKX6-1", "UCN3", "NEUROD1", "NKX2-2", "PAX6"],
        "lipid": ["PPARA", "PPARD", "HNF4A"]
    },
    "stress": {
        "upr_adaptive": ["XBP1", "ATF6", "ERN1", "EIF2AK3"],
        "upr_terminal": ["DDIT3", "ATF4", "TRIB3", "BBC3", "GADD34"],
        "oxidative": ["NFE2L2", "SOD1", "SOD2", "GPX1", "CAT"],
        "inflammatory": ["NFKB1", "TNF", "IL1B", "CCL2"]
    },
    "dedifferentiation": {
        "progenitor": ["ALDH1A3", "NEUROG3", "SOX9", "HES1", "GASTRIN"],
        "disallowed": ["LDHA", "HK1"]
    }
}


def flatten_signatures(signatures):
    """Flatten nested signature dict to single list."""
    genes = []
    for pillar, modules in signatures.items():
        for module, gene_list in modules.items():
            genes.extend(gene_list)
    return list(set(genes))


def load_processed_data():
    """Load processed h5ad files from betacell project."""
    print("=" * 60)
    print("Loading processed single-cell data")
    print("=" * 60)

    datasets = {}

    # 1. Segerstolpe (primary dataset with T2D labels)
    seg_path = BETACELL_DIR / "data/processed/segerstolpe_processed.h5ad"
    if seg_path.exists():
        print(f"\nLoading Segerstolpe: {seg_path}")
        datasets["segerstolpe"] = sc.read_h5ad(seg_path)
        print(f"  Cells: {datasets['segerstolpe'].n_obs}")
        print(f"  Genes: {datasets['segerstolpe'].n_vars}")
        print(f"  Obs columns: {list(datasets['segerstolpe'].obs.columns)}")

    # 2. GSE221156 Atlas (large dataset) - skip for now due to size/compatibility
    atlas_path = BETACELL_DIR / "data/processed/gse221156_atlas_processed.h5ad"
    if atlas_path.exists():
        print(f"\nSkipping GSE221156 Atlas (large file): {atlas_path}")
        print("  Using Segerstolpe as primary dataset for training")
        # Uncomment below to load atlas if needed (requires compatible pandas/anndata)
        # try:
        #     datasets["atlas"] = sc.read_h5ad(atlas_path)
        #     print(f"  Cells: {datasets['atlas'].n_obs}")
        # except Exception as e:
        #     print(f"  Skipped due to error: {e}")

    # 3. Beta trajectory (for pseudotime)
    traj_path = BETACELL_DIR / "results/celltype_deep_analysis/beta_trajectory.h5ad"
    if traj_path.exists():
        print(f"\nLoading Beta trajectory: {traj_path}")
        datasets["trajectory"] = sc.read_h5ad(traj_path)
        print(f"  Cells: {datasets['trajectory'].n_obs}")

    return datasets


def filter_to_beta_cells(adata, cell_type_col=None):
    """Filter AnnData to beta cells only."""
    # Try common column names
    possible_cols = ["cell_type", "celltype", "CellType", "cell_ontology_class"]

    if cell_type_col is None:
        for col in possible_cols:
            if col in adata.obs.columns:
                cell_type_col = col
                break

    if cell_type_col is None:
        print("  Warning: No cell type column found, returning all cells")
        return adata

    # Filter to beta cells
    beta_patterns = ["beta", "Beta", "BETA", "β"]
    mask = adata.obs[cell_type_col].astype(str).str.contains("|".join(beta_patterns), case=False, na=False)

    n_beta = mask.sum()
    print(f"  Found {n_beta} beta cells out of {adata.n_obs} total")

    return adata[mask].copy()


def score_gene_module(adata, genes, score_name):
    """Score cells for a gene module."""
    genes_found = [g for g in genes if g in adata.var_names]
    if len(genes_found) == 0:
        adata.obs[score_name] = 0
        return 0

    sc.tl.score_genes(adata, genes_found, score_name=score_name)
    return len(genes_found)


def compute_literature_scores(adata):
    """Compute pillar scores based on literature signatures."""
    print("\nComputing literature-based scores...")

    all_genes = flatten_signatures(LITERATURE_SIGNATURES)
    genes_found = [g for g in all_genes if g in adata.var_names]
    print(f"  Literature genes found: {len(genes_found)}/{len(all_genes)}")

    # Score each pillar
    for pillar, modules in LITERATURE_SIGNATURES.items():
        pillar_genes = []
        for module, genes in modules.items():
            pillar_genes.extend(genes)
            n = score_gene_module(adata, genes, f"{pillar}_{module}")

        # Aggregate pillar score
        module_cols = [f"{pillar}_{m}" for m in modules.keys()]
        module_cols = [c for c in module_cols if c in adata.obs.columns]
        if module_cols:
            adata.obs[f"PILLAR_{pillar.upper()}"] = adata.obs[module_cols].mean(axis=1)
            print(f"  {pillar}: {len(pillar_genes)} genes")

    return adata


def compute_composite_workload_index(adata):
    """Compute CWI from pillar scores."""
    print("\nComputing Composite Workload Index (CWI)...")

    # Get pillar scores
    biosyn = adata.obs.get("PILLAR_BIOSYNTHETIC", pd.Series(0, index=adata.obs.index))
    metab = adata.obs.get("PILLAR_METABOLIC", pd.Series(0, index=adata.obs.index))
    stress = adata.obs.get("PILLAR_STRESS", pd.Series(0, index=adata.obs.index))
    dediff = adata.obs.get("PILLAR_DEDIFFERENTIATION", pd.Series(0, index=adata.obs.index))

    # Normalize each pillar to 0-1
    def normalize(x):
        x_min, x_max = x.min(), x.max()
        if x_max - x_min == 0:
            return x * 0
        return (x - x_min) / (x_max - x_min)

    biosyn_norm = normalize(biosyn) + 0.1
    metab_norm = normalize(metab) + 0.1
    stress_norm = normalize(stress)
    dediff_norm = normalize(dediff)

    # CWI = (Demand / Capacity) × (1 + Stress) × (1 + Dediff_penalty)
    cwi = (biosyn_norm / metab_norm) * (1 + stress_norm) * (1 + dediff_norm * 0.5)

    # Normalize to 0-5 range
    cwi_scaled = normalize(cwi) * 5
    adata.obs["CWI_literature"] = cwi_scaled

    print(f"  CWI range: {cwi_scaled.min():.2f} - {cwi_scaled.max():.2f}")
    print(f"  CWI mean: {cwi_scaled.mean():.2f}")

    return adata


def assign_workload_states(adata):
    """Assign cells to workload states based on CWI."""
    cwi = adata.obs["CWI_literature"]

    states = pd.cut(
        cwi,
        bins=[-np.inf, 0.5, 1.0, 1.5, 2.5, np.inf],
        labels=["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]
    )

    adata.obs["workload_state"] = states

    print("\nWorkload state distribution:")
    print(states.value_counts())

    return adata


def prepare_training_data(adata, output_path):
    """Prepare and save training data for deep learning."""
    print("\n" + "=" * 60)
    print("Preparing training data")
    print("=" * 60)

    # Get literature genes + top variable genes
    lit_genes = flatten_signatures(LITERATURE_SIGNATURES)
    lit_genes_found = [g for g in lit_genes if g in adata.var_names]

    # Add highly variable genes
    if "highly_variable" in adata.var.columns:
        hv_genes = adata.var_names[adata.var["highly_variable"]].tolist()
    else:
        # Handle potential inf/nan values before HVG computation
        try:
            # Make a copy and clean data
            adata_clean = adata.copy()
            if hasattr(adata_clean.X, 'toarray'):
                X = adata_clean.X.toarray()
            else:
                X = adata_clean.X.copy()
            # Replace inf with large values and nan with 0
            X = np.nan_to_num(X, nan=0.0, posinf=np.finfo(np.float32).max, neginf=0.0)
            adata_clean.X = X
            sc.pp.highly_variable_genes(adata_clean, n_top_genes=2000, subset=False)
            hv_genes = adata_clean.var_names[adata_clean.var["highly_variable"]].tolist()
        except Exception as e:
            print(f"  Warning: Could not compute HVGs ({e}), using variance-based selection")
            # Fallback: use variance-based selection
            if hasattr(adata.X, 'toarray'):
                X = adata.X.toarray()
            else:
                X = adata.X
            X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
            variances = np.var(X, axis=0)
            top_idx = np.argsort(variances)[-2000:]
            hv_genes = adata.var_names[top_idx].tolist()

    # Combine gene sets
    training_genes = list(set(lit_genes_found + hv_genes[:500]))
    training_genes = [g for g in training_genes if g in adata.var_names]
    print(f"Training genes: {len(training_genes)}")

    # Subset to training genes
    adata_train = adata[:, training_genes].copy()

    # Add condition labels
    if "condition" in adata.obs.columns:
        adata_train.obs["condition_binary"] = (adata.obs["condition"].str.contains("T2D", case=False, na=False)).astype(int)
    elif "disease" in adata.obs.columns:
        adata_train.obs["condition_binary"] = (adata.obs["disease"].str.contains("T2D", case=False, na=False)).astype(int)
    else:
        print("  Warning: No condition column found")
        adata_train.obs["condition_binary"] = 0

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata_train.write_h5ad(output_path)
    print(f"Saved training data to {output_path}")
    print(f"  Cells: {adata_train.n_obs}")
    print(f"  Genes: {adata_train.n_vars}")

    # Save gene list
    gene_df = pd.DataFrame({
        "gene": training_genes,
        "is_literature": [g in lit_genes_found for g in training_genes]
    })
    gene_df.to_csv(output_path.parent / "training_genes.csv", index=False)

    # Save labels
    labels_df = adata_train.obs[["CWI_literature", "workload_state", "condition_binary"]].copy()
    labels_df.to_csv(output_path.parent / "training_labels.csv")

    return adata_train


def main():
    print("=" * 60)
    print("BETA-CELL WORKLOAD - DATA PREPARATION")
    print("=" * 60)

    # Load data
    datasets = load_processed_data()

    if "segerstolpe" not in datasets:
        print("ERROR: Segerstolpe dataset not found")
        return

    # Use Segerstolpe as primary (has T2D labels)
    adata = datasets["segerstolpe"]

    # Filter to beta cells
    print("\nFiltering to beta cells...")
    adata_beta = filter_to_beta_cells(adata)

    if adata_beta.n_obs < 100:
        print("Warning: Very few beta cells found, using all cells")
        adata_beta = adata

    # Compute literature-based scores
    adata_beta = compute_literature_scores(adata_beta)
    adata_beta = compute_composite_workload_index(adata_beta)
    adata_beta = assign_workload_states(adata_beta)

    # Prepare training data
    output_path = RESULTS_DIR / "training_data.h5ad"
    adata_train = prepare_training_data(adata_beta, output_path)

    # Summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total cells for training: {adata_train.n_obs}")
    print(f"Total genes: {adata_train.n_vars}")
    print(f"Literature genes included: {sum(adata_train.var_names.isin(flatten_signatures(LITERATURE_SIGNATURES)))}")

    if "condition_binary" in adata_train.obs:
        print(f"\nCondition distribution:")
        print(f"  T2D: {(adata_train.obs['condition_binary'] == 1).sum()}")
        print(f"  Normal: {(adata_train.obs['condition_binary'] == 0).sum()}")

    print("\nFiles saved:")
    print(f"  - {RESULTS_DIR}/training_data.h5ad")
    print(f"  - {RESULTS_DIR}/training_genes.csv")
    print(f"  - {RESULTS_DIR}/training_labels.csv")


if __name__ == "__main__":
    main()
