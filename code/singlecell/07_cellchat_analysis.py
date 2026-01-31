#!/usr/bin/env python3
"""
07_cellchat_analysis.py - Cell-Cell Communication Analysis
Beta-Cell Workload and Islet Cell Interactions

Implements:
- CellPhoneDB-style ligand-receptor analysis
- Communication network visualization
- Signaling pathway analysis
- Workload-dependent communication changes
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from itertools import product
import warnings
warnings.filterwarnings('ignore')

# Configuration
RESULTS_DIR = Path("results/singlecell/cellchat")
FIGURES_DIR = Path("figures/singlecell/cellchat")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Cell type colors for islets
CELL_TYPE_COLORS = {
    "Beta": "#e74c3c",
    "Alpha": "#3498db",
    "Delta": "#2ecc71",
    "PP": "#9b59b6",
    "Acinar": "#f39c12",
    "Ductal": "#1abc9c",
    "Endothelial": "#e91e63",
    "Immune": "#607d8b"
}

# Workload state colors
STATE_COLORS = {
    "S1_Resting": "#2ecc71",
    "S2_Active": "#3498db",
    "S3_Stressed": "#f39c12",
    "S4_Exhausted": "#e74c3c",
    "S5_Failing": "#8e44ad"
}


def load_ligand_receptor_database():
    """
    Load ligand-receptor interaction database.
    Using a curated subset of CellPhoneDB/CellChat interactions.
    """
    # Curated ligand-receptor pairs relevant to islet biology
    lr_pairs = [
        # Insulin signaling
        {"ligand": "INS", "receptor": "INSR", "pathway": "Insulin", "annotation": "Autocrine/paracrine insulin"},
        {"ligand": "INS", "receptor": "IGF1R", "pathway": "Insulin", "annotation": "Cross-activation"},
        {"ligand": "IGF1", "receptor": "IGF1R", "pathway": "IGF", "annotation": "Growth signaling"},
        {"ligand": "IGF2", "receptor": "IGF1R", "pathway": "IGF", "annotation": "Growth signaling"},

        # Glucagon signaling
        {"ligand": "GCG", "receptor": "GCGR", "pathway": "Glucagon", "annotation": "Counter-regulatory"},
        {"ligand": "GCG", "receptor": "GLP1R", "pathway": "Incretin", "annotation": "GLP1 effect"},

        # Somatostatin
        {"ligand": "SST", "receptor": "SSTR1", "pathway": "Somatostatin", "annotation": "Inhibitory"},
        {"ligand": "SST", "receptor": "SSTR2", "pathway": "Somatostatin", "annotation": "Inhibitory"},
        {"ligand": "SST", "receptor": "SSTR5", "pathway": "Somatostatin", "annotation": "Inhibitory"},

        # Growth factors
        {"ligand": "VEGFA", "receptor": "KDR", "pathway": "VEGF", "annotation": "Angiogenesis"},
        {"ligand": "VEGFA", "receptor": "FLT1", "pathway": "VEGF", "annotation": "Angiogenesis"},
        {"ligand": "EGF", "receptor": "EGFR", "pathway": "EGF", "annotation": "Proliferation"},
        {"ligand": "HGF", "receptor": "MET", "pathway": "HGF", "annotation": "Survival"},
        {"ligand": "FGF2", "receptor": "FGFR1", "pathway": "FGF", "annotation": "Growth"},

        # WNT signaling
        {"ligand": "WNT3A", "receptor": "FZD1", "pathway": "WNT", "annotation": "Beta-cell mass"},
        {"ligand": "WNT4", "receptor": "FZD4", "pathway": "WNT", "annotation": "Maturation"},

        # Notch signaling
        {"ligand": "DLL1", "receptor": "NOTCH1", "pathway": "NOTCH", "annotation": "Differentiation"},
        {"ligand": "JAG1", "receptor": "NOTCH2", "pathway": "NOTCH", "annotation": "Fate decision"},

        # Inflammatory
        {"ligand": "IL1B", "receptor": "IL1R1", "pathway": "IL1", "annotation": "Inflammation"},
        {"ligand": "IL6", "receptor": "IL6R", "pathway": "IL6", "annotation": "Inflammation"},
        {"ligand": "TNF", "receptor": "TNFRSF1A", "pathway": "TNF", "annotation": "Apoptosis"},
        {"ligand": "IFNG", "receptor": "IFNGR1", "pathway": "IFN", "annotation": "Immune"},

        # TGF-beta family
        {"ligand": "TGFB1", "receptor": "TGFBR1", "pathway": "TGFb", "annotation": "Fibrosis"},
        {"ligand": "BMP4", "receptor": "BMPR1A", "pathway": "BMP", "annotation": "Development"},
        {"ligand": "INHBA", "receptor": "ACVR1B", "pathway": "Activin", "annotation": "Beta-cell mass"},

        # Chemokines
        {"ligand": "CXCL12", "receptor": "CXCR4", "pathway": "Chemokine", "annotation": "Migration"},
        {"ligand": "CCL2", "receptor": "CCR2", "pathway": "Chemokine", "annotation": "Immune recruitment"},

        # ECM interactions
        {"ligand": "COL1A1", "receptor": "ITGA1", "pathway": "ECM", "annotation": "Adhesion"},
        {"ligand": "FN1", "receptor": "ITGB1", "pathway": "ECM", "annotation": "Adhesion"},
        {"ligand": "LAMA1", "receptor": "ITGA6", "pathway": "ECM", "annotation": "Basement membrane"},
    ]

    return pd.DataFrame(lr_pairs)


def calculate_communication_scores(adata, lr_database, cell_type_col='cell_type'):
    """
    Calculate ligand-receptor communication scores between cell types.
    """
    print("\n" + "="*60)
    print("Calculating Cell-Cell Communication Scores")
    print("="*60)

    cell_types = adata.obs[cell_type_col].unique()

    # Get mean expression per cell type
    mean_expr = {}
    for ct in cell_types:
        mask = adata.obs[cell_type_col] == ct
        expr = adata[mask].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()
        mean_expr[ct] = pd.Series(expr.mean(axis=0), index=adata.var_names)

    # Calculate LR scores
    results = []
    for _, row in lr_database.iterrows():
        ligand = row['ligand']
        receptor = row['receptor']

        if ligand not in adata.var_names or receptor not in adata.var_names:
            continue

        for source, target in product(cell_types, cell_types):
            lig_expr = mean_expr[source].get(ligand, 0)
            rec_expr = mean_expr[target].get(receptor, 0)

            # Communication score = ligand * receptor expression
            score = np.sqrt(lig_expr * rec_expr)

            results.append({
                'source': source,
                'target': target,
                'ligand': ligand,
                'receptor': receptor,
                'pathway': row['pathway'],
                'annotation': row['annotation'],
                'ligand_expr': lig_expr,
                'receptor_expr': rec_expr,
                'score': score
            })

    comm_df = pd.DataFrame(results)
    print(f"Calculated {len(comm_df)} communication scores")
    print(f"Cell types: {list(cell_types)}")
    print(f"LR pairs analyzed: {comm_df[['ligand', 'receptor']].drop_duplicates().shape[0]}")

    return comm_df


def calculate_workload_dependent_communication(adata, lr_database, state_col='Workload_State'):
    """
    Calculate how communication changes across workload states.
    Focus on beta-cell communications.
    """
    print("\n" + "="*60)
    print("Analyzing Workload-Dependent Communication")
    print("="*60)

    if state_col not in adata.obs.columns:
        print(f"Column {state_col} not found")
        return None

    states = adata.obs[state_col].unique()

    # Calculate mean expression per state
    state_expr = {}
    for state in states:
        mask = adata.obs[state_col] == state
        expr = adata[mask].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()
        state_expr[state] = pd.Series(expr.mean(axis=0), index=adata.var_names)

    # Calculate LR scores per state
    results = []
    for _, row in lr_database.iterrows():
        ligand = row['ligand']
        receptor = row['receptor']

        if ligand not in adata.var_names or receptor not in adata.var_names:
            continue

        for state in states:
            lig_expr = state_expr[state].get(ligand, 0)
            rec_expr = state_expr[state].get(receptor, 0)
            score = np.sqrt(lig_expr * rec_expr)

            results.append({
                'state': state,
                'ligand': ligand,
                'receptor': receptor,
                'pathway': row['pathway'],
                'ligand_expr': lig_expr,
                'receptor_expr': rec_expr,
                'score': score
            })

    state_comm_df = pd.DataFrame(results)

    # Calculate fold changes vs resting state
    if 'S1_Resting' in states:
        baseline = state_comm_df[state_comm_df['state'] == 'S1_Resting'].set_index(['ligand', 'receptor'])['score']
        state_comm_df['fold_change'] = state_comm_df.apply(
            lambda x: x['score'] / (baseline.get((x['ligand'], x['receptor']), 0.001) + 0.001),
            axis=1
        )

    print(f"Analyzed {len(states)} workload states")

    return state_comm_df


def plot_communication_heatmap(comm_df, save_path=None):
    """
    Plot heatmap of cell-cell communication.
    """
    print("\nPlotting communication heatmap...")

    # Aggregate by source-target pairs
    comm_matrix = comm_df.groupby(['source', 'target'])['score'].sum().unstack(fill_value=0)

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(comm_matrix, cmap='Reds', annot=True, fmt='.1f', ax=ax)
    ax.set_title('Cell-Cell Communication Strength')
    ax.set_xlabel('Target Cell Type')
    ax.set_ylabel('Source Cell Type')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_pathway_communication(comm_df, save_path=None):
    """
    Plot communication by signaling pathway.
    """
    print("\nPlotting pathway communication...")

    # Aggregate by pathway
    pathway_scores = comm_df.groupby('pathway')['score'].sum().sort_values(ascending=True)

    fig, ax = plt.subplots(figsize=(10, 8))
    pathway_scores.plot(kind='barh', ax=ax, color='steelblue')
    ax.set_xlabel('Total Communication Score')
    ax.set_ylabel('Signaling Pathway')
    ax.set_title('Communication by Signaling Pathway')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_workload_communication_changes(state_comm_df, save_path=None):
    """
    Plot how communication changes across workload states.
    """
    print("\nPlotting workload-dependent communication changes...")

    if state_comm_df is None:
        return

    # Top changing LR pairs
    if 'fold_change' in state_comm_df.columns:
        # Get max fold change for each LR pair
        max_fc = state_comm_df.groupby(['ligand', 'receptor'])['fold_change'].max()
        top_pairs = max_fc.nlargest(15).index.tolist()

        # Filter to top pairs
        plot_df = state_comm_df[
            state_comm_df.apply(lambda x: (x['ligand'], x['receptor']) in top_pairs, axis=1)
        ]
        plot_df['lr_pair'] = plot_df['ligand'] + '-' + plot_df['receptor']

        fig, ax = plt.subplots(figsize=(14, 8))

        # Pivot for heatmap
        pivot_df = plot_df.pivot_table(
            index='lr_pair', columns='state', values='score', aggfunc='mean'
        )

        # Reorder columns by state progression
        state_order = ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]
        pivot_df = pivot_df[[c for c in state_order if c in pivot_df.columns]]

        sns.heatmap(pivot_df, cmap='RdYlBu_r', annot=True, fmt='.2f', ax=ax)
        ax.set_title('Ligand-Receptor Communication Across Workload States')
        ax.set_xlabel('Workload State')
        ax.set_ylabel('Ligand-Receptor Pair')

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()


def plot_communication_network(comm_df, threshold=0.5, save_path=None):
    """
    Plot communication as a network graph.
    """
    print("\nPlotting communication network...")

    try:
        import networkx as nx
    except ImportError:
        print("NetworkX not installed, skipping network plot")
        return

    # Aggregate and filter
    comm_agg = comm_df.groupby(['source', 'target'])['score'].sum().reset_index()
    comm_agg = comm_agg[comm_agg['score'] > threshold]

    if len(comm_agg) == 0:
        print("No communications above threshold")
        return

    # Create network
    G = nx.DiGraph()

    for _, row in comm_agg.iterrows():
        G.add_edge(row['source'], row['target'], weight=row['score'])

    # Plot
    fig, ax = plt.subplots(figsize=(12, 10))

    pos = nx.spring_layout(G, k=2, iterations=50)

    # Node colors
    node_colors = [CELL_TYPE_COLORS.get(node, '#808080') for node in G.nodes()]

    # Edge widths
    edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
    max_weight = max(edge_weights) if edge_weights else 1
    edge_widths = [3 * w / max_weight for w in edge_weights]

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=2000, alpha=0.8, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax)
    nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.6,
                          edge_color='gray', arrows=True,
                          arrowsize=20, connectionstyle='arc3,rad=0.1', ax=ax)

    ax.set_title('Cell-Cell Communication Network')
    ax.axis('off')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_demo_islet_data():
    """
    Create synthetic islet data with multiple cell types.
    """
    np.random.seed(42)

    cell_types = ['Beta', 'Alpha', 'Delta', 'PP', 'Endothelial']
    n_per_type = [300, 150, 50, 30, 70]
    n_genes = 1500

    all_data = []
    all_obs = []

    for ct, n in zip(cell_types, n_per_type):
        X = np.random.negative_binomial(5, 0.3, size=(n, n_genes)).astype(np.float32)

        obs = pd.DataFrame({
            'cell_type': ct,
            'cell_id': [f"{ct}_{i}" for i in range(n)]
        })

        all_data.append(X)
        all_obs.append(obs)

    X = np.vstack(all_data)
    obs = pd.concat(all_obs, ignore_index=True)

    # Gene names
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    lr_genes = ['INS', 'INSR', 'IGF1R', 'GCG', 'GCGR', 'SST', 'SSTR2',
                'VEGFA', 'KDR', 'IL1B', 'IL1R1', 'TNF', 'TNFRSF1A',
                'TGFB1', 'TGFBR1', 'WNT3A', 'FZD1', 'CXCL12', 'CXCR4']
    for i, gene in enumerate(lr_genes):
        if i < len(gene_names):
            gene_names[i] = gene

    adata = sc.AnnData(X)
    adata.var_names = gene_names
    adata.obs = obs
    adata.obs_names = obs['cell_id']

    # Add workload states for beta cells
    beta_mask = adata.obs['cell_type'] == 'Beta'
    cwi = np.random.beta(2, 5, size=beta_mask.sum())
    adata.obs.loc[beta_mask, 'CWI'] = cwi

    state_order = ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]
    conditions = [cwi < 0.2, (cwi >= 0.2) & (cwi < 0.4), (cwi >= 0.4) & (cwi < 0.6),
                  (cwi >= 0.6) & (cwi < 0.8), cwi >= 0.8]
    adata.obs.loc[beta_mask, 'Workload_State'] = np.select(conditions, state_order, default='S3_Stressed')

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    return adata


def main():
    """Run complete cell-cell communication analysis."""
    print("="*60)
    print("CELL-CELL COMMUNICATION ANALYSIS PIPELINE")
    print("="*60)

    # Load LR database
    lr_database = load_ligand_receptor_database()
    print(f"Loaded {len(lr_database)} ligand-receptor pairs")

    # Load data
    print("\nLoading data...")
    try:
        data_paths = [
            "results/deep_learning/training_data.h5ad",
            "data/processed/islet_cells.h5ad"
        ]
        adata = None
        for path in data_paths:
            if Path(path).exists():
                adata = sc.read_h5ad(path)
                print(f"Loaded from {path}")
                break

        if adata is None:
            print("Creating demo islet data...")
            adata = create_demo_islet_data()

    except Exception as e:
        print(f"Error: {e}")
        adata = create_demo_islet_data()

    print(f"Data shape: {adata.shape}")

    # Determine cell type column
    cell_type_col = None
    for col in ['cell_type', 'celltype', 'Cell_Type', 'Workload_State']:
        if col in adata.obs.columns:
            cell_type_col = col
            break

    if cell_type_col is None:
        print("No cell type column found, using clusters")
        sc.tl.leiden(adata, resolution=0.5)
        cell_type_col = 'leiden'

    # Calculate communication
    comm_df = calculate_communication_scores(adata, lr_database, cell_type_col)

    # Workload-dependent communication
    state_comm_df = None
    if 'Workload_State' in adata.obs.columns:
        beta_adata = adata[adata.obs['cell_type'] == 'Beta'] if 'cell_type' in adata.obs.columns else adata
        state_comm_df = calculate_workload_dependent_communication(beta_adata, lr_database)

    # Generate plots
    plot_communication_heatmap(comm_df, FIGURES_DIR / "communication_heatmap.png")
    plot_pathway_communication(comm_df, FIGURES_DIR / "pathway_communication.png")
    plot_communication_network(comm_df, threshold=0.3, save_path=FIGURES_DIR / "communication_network.png")

    if state_comm_df is not None:
        plot_workload_communication_changes(state_comm_df, FIGURES_DIR / "workload_communication.png")
        state_comm_df.to_csv(RESULTS_DIR / "workload_communication.csv", index=False)

    # Save results
    comm_df.to_csv(RESULTS_DIR / "communication_scores.csv", index=False)
    lr_database.to_csv(RESULTS_DIR / "ligand_receptor_database.csv", index=False)

    print("\n" + "="*60)
    print("CELL-CELL COMMUNICATION ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nResults saved to: {RESULTS_DIR}")
    print(f"Figures saved to: {FIGURES_DIR}")


if __name__ == "__main__":
    main()
