#!/usr/bin/env python3
"""
06_rna_velocity.py - RNA Velocity Analysis
Beta-Cell Workload Progression Analysis

Implements:
- scVelo RNA velocity analysis
- Velocity-based pseudotime
- Latent time estimation
- Velocity-based trajectory visualization
"""

import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
RESULTS_DIR = Path("results/singlecell/velocity")
FIGURES_DIR = Path("figures/singlecell/velocity")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Workload state colors
STATE_COLORS = {
    "S1_Resting": "#2ecc71",
    "S2_Active": "#3498db",
    "S3_Stressed": "#f39c12",
    "S4_Exhausted": "#e74c3c",
    "S5_Failing": "#8e44ad"
}

# scVelo settings
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.settings.set_figure_params('scvelo')


def load_and_prepare_data(loom_path=None, h5ad_path=None):
    """
    Load data with spliced/unspliced counts for velocity analysis.
    """
    print("="*60)
    print("Loading Data for RNA Velocity")
    print("="*60)

    if loom_path and Path(loom_path).exists():
        # Load from loom file (contains spliced/unspliced)
        adata = scv.read(loom_path)
        print(f"Loaded from loom: {adata.shape}")
    elif h5ad_path and Path(h5ad_path).exists():
        # Load h5ad and check for velocity layers
        adata = sc.read_h5ad(h5ad_path)
        print(f"Loaded from h5ad: {adata.shape}")
    else:
        # Create demo data
        print("Creating synthetic velocity data for demonstration...")
        adata = create_demo_velocity_data()

    return adata


def preprocess_for_velocity(adata):
    """
    Preprocess data for RNA velocity analysis.
    """
    print("\n" + "="*60)
    print("Preprocessing for Velocity Analysis")
    print("="*60)

    # Filter genes and cells
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

    # Compute moments for velocity estimation
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    print(f"Preprocessed data: {adata.shape}")

    return adata


def run_velocity_analysis(adata, mode='stochastic'):
    """
    Run RNA velocity analysis.

    Args:
        adata: AnnData object with spliced/unspliced layers
        mode: 'deterministic', 'stochastic', or 'dynamical'
    """
    print("\n" + "="*60)
    print(f"Running RNA Velocity ({mode} mode)")
    print("="*60)

    if mode == 'dynamical':
        # Dynamical model (most accurate, but slower)
        scv.tl.recover_dynamics(adata, n_jobs=-1)
        scv.tl.velocity(adata, mode='dynamical')
    else:
        # Stochastic or deterministic
        scv.tl.velocity(adata, mode=mode)

    # Compute velocity graph
    scv.tl.velocity_graph(adata, n_jobs=-1)

    print("Velocity analysis complete")

    return adata


def compute_latent_time(adata):
    """
    Compute latent time from dynamical model.
    """
    print("\n" + "="*60)
    print("Computing Latent Time")
    print("="*60)

    if 'fit_t' not in adata.layers:
        print("Running dynamical model for latent time...")
        scv.tl.recover_dynamics(adata, n_jobs=-1)

    scv.tl.latent_time(adata)

    print(f"Latent time range: {adata.obs['latent_time'].min():.3f} - {adata.obs['latent_time'].max():.3f}")

    return adata


def plot_velocity_results(adata):
    """
    Generate comprehensive velocity visualizations.
    """
    print("\n" + "="*60)
    print("Generating Velocity Plots")
    print("="*60)

    # 1. Velocity stream plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    scv.pl.velocity_embedding_stream(
        adata, basis='umap', color='Workload_State' if 'Workload_State' in adata.obs else None,
        palette=STATE_COLORS if 'Workload_State' in adata.obs else None,
        ax=axes[0, 0], show=False, title='Velocity Streamlines'
    )

    # 2. Velocity arrows
    scv.pl.velocity_embedding(
        adata, basis='umap', arrow_length=3, arrow_size=2,
        color='Workload_State' if 'Workload_State' in adata.obs else None,
        palette=STATE_COLORS if 'Workload_State' in adata.obs else None,
        ax=axes[0, 1], show=False, title='Velocity Arrows'
    )

    # 3. Velocity length (cell dynamism)
    scv.tl.velocity_confidence(adata)
    sc.pl.umap(adata, color='velocity_length', cmap='coolwarm',
               ax=axes[1, 0], show=False, title='Velocity Length (Cell Dynamism)')

    # 4. Velocity confidence
    sc.pl.umap(adata, color='velocity_confidence', cmap='coolwarm',
               ax=axes[1, 1], show=False, title='Velocity Confidence')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "velocity_overview.png", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "velocity_overview.pdf", bbox_inches='tight')
    plt.close()

    # 5. Latent time if available
    if 'latent_time' in adata.obs:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        sc.pl.umap(adata, color='latent_time', cmap='viridis',
                   ax=axes[0], show=False, title='Latent Time')

        if 'Workload_State' in adata.obs:
            sc.pl.umap(adata, color='Workload_State', palette=STATE_COLORS,
                       ax=axes[1], show=False, title='Workload States')

        if 'CWI' in adata.obs:
            sc.pl.umap(adata, color='CWI', cmap='RdYlBu_r',
                       ax=axes[2], show=False, title='CWI')

        plt.tight_layout()
        plt.savefig(FIGURES_DIR / "latent_time.png", dpi=300, bbox_inches='tight')
        plt.savefig(FIGURES_DIR / "latent_time.pdf", bbox_inches='tight')
        plt.close()

    print("Velocity plots saved")


def plot_velocity_genes(adata, genes=None):
    """
    Plot velocity for specific genes.
    """
    print("\n" + "="*60)
    print("Plotting Gene-Specific Velocity")
    print("="*60)

    if genes is None:
        # Key beta-cell genes
        genes = ['INS', 'GCK', 'PDX1', 'MAFA', 'HSPA5', 'DDIT3']

    available_genes = [g for g in genes if g in adata.var_names]

    if len(available_genes) == 0:
        print("No velocity genes found")
        return

    # Phase portraits
    scv.pl.velocity(adata, var_names=available_genes, ncols=3,
                    save=str(FIGURES_DIR / "velocity_phase_portraits.png"))

    # Scatter plots with velocity
    for gene in available_genes:
        try:
            fig, ax = plt.subplots(figsize=(6, 5))
            scv.pl.scatter(adata, gene, color='Workload_State' if 'Workload_State' in adata.obs else None,
                          palette=STATE_COLORS if 'Workload_State' in adata.obs else None,
                          ax=ax, show=False)
            plt.tight_layout()
            plt.savefig(FIGURES_DIR / f"velocity_{gene}.png", dpi=300, bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"Could not plot {gene}: {e}")

    print(f"Gene velocity plots saved for {len(available_genes)} genes")


def analyze_velocity_pseudotime_correlation(adata):
    """
    Analyze correlation between velocity-based metrics and CWI/workload.
    """
    print("\n" + "="*60)
    print("Velocity-Workload Correlation Analysis")
    print("="*60)

    results = {}

    # Correlation with CWI
    if 'CWI' in adata.obs:
        if 'latent_time' in adata.obs:
            corr = np.corrcoef(adata.obs['CWI'], adata.obs['latent_time'])[0, 1]
            results['cwi_latent_time_correlation'] = corr
            print(f"CWI vs Latent Time correlation: r={corr:.3f}")

        if 'velocity_length' in adata.obs:
            corr = np.corrcoef(adata.obs['CWI'], adata.obs['velocity_length'])[0, 1]
            results['cwi_velocity_length_correlation'] = corr
            print(f"CWI vs Velocity Length correlation: r={corr:.3f}")

        if 'velocity_confidence' in adata.obs:
            corr = np.corrcoef(adata.obs['CWI'], adata.obs['velocity_confidence'])[0, 1]
            results['cwi_velocity_confidence_correlation'] = corr
            print(f"CWI vs Velocity Confidence correlation: r={corr:.3f}")

    # Summary by workload state
    if 'Workload_State' in adata.obs:
        state_stats = []
        for state in adata.obs['Workload_State'].unique():
            mask = adata.obs['Workload_State'] == state
            state_data = {
                'state': state,
                'n_cells': mask.sum()
            }
            if 'latent_time' in adata.obs:
                state_data['latent_time_mean'] = adata.obs.loc[mask, 'latent_time'].mean()
            if 'velocity_length' in adata.obs:
                state_data['velocity_length_mean'] = adata.obs.loc[mask, 'velocity_length'].mean()
            state_stats.append(state_data)

        state_df = pd.DataFrame(state_stats)
        state_df.to_csv(RESULTS_DIR / "velocity_by_state.csv", index=False)
        print("\nVelocity metrics by workload state:")
        print(state_df.to_string(index=False))

    # Save results
    pd.DataFrame([results]).T.to_csv(RESULTS_DIR / "velocity_correlations.csv")

    return results


def create_demo_velocity_data():
    """
    Create synthetic velocity data for demonstration.
    """
    np.random.seed(42)
    n_cells = 500
    n_genes = 1000

    # Spliced counts (mature mRNA)
    spliced = np.random.negative_binomial(10, 0.3, size=(n_cells, n_genes)).astype(np.float32)

    # Unspliced counts (nascent mRNA) - correlated but noisier
    unspliced = np.random.negative_binomial(3, 0.3, size=(n_cells, n_genes)).astype(np.float32)

    # Add temporal structure
    pseudotime = np.linspace(0, 1, n_cells)

    # Genes that increase over time (unspliced leads spliced)
    for i in range(100):
        unspliced[:, i] *= (0.3 + 0.7 * pseudotime)
        spliced[:, i] *= (0.2 + 0.5 * np.roll(pseudotime, 20))

    # Genes that decrease (spliced leads unspliced decline)
    for i in range(100, 200):
        spliced[:, i] *= (1 - 0.5 * pseudotime)
        unspliced[:, i] *= (1 - 0.3 * np.roll(pseudotime, -20))

    # Create AnnData
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    marker_genes = ['INS', 'GCK', 'PDX1', 'MAFA', 'UCN3', 'HSPA5', 'XBP1', 'ATF4', 'DDIT3', 'ALDH1A3']
    for i, gene in enumerate(marker_genes):
        if i < len(gene_names):
            gene_names[i] = gene

    adata = sc.AnnData(spliced)
    adata.layers['spliced'] = spliced
    adata.layers['unspliced'] = unspliced
    adata.var_names = gene_names
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]

    # Add CWI
    adata.obs['CWI'] = pseudotime + np.random.normal(0, 0.1, n_cells)
    adata.obs['CWI'] = np.clip(adata.obs['CWI'], 0, 1)

    # Workload states
    state_order = ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]
    conditions = [
        adata.obs['CWI'] < 0.2,
        (adata.obs['CWI'] >= 0.2) & (adata.obs['CWI'] < 0.4),
        (adata.obs['CWI'] >= 0.4) & (adata.obs['CWI'] < 0.6),
        (adata.obs['CWI'] >= 0.6) & (adata.obs['CWI'] < 0.8),
        adata.obs['CWI'] >= 0.8
    ]
    adata.obs['Workload_State'] = np.select(conditions, state_order, default='S3_Stressed')
    adata.obs['Workload_State'] = pd.Categorical(adata.obs['Workload_State'], categories=state_order)

    # Preprocess
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    return adata


def main():
    """Run complete RNA velocity analysis pipeline."""
    print("="*60)
    print("RNA VELOCITY ANALYSIS PIPELINE")
    print("="*60)

    # Load data
    adata = load_and_prepare_data()

    # Preprocess
    adata = preprocess_for_velocity(adata)

    # Run velocity
    adata = run_velocity_analysis(adata, mode='stochastic')

    # Try dynamical model for latent time
    try:
        adata = compute_latent_time(adata)
    except Exception as e:
        print(f"Dynamical model failed: {e}")
        print("Continuing with stochastic velocity...")

    # Plot results
    plot_velocity_results(adata)
    plot_velocity_genes(adata)

    # Analyze correlations
    analyze_velocity_pseudotime_correlation(adata)

    # Save data
    adata.write_h5ad(RESULTS_DIR / "velocity_analyzed.h5ad")

    # Save cell metadata
    velocity_cols = ['velocity_length', 'velocity_confidence']
    if 'latent_time' in adata.obs:
        velocity_cols.append('latent_time')
    if 'CWI' in adata.obs:
        velocity_cols.append('CWI')
    if 'Workload_State' in adata.obs:
        velocity_cols.append('Workload_State')

    available_cols = [c for c in velocity_cols if c in adata.obs.columns]
    adata.obs[available_cols].to_csv(RESULTS_DIR / "cell_velocity_metadata.csv")

    print("\n" + "="*60)
    print("RNA VELOCITY ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nResults saved to: {RESULTS_DIR}")
    print(f"Figures saved to: {FIGURES_DIR}")


if __name__ == "__main__":
    main()
