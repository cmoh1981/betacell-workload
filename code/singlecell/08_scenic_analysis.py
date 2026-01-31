#!/usr/bin/env python3
"""
08_scenic_analysis.py - Gene Regulatory Network Analysis
Beta-Cell Transcription Factor Activity Analysis

Implements:
- SCENIC-style regulon inference
- AUCell activity scoring
- TF activity across workload states
- Key beta-cell TF analysis (PDX1, MAFA, NKX6.1)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Configuration
RESULTS_DIR = Path("results/singlecell/scenic")
FIGURES_DIR = Path("figures/singlecell/scenic")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Key beta-cell transcription factors
BETA_CELL_TFS = [
    'PDX1', 'MAFA', 'NKX6-1', 'NKX2-2', 'PAX6', 'NEUROD1',
    'FOXA2', 'HNF1A', 'HNF4A', 'ISL1', 'PAX4', 'RFX6'
]

# Stress-related TFs
STRESS_TFS = [
    'ATF4', 'ATF6', 'XBP1', 'DDIT3', 'CREB3L1', 'NFE2L2',
    'FOXO1', 'FOXO3', 'HIF1A', 'NFKB1', 'RELA', 'JUN'
]

# State colors
STATE_COLORS = {
    "S1_Resting": "#2ecc71",
    "S2_Active": "#3498db",
    "S3_Stressed": "#f39c12",
    "S4_Exhausted": "#e74c3c",
    "S5_Failing": "#8e44ad"
}


def load_tf_target_database():
    """
    Load TF-target gene relationships.
    Using curated beta-cell relevant TF targets.
    """
    # Curated TF-target relationships for beta-cells
    tf_targets = {
        'PDX1': ['INS', 'GCK', 'SLC2A2', 'MAFA', 'NKX6-1', 'IAPP', 'UCN3', 'G6PC2'],
        'MAFA': ['INS', 'GCK', 'SLC2A2', 'PDX1', 'GLUT2', 'SCG5', 'PCSK1'],
        'NKX6-1': ['INS', 'GCK', 'MAFA', 'PDX1', 'UCN3', 'SYT4'],
        'NKX2-2': ['INS', 'MAFA', 'PAX6', 'ARX', 'NEUROD1'],
        'NEUROD1': ['INS', 'GCK', 'PAX6', 'ABCC8', 'SLC30A8'],
        'PAX6': ['INS', 'GCG', 'SST', 'PCSK2', 'MAFA'],
        'FOXA2': ['PDX1', 'HNF4A', 'GCK', 'SLC2A2', 'INS'],
        'HNF1A': ['INS', 'GCK', 'GLUT2', 'HNF4A', 'SLC2A2'],
        'HNF4A': ['HNF1A', 'PDX1', 'INS', 'GCK', 'FOXA2'],
        'ATF4': ['DDIT3', 'ASNS', 'TRIB3', 'HSPA5', 'XBP1'],
        'ATF6': ['HSPA5', 'CALR', 'PDIA4', 'XBP1', 'DDIT3'],
        'XBP1': ['HSPA5', 'DNAJB9', 'EDEM1', 'SEC61A1', 'CALR'],
        'DDIT3': ['GADD34', 'ERO1A', 'BCL2', 'TRIB3', 'ATF4'],
        'FOXO1': ['INS', 'PDX1', 'G6PC', 'PEPCK', 'IGFBP1'],
        'FOXO3': ['SOD2', 'BNIP3', 'BCL6', 'GADD45A', 'ATM'],
        'NFE2L2': ['GCLC', 'NQO1', 'HMOX1', 'GSR', 'TXN'],
    }

    # Convert to DataFrame
    records = []
    for tf, targets in tf_targets.items():
        for target in targets:
            records.append({'TF': tf, 'target': target, 'direction': 'activation'})

    return pd.DataFrame(records)


def calculate_regulon_activity(adata, tf_targets_df):
    """
    Calculate regulon activity using AUCell-like approach.
    Simplified implementation based on target gene expression.
    """
    print("\n" + "="*60)
    print("Calculating Regulon Activity (AUCell-style)")
    print("="*60)

    # Get unique TFs
    tfs = tf_targets_df['TF'].unique()

    # Calculate activity for each TF
    activity_scores = pd.DataFrame(index=adata.obs_names)

    for tf in tfs:
        # Get target genes
        targets = tf_targets_df[tf_targets_df['TF'] == tf]['target'].tolist()
        available_targets = [t for t in targets if t in adata.var_names]

        if len(available_targets) == 0:
            continue

        # Get expression of targets
        target_expr = adata[:, available_targets].X
        if hasattr(target_expr, 'toarray'):
            target_expr = target_expr.toarray()

        # AUCell-like scoring: rank-based enrichment
        # Simplified: mean z-score of target gene expression
        scaler = StandardScaler()
        if target_expr.shape[0] > 1:
            scaled_expr = scaler.fit_transform(target_expr)
            activity = scaled_expr.mean(axis=1)
        else:
            activity = target_expr.mean(axis=1)

        activity_scores[f'{tf}_activity'] = activity

        print(f"  {tf}: {len(available_targets)} targets, mean activity: {activity.mean():.3f}")

    # Add TF expression as additional feature
    for tf in tfs:
        if tf in adata.var_names:
            expr = adata[:, tf].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            else:
                expr = expr.flatten()
            activity_scores[f'{tf}_expression'] = expr

    print(f"\nCalculated activity for {len([c for c in activity_scores.columns if '_activity' in c])} TFs")

    return activity_scores


def analyze_tf_activity_by_workload(adata, activity_scores, state_col='Workload_State'):
    """
    Analyze how TF activity changes across workload states.
    """
    print("\n" + "="*60)
    print("TF Activity Across Workload States")
    print("="*60)

    if state_col not in adata.obs.columns:
        print(f"Column {state_col} not found")
        return None

    # Combine activity with metadata
    combined = pd.concat([activity_scores, adata.obs[[state_col]]], axis=1)

    # Calculate mean activity per state
    activity_cols = [c for c in activity_scores.columns if '_activity' in c]

    state_means = combined.groupby(state_col)[activity_cols].mean()
    state_stds = combined.groupby(state_col)[activity_cols].std()

    # Statistical testing
    stats_results = []
    state_order = ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]
    available_states = [s for s in state_order if s in combined[state_col].values]

    for col in activity_cols:
        tf_name = col.replace('_activity', '')

        # ANOVA across states
        groups = [combined[combined[state_col] == s][col].dropna() for s in available_states]
        groups = [g for g in groups if len(g) > 0]

        if len(groups) >= 2:
            f_stat, p_value = stats.f_oneway(*groups)

            # Trend test (correlation with state order)
            state_numeric = combined[state_col].map({s: i for i, s in enumerate(available_states)})
            valid_idx = ~state_numeric.isna() & ~combined[col].isna()
            if valid_idx.sum() > 10:
                corr, corr_p = stats.spearmanr(state_numeric[valid_idx], combined.loc[valid_idx, col])
            else:
                corr, corr_p = np.nan, np.nan

            stats_results.append({
                'TF': tf_name,
                'anova_f': f_stat,
                'anova_p': p_value,
                'trend_correlation': corr,
                'trend_p': corr_p,
                'mean_resting': state_means.loc[available_states[0], col] if available_states[0] in state_means.index else np.nan,
                'mean_failing': state_means.loc[available_states[-1], col] if available_states[-1] in state_means.index else np.nan
            })

    stats_df = pd.DataFrame(stats_results)
    stats_df['fold_change'] = stats_df['mean_failing'] / (stats_df['mean_resting'] + 0.001)
    stats_df = stats_df.sort_values('anova_p')

    print("\nTop differentially active TFs:")
    print(stats_df.head(10).to_string(index=False))

    return stats_df, state_means


def plot_tf_activity_heatmap(state_means, save_path=None):
    """
    Plot heatmap of TF activity across states.
    """
    print("\nPlotting TF activity heatmap...")

    # Prepare data
    plot_df = state_means.T
    plot_df.index = [idx.replace('_activity', '') for idx in plot_df.index]

    # Reorder columns
    state_order = ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]
    plot_df = plot_df[[c for c in state_order if c in plot_df.columns]]

    # Z-score normalize
    plot_df_z = (plot_df.T - plot_df.T.mean()) / plot_df.T.std()
    plot_df_z = plot_df_z.T

    fig, ax = plt.subplots(figsize=(10, 12))
    sns.heatmap(plot_df_z, cmap='RdBu_r', center=0, annot=False,
                xticklabels=True, yticklabels=True, ax=ax)
    ax.set_title('TF Activity Across Workload States\n(Z-score normalized)')
    ax.set_xlabel('Workload State')
    ax.set_ylabel('Transcription Factor')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.savefig(str(save_path).replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()


def plot_key_tf_trajectories(adata, activity_scores, state_col='Workload_State', save_path=None):
    """
    Plot trajectories of key TFs across workload progression.
    """
    print("\nPlotting key TF trajectories...")

    # Key TFs to highlight
    key_tfs = ['PDX1', 'MAFA', 'NKX6-1', 'ATF4', 'XBP1', 'DDIT3', 'FOXO1']

    activity_cols = [f'{tf}_activity' for tf in key_tfs if f'{tf}_activity' in activity_scores.columns]

    if len(activity_cols) == 0:
        print("No key TFs found in activity scores")
        return

    n_tfs = len(activity_cols)
    n_cols = 3
    n_rows = (n_tfs + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    axes = axes.flatten()

    combined = pd.concat([activity_scores, adata.obs[[state_col]]], axis=1)

    state_order = ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]

    for i, col in enumerate(activity_cols):
        tf_name = col.replace('_activity', '')

        # Violin plot
        plot_data = combined[[state_col, col]].dropna()
        plot_data[state_col] = pd.Categorical(plot_data[state_col], categories=state_order, ordered=True)
        plot_data = plot_data.sort_values(state_col)

        sns.violinplot(data=plot_data, x=state_col, y=col,
                      palette=STATE_COLORS, ax=axes[i])
        axes[i].set_title(f'{tf_name} Activity')
        axes[i].set_xlabel('')
        axes[i].set_ylabel('Activity Score')
        axes[i].tick_params(axis='x', rotation=45)

    # Hide empty subplots
    for j in range(len(activity_cols), len(axes)):
        axes[j].set_visible(False)

    plt.suptitle('Key Transcription Factor Activity Across Beta-Cell Workload States', y=1.02)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.savefig(str(save_path).replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()


def plot_identity_vs_stress_tf(activity_scores, adata, state_col='Workload_State', save_path=None):
    """
    Plot identity TFs vs stress TFs relationship.
    """
    print("\nPlotting identity vs stress TF relationship...")

    # Identity score (mean of PDX1, MAFA, NKX6-1)
    identity_cols = [c for c in activity_scores.columns
                    if any(tf in c for tf in ['PDX1', 'MAFA', 'NKX6']) and '_activity' in c]

    # Stress score (mean of ATF4, XBP1, DDIT3)
    stress_cols = [c for c in activity_scores.columns
                  if any(tf in c for tf in ['ATF4', 'XBP1', 'DDIT3']) and '_activity' in c]

    if len(identity_cols) == 0 or len(stress_cols) == 0:
        print("Not enough TFs for identity/stress comparison")
        return

    identity_score = activity_scores[identity_cols].mean(axis=1)
    stress_score = activity_scores[stress_cols].mean(axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Scatter plot
    if state_col in adata.obs.columns:
        colors = adata.obs[state_col].map(STATE_COLORS)
        scatter = axes[0].scatter(identity_score, stress_score, c=colors, alpha=0.6, s=20)

        # Add legend
        handles = [plt.Line2D([0], [0], marker='o', color='w',
                             markerfacecolor=STATE_COLORS[s], markersize=10, label=s)
                  for s in STATE_COLORS.keys()]
        axes[0].legend(handles=handles, title='Workload State', loc='upper right')
    else:
        axes[0].scatter(identity_score, stress_score, alpha=0.6, s=20)

    # Add correlation
    corr, p = stats.pearsonr(identity_score, stress_score)
    axes[0].set_xlabel('Identity TF Activity\n(PDX1, MAFA, NKX6-1)')
    axes[0].set_ylabel('Stress TF Activity\n(ATF4, XBP1, DDIT3)')
    axes[0].set_title(f'Identity vs Stress TF Activity\n(r={corr:.3f}, p={p:.2e})')

    # Trajectory plot (mean per state)
    if state_col in adata.obs.columns:
        combined = pd.DataFrame({
            'identity': identity_score,
            'stress': stress_score,
            'state': adata.obs[state_col]
        })

        state_order = ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]
        state_means = combined.groupby('state')[['identity', 'stress']].mean()
        state_means = state_means.reindex([s for s in state_order if s in state_means.index])

        # Plot trajectory
        for i, state in enumerate(state_means.index):
            axes[1].scatter(state_means.loc[state, 'identity'],
                          state_means.loc[state, 'stress'],
                          c=STATE_COLORS.get(state, 'gray'), s=200, zorder=5)
            axes[1].annotate(state.replace('S', '').split('_')[0],
                           (state_means.loc[state, 'identity'],
                            state_means.loc[state, 'stress']),
                           fontsize=12, ha='center', va='bottom')

        # Draw trajectory line
        axes[1].plot(state_means['identity'], state_means['stress'],
                    'k--', alpha=0.5, linewidth=2)

        axes[1].set_xlabel('Identity TF Activity')
        axes[1].set_ylabel('Stress TF Activity')
        axes[1].set_title('TF Activity Trajectory\n(State Means)')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.savefig(str(save_path).replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()


def create_demo_data():
    """Create demo data for testing."""
    np.random.seed(42)
    n_cells = 500
    n_genes = 1500

    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.float32)

    # Gene names with TFs and targets
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    key_genes = BETA_CELL_TFS + STRESS_TFS + ['INS', 'GCK', 'SLC2A2', 'HSPA5', 'CALR']
    for i, gene in enumerate(key_genes):
        if i < len(gene_names):
            gene_names[i] = gene

    adata = sc.AnnData(X)
    adata.var_names = gene_names
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]

    # Add workload states
    cwi = np.random.beta(2, 5, size=n_cells)
    adata.obs['CWI'] = cwi

    state_order = ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]
    conditions = [cwi < 0.2, (cwi >= 0.2) & (cwi < 0.4), (cwi >= 0.4) & (cwi < 0.6),
                  (cwi >= 0.6) & (cwi < 0.8), cwi >= 0.8]
    adata.obs['Workload_State'] = np.select(conditions, state_order, default='S3_Stressed')
    adata.obs['Workload_State'] = pd.Categorical(adata.obs['Workload_State'], categories=state_order)

    # Preprocess
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    return adata


def main():
    """Run complete SCENIC-style analysis pipeline."""
    print("="*60)
    print("GENE REGULATORY NETWORK ANALYSIS PIPELINE")
    print("="*60)

    # Load TF-target database
    tf_targets_df = load_tf_target_database()
    print(f"Loaded {len(tf_targets_df)} TF-target relationships")

    # Load data
    print("\nLoading data...")
    try:
        data_paths = [
            "results/deep_learning/training_data.h5ad",
            "results/singlecell/trajectory/trajectory_analyzed.h5ad"
        ]
        adata = None
        for path in data_paths:
            if Path(path).exists():
                adata = sc.read_h5ad(path)
                print(f"Loaded from {path}")
                break

        if adata is None:
            print("Creating demo data...")
            adata = create_demo_data()

    except Exception as e:
        print(f"Error: {e}")
        adata = create_demo_data()

    print(f"Data shape: {adata.shape}")

    # Calculate regulon activity
    activity_scores = calculate_regulon_activity(adata, tf_targets_df)

    # Analyze by workload state
    stats_df, state_means = analyze_tf_activity_by_workload(adata, activity_scores)

    # Generate plots
    if state_means is not None:
        plot_tf_activity_heatmap(state_means, FIGURES_DIR / "tf_activity_heatmap.png")

    plot_key_tf_trajectories(adata, activity_scores, save_path=FIGURES_DIR / "key_tf_trajectories.png")
    plot_identity_vs_stress_tf(activity_scores, adata, save_path=FIGURES_DIR / "identity_vs_stress_tf.png")

    # Save results
    activity_scores.to_csv(RESULTS_DIR / "tf_activity_scores.csv")
    tf_targets_df.to_csv(RESULTS_DIR / "tf_target_database.csv", index=False)
    if stats_df is not None:
        stats_df.to_csv(RESULTS_DIR / "tf_differential_activity.csv", index=False)
    if state_means is not None:
        state_means.to_csv(RESULTS_DIR / "tf_activity_by_state.csv")

    print("\n" + "="*60)
    print("GENE REGULATORY NETWORK ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nResults saved to: {RESULTS_DIR}")
    print(f"Figures saved to: {FIGURES_DIR}")


if __name__ == "__main__":
    main()
