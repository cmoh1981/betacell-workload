"""
Multi-Omics Data Preparation for mixOmics Integration
======================================================

Prepares scRNA-seq workload data and multi-omics datasets for
DIABLO analysis in R.

Author: Beta-Cell Workload Analysis Pipeline
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings

warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
RESULTS_DIR = WORKLOAD_DIR / "results" / "multiomics"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_workload_scores() -> pd.DataFrame:
    """Load composite workload index from X-intNMF analysis."""
    # Try multiple possible locations
    possible_paths = [
        WORKLOAD_DIR / "results" / "xintnmf" / "composite_workload_index.csv",
        WORKLOAD_DIR / "results" / "ensemble" / "ensemble_cwi.csv",
        WORKLOAD_DIR / "results" / "deep_learning" / "cwi_scores.csv",
    ]

    for path in possible_paths:
        if path.exists():
            print(f"Loading workload scores from {path}")
            df = pd.read_csv(path)

            # Standardize column names to expected format
            df = standardize_cwi_columns(df)
            return df

    print("No workload scores found, generating simulated data")
    return generate_simulated_workload_data()


def standardize_cwi_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize column names from different CWI file formats."""
    # Map actual columns to expected column names
    column_mapping = {
        'cwi_predicted': 'CWI',
        'cwi_literature': 'literature_cwi',
        'workload_state': 'workload_state',
    }

    # Rename columns that exist
    rename_dict = {k: v for k, v in column_mapping.items() if k in df.columns}
    df = df.rename(columns=rename_dict)

    # Convert condition from 0/1 to ND/T2D if numeric
    if 'condition' in df.columns:
        if df['condition'].dtype in ['int64', 'float64']:
            df['condition'] = df['condition'].map({0: 'ND', 1: 'T2D'})

    # Normalize CWI to 0-1 range if it's not already
    if 'CWI' in df.columns:
        cwi_min, cwi_max = df['CWI'].min(), df['CWI'].max()
        if cwi_max > 1.5:  # Likely not normalized
            df['CWI'] = (df['CWI'] - cwi_min) / (cwi_max - cwi_min)
            print(f"   Normalized CWI from [{cwi_min:.2f}, {cwi_max:.2f}] to [0, 1]")

    # Generate component scores from workload_state if individual scores missing
    if 'workload_state' in df.columns and 'demand_score' not in df.columns:
        df = derive_component_scores(df)

    return df


def derive_component_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Derive component scores from workload state categories."""
    np.random.seed(42)
    n = len(df)

    # Base scores derived from workload state
    state_mapping = {
        'S1_Healthy': {'demand': 0.3, 'capacity': 0.8, 'stress': 0.1, 'dediff': 0.05},
        'S2_Active': {'demand': 0.5, 'capacity': 0.7, 'stress': 0.2, 'dediff': 0.1},
        'S3_Stressed': {'demand': 0.6, 'capacity': 0.5, 'stress': 0.5, 'dediff': 0.2},
        'S4_Exhausted': {'demand': 0.7, 'capacity': 0.3, 'stress': 0.7, 'dediff': 0.4},
        'S5_Failing': {'demand': 0.8, 'capacity': 0.2, 'stress': 0.8, 'dediff': 0.6},
    }

    # Initialize score columns
    df['demand_score'] = 0.5
    df['capacity_score'] = 0.5
    df['stress_score'] = 0.3
    df['dediff_score'] = 0.2

    # Assign based on state with some noise
    for state, scores in state_mapping.items():
        mask = df['workload_state'] == state
        noise = 0.1
        df.loc[mask, 'demand_score'] = scores['demand'] + np.random.normal(0, noise, mask.sum())
        df.loc[mask, 'capacity_score'] = scores['capacity'] + np.random.normal(0, noise, mask.sum())
        df.loc[mask, 'stress_score'] = scores['stress'] + np.random.normal(0, noise, mask.sum())
        df.loc[mask, 'dediff_score'] = scores['dediff'] + np.random.normal(0, noise, mask.sum())

    # Clip to valid range
    for col in ['demand_score', 'capacity_score', 'stress_score', 'dediff_score']:
        df[col] = df[col].clip(0, 1)

    print(f"   Derived component scores from workload states")
    return df


def generate_simulated_workload_data(n_cells: int = 500) -> pd.DataFrame:
    """Generate simulated workload data for demonstration."""
    np.random.seed(42)

    # Create base data
    cwi = pd.DataFrame({
        'cell_id': [f"cell_{i}" for i in range(n_cells)],
        'condition': np.random.choice(['ND', 'T2D'], n_cells, p=[0.6, 0.4]),
        'CWI': np.random.normal(0.5, 0.2, n_cells),
        'demand_score': np.random.normal(0.4, 0.15, n_cells),
        'capacity_score': np.random.normal(0.6, 0.2, n_cells),
        'stress_score': np.random.normal(0.3, 0.15, n_cells),
        'dediff_score': np.random.normal(0.2, 0.1, n_cells)
    })

    # Make T2D have higher workload
    t2d_mask = cwi['condition'] == 'T2D'
    cwi.loc[t2d_mask, 'CWI'] += 0.3
    cwi.loc[t2d_mask, 'stress_score'] += 0.2

    # Clip to valid range
    for col in ['CWI', 'demand_score', 'capacity_score', 'stress_score', 'dediff_score']:
        cwi[col] = cwi[col].clip(0, 1)

    return cwi


def generate_multiomics_matrices(cwi: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Generate multi-omics matrices correlated with workload scores.

    In real analysis, these would be loaded from actual data sources.
    """
    np.random.seed(42)
    n_cells = len(cwi)

    # Workload-related genes (from our analysis)
    workload_genes = {
        'capacity': ['PDX1', 'NKX6-1', 'MAFA', 'GCK', 'SLC2A2', 'INS', 'IAPP',
                    'PCSK1', 'PCSK2', 'CPE', 'ABCC8', 'KCNJ11'],
        'stress': ['DDIT3', 'ATF4', 'XBP1', 'HSPA5', 'ERN1', 'EIF2AK3', 'ATF6',
                  'TRIB3', 'CALR', 'CANX', 'HSP90B1'],
        'dediff': ['ALDH1A3', 'SOX9', 'HES1', 'NEUROG3']
    }

    all_workload_genes = sum(workload_genes.values(), [])
    other_genes = [f"Gene_{i}" for i in range(200)]
    all_genes = all_workload_genes + other_genes

    # RNA-seq matrix
    X_rna = pd.DataFrame(
        np.random.normal(5, 2, (n_cells, len(all_genes))),
        index=cwi['cell_id'],
        columns=all_genes
    )

    # Correlate with workload
    for gene in workload_genes['capacity']:
        if gene in X_rna.columns:
            X_rna[gene] = X_rna[gene] - cwi['CWI'].values * 2

    for gene in workload_genes['stress']:
        if gene in X_rna.columns:
            X_rna[gene] = X_rna[gene] + cwi['stress_score'].values * 3

    for gene in workload_genes['dediff']:
        if gene in X_rna.columns:
            X_rna[gene] = X_rna[gene] + cwi['dediff_score'].values * 2

    # Proteomics matrix (subset of genes + additional proteins)
    protein_genes = all_workload_genes[:20]
    other_proteins = [f"Protein_{i}" for i in range(80)]
    all_proteins = protein_genes + other_proteins

    X_protein = pd.DataFrame(
        np.random.normal(10, 3, (n_cells, len(all_proteins))),
        index=cwi['cell_id'],
        columns=all_proteins
    )

    # Correlate with workload
    for i, gene in enumerate(protein_genes):
        if gene in workload_genes['capacity']:
            X_protein[gene] = X_protein[gene] - cwi['CWI'].values * 1.5
        elif gene in workload_genes['stress']:
            X_protein[gene] = X_protein[gene] + cwi['stress_score'].values * 2

    # Metabolomics matrix
    metabolites = [
        'Glucose', 'Pyruvate', 'Lactate', 'ATP', 'ADP', 'NADH', 'NAD',
        'Glutamate', 'Glutamine', 'Proinsulin', 'C-peptide', 'Palmitate',
        'Oleate', 'Ceramide', 'Sphingomyelin'
    ] + [f"Metabolite_{i}" for i in range(46)]

    X_metab = pd.DataFrame(
        np.random.normal(8, 2, (n_cells, len(metabolites))),
        index=cwi['cell_id'],
        columns=metabolites
    )

    # Correlate key metabolites
    X_metab['Glucose'] = X_metab['Glucose'] + cwi['CWI'].values * 2
    X_metab['ATP'] = X_metab['ATP'] - cwi['stress_score'].values * 2
    X_metab['Proinsulin'] = X_metab['Proinsulin'] + cwi['demand_score'].values * 2
    X_metab['Palmitate'] = X_metab['Palmitate'] + cwi['stress_score'].values * 1.5

    return {
        'RNA': X_rna,
        'Protein': X_protein,
        'Metabolomics': X_metab
    }


def create_phenotype_file(cwi: pd.DataFrame) -> pd.DataFrame:
    """Create phenotype/metadata file for mixOmics."""
    phenotype = pd.DataFrame({
        'Sample': cwi['cell_id'],
        'Condition': cwi['condition'],
        'CWI': cwi['CWI'],
        'Workload_Group': np.where(cwi['CWI'] > cwi['CWI'].median(), 'High', 'Low'),
        'Demand': cwi['demand_score'],
        'Capacity': cwi['capacity_score'],
        'Stress': cwi['stress_score'],
        'Dediff': cwi['dediff_score']
    })
    return phenotype


def export_for_mixomics(
    multiomics: Dict[str, pd.DataFrame],
    phenotype: pd.DataFrame,
    output_dir: Path
) -> None:
    """Export data matrices for R/mixOmics analysis."""
    print("\nExporting data for mixOmics...")

    # Export each omics matrix
    for name, matrix in multiomics.items():
        output_file = output_dir / f"{name.lower()}_matrix.csv"
        matrix.to_csv(output_file)
        print(f"  {name}: {matrix.shape[0]} samples x {matrix.shape[1]} features -> {output_file.name}")

    # Export phenotype
    phenotype_file = output_dir / "phenotype.csv"
    phenotype.to_csv(phenotype_file, index=False)
    print(f"  Phenotype: {len(phenotype)} samples -> {phenotype_file.name}")


def analyze_data_quality(multiomics: Dict[str, pd.DataFrame]) -> Dict:
    """Analyze data quality metrics."""
    quality = {}

    for name, matrix in multiomics.items():
        quality[name] = {
            'n_samples': matrix.shape[0],
            'n_features': matrix.shape[1],
            'missing_values': matrix.isna().sum().sum(),
            'mean_expression': matrix.mean().mean(),
            'std_expression': matrix.std().mean(),
            'zero_variance_features': (matrix.std() == 0).sum()
        }

    return quality


def main():
    """Main data preparation pipeline."""
    print("=" * 60)
    print("MULTI-OMICS DATA PREPARATION")
    print("=" * 60)

    # Load workload scores
    print("\n1. Loading workload scores...")
    cwi = load_workload_scores()
    print(f"   Loaded {len(cwi)} cells")

    # Generate multi-omics data
    print("\n2. Generating multi-omics matrices...")
    multiomics = generate_multiomics_matrices(cwi)

    for name, matrix in multiomics.items():
        print(f"   {name}: {matrix.shape}")

    # Create phenotype file
    print("\n3. Creating phenotype file...")
    phenotype = create_phenotype_file(cwi)
    print(f"   High workload: {(phenotype['Workload_Group'] == 'High').sum()}")
    print(f"   Low workload: {(phenotype['Workload_Group'] == 'Low').sum()}")

    # Analyze data quality
    print("\n4. Data quality analysis...")
    quality = analyze_data_quality(multiomics)

    for name, metrics in quality.items():
        print(f"\n   {name}:")
        print(f"     Samples: {metrics['n_samples']}")
        print(f"     Features: {metrics['n_features']}")
        print(f"     Missing values: {metrics['missing_values']}")
        print(f"     Zero-variance features: {metrics['zero_variance_features']}")

    # Export for mixOmics
    print("\n5. Exporting for mixOmics...")
    export_for_mixomics(multiomics, phenotype, RESULTS_DIR)

    # Save workload scores
    cwi.to_csv(RESULTS_DIR / "workload_scores.csv", index=False)

    # Create summary
    summary = {
        'n_cells': len(cwi),
        'n_t2d': (cwi['condition'] == 'T2D').sum(),
        'n_nd': (cwi['condition'] == 'ND').sum(),
        'rna_features': multiomics['RNA'].shape[1],
        'protein_features': multiomics['Protein'].shape[1],
        'metabolite_features': multiomics['Metabolomics'].shape[1],
        'high_workload_cells': (phenotype['Workload_Group'] == 'High').sum(),
        'low_workload_cells': (phenotype['Workload_Group'] == 'Low').sum()
    }

    pd.DataFrame([summary]).to_csv(RESULTS_DIR / "data_summary.csv", index=False)

    print("\n" + "=" * 60)
    print("DATA PREPARATION COMPLETE")
    print("=" * 60)
    print(f"\nOutput directory: {RESULTS_DIR}")
    print("\nFiles created:")
    print("  - rna_matrix.csv")
    print("  - protein_matrix.csv")
    print("  - metabolomics_matrix.csv")
    print("  - phenotype.csv")
    print("  - workload_scores.csv")
    print("  - data_summary.csv")

    print("\nNext steps:")
    print("  1. Run R script: Rscript 01_mixomics_integration.R")
    print("  2. Or use real multi-omics data by replacing the generate_* functions")

    return multiomics, phenotype


if __name__ == "__main__":
    multiomics, phenotype = main()
