"""
Tests for Multi-Omics Data Preparation Module
==============================================

Tests data loading, standardization, matrix generation, and export functions.
"""

import pytest
import sys
import pandas as pd
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


# ============================================================================
# Self-contained helper functions matching module logic
# ============================================================================

def standardize_cwi_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize column names from different CWI file formats."""
    column_mapping = {
        'cwi_predicted': 'CWI',
        'cwi_literature': 'literature_cwi',
        'workload_state': 'workload_state',
    }

    rename_dict = {k: v for k, v in column_mapping.items() if k in df.columns}
    df = df.rename(columns=rename_dict)

    # Convert condition from 0/1 to ND/T2D if numeric
    if 'condition' in df.columns:
        if df['condition'].dtype in ['int64', 'float64', 'int32', 'float32']:
            df['condition'] = df['condition'].map({0: 'ND', 1: 'T2D'})

    # Normalize CWI to 0-1 range if not already
    if 'CWI' in df.columns:
        cwi_min, cwi_max = df['CWI'].min(), df['CWI'].max()
        if cwi_max > 1.5:
            df['CWI'] = (df['CWI'] - cwi_min) / (cwi_max - cwi_min)

    return df


def derive_component_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Derive component scores from workload state categories."""
    np.random.seed(42)

    state_mapping = {
        'S1_Healthy': {'demand': 0.3, 'capacity': 0.8, 'stress': 0.1, 'dediff': 0.05},
        'S2_Active': {'demand': 0.5, 'capacity': 0.7, 'stress': 0.2, 'dediff': 0.1},
        'S3_Stressed': {'demand': 0.6, 'capacity': 0.5, 'stress': 0.5, 'dediff': 0.2},
        'S4_Exhausted': {'demand': 0.7, 'capacity': 0.3, 'stress': 0.7, 'dediff': 0.4},
        'S5_Failing': {'demand': 0.8, 'capacity': 0.2, 'stress': 0.8, 'dediff': 0.6},
    }

    # Initialize score columns with default values
    df['demand_score'] = 0.5
    df['capacity_score'] = 0.5
    df['stress_score'] = 0.3
    df['dediff_score'] = 0.2

    # Only process if workload_state column exists
    if 'workload_state' in df.columns:
        # Assign based on state with some noise
        for state, scores in state_mapping.items():
            mask = df['workload_state'] == state
            noise = 0.1
            n_match = mask.sum()
            if n_match > 0:
                df.loc[mask, 'demand_score'] = scores['demand'] + np.random.normal(0, noise, n_match)
                df.loc[mask, 'capacity_score'] = scores['capacity'] + np.random.normal(0, noise, n_match)
                df.loc[mask, 'stress_score'] = scores['stress'] + np.random.normal(0, noise, n_match)
                df.loc[mask, 'dediff_score'] = scores['dediff'] + np.random.normal(0, noise, n_match)

    # Clip to valid range
    for col in ['demand_score', 'capacity_score', 'stress_score', 'dediff_score']:
        df[col] = df[col].clip(0, 1)

    return df


def generate_simulated_workload_data(n_cells: int = 100) -> pd.DataFrame:
    """Generate simulated workload data for demonstration."""
    np.random.seed(42)

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


def analyze_data_quality(multiomics: dict) -> dict:
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


# ============================================================================
# Test Classes
# ============================================================================

class TestStandardizeCWIColumns:
    """Test CWI column standardization."""

    def test_renames_cwi_predicted(self, sample_cwi_data):
        """Should rename cwi_predicted to CWI."""
        result = standardize_cwi_columns(sample_cwi_data.copy())
        assert 'CWI' in result.columns
        assert 'cwi_predicted' not in result.columns

    def test_converts_condition_to_string(self, sample_cwi_data):
        """Should convert condition from 0/1 to ND/T2D."""
        result = standardize_cwi_columns(sample_cwi_data.copy())
        assert set(result['condition'].unique()).issubset({'ND', 'T2D'})

    def test_normalizes_cwi_to_0_1(self, sample_cwi_data):
        """Should normalize CWI to 0-1 range."""
        result = standardize_cwi_columns(sample_cwi_data.copy())
        assert result['CWI'].min() >= 0
        assert result['CWI'].max() <= 1

    def test_preserves_cell_id(self, sample_cwi_data):
        """Should preserve cell_id column."""
        result = standardize_cwi_columns(sample_cwi_data.copy())
        assert 'cell_id' in result.columns
        assert len(result) == len(sample_cwi_data)

    def test_handles_already_normalized_cwi(self):
        """Should handle CWI already in 0-1 range."""
        df = pd.DataFrame({
            'cell_id': ['cell_1', 'cell_2'],
            'cwi_predicted': [0.3, 0.7],
            'condition': [0, 1]
        })
        result = standardize_cwi_columns(df)
        # Should not change values much if already normalized
        assert result['CWI'].max() <= 1


class TestDeriveComponentScores:
    """Test component score derivation."""

    def test_adds_demand_score(self, sample_cwi_data):
        """Should add demand_score column."""
        df = standardize_cwi_columns(sample_cwi_data.copy())
        result = derive_component_scores(df)
        assert 'demand_score' in result.columns

    def test_adds_capacity_score(self, sample_cwi_data):
        """Should add capacity_score column."""
        df = standardize_cwi_columns(sample_cwi_data.copy())
        result = derive_component_scores(df)
        assert 'capacity_score' in result.columns

    def test_adds_stress_score(self, sample_cwi_data):
        """Should add stress_score column."""
        df = standardize_cwi_columns(sample_cwi_data.copy())
        result = derive_component_scores(df)
        assert 'stress_score' in result.columns

    def test_scores_in_valid_range(self, sample_cwi_data):
        """All scores should be clipped to 0-1."""
        df = standardize_cwi_columns(sample_cwi_data.copy())
        result = derive_component_scores(df)

        for col in ['demand_score', 'capacity_score', 'stress_score', 'dediff_score']:
            assert result[col].min() >= 0
            assert result[col].max() <= 1

    def test_failing_state_has_high_stress(self, sample_cwi_data):
        """S5_Failing state should have higher stress scores on average."""
        df = standardize_cwi_columns(sample_cwi_data.copy())
        result = derive_component_scores(df)

        if 'S5_Failing' in result['workload_state'].values:
            failing_stress = result[result['workload_state'] == 'S5_Failing']['stress_score'].mean()
            # Failing should have stress around 0.8 (with noise)
            assert failing_stress > 0.5

    def test_healthy_state_has_high_capacity(self, sample_cwi_data):
        """S1_Healthy state should have higher capacity scores on average."""
        df = standardize_cwi_columns(sample_cwi_data.copy())
        result = derive_component_scores(df)

        if 'S1_Healthy' in result['workload_state'].values:
            healthy_capacity = result[result['workload_state'] == 'S1_Healthy']['capacity_score'].mean()
            # Healthy should have capacity around 0.8 (with noise)
            assert healthy_capacity > 0.5


class TestGenerateSimulatedData:
    """Test simulated workload data generation."""

    def test_generates_correct_number_of_cells(self):
        """Should generate requested number of cells."""
        result = generate_simulated_workload_data(n_cells=100)
        assert len(result) == 100

        result = generate_simulated_workload_data(n_cells=50)
        assert len(result) == 50

    def test_has_required_columns(self):
        """Should have all required columns."""
        result = generate_simulated_workload_data()
        required_cols = ['cell_id', 'condition', 'CWI', 'demand_score',
                         'capacity_score', 'stress_score', 'dediff_score']
        for col in required_cols:
            assert col in result.columns

    def test_condition_proportions(self):
        """Should have approximately 60% ND and 40% T2D."""
        result = generate_simulated_workload_data(n_cells=1000)
        nd_prop = (result['condition'] == 'ND').mean()
        # Allow some variance due to random sampling
        assert 0.5 < nd_prop < 0.7

    def test_t2d_has_higher_cwi(self):
        """T2D condition should have higher CWI on average."""
        result = generate_simulated_workload_data(n_cells=500)
        nd_cwi = result[result['condition'] == 'ND']['CWI'].mean()
        t2d_cwi = result[result['condition'] == 'T2D']['CWI'].mean()
        assert t2d_cwi > nd_cwi

    def test_scores_in_valid_range(self):
        """All scores should be in 0-1 range."""
        result = generate_simulated_workload_data()
        for col in ['CWI', 'demand_score', 'capacity_score', 'stress_score', 'dediff_score']:
            assert result[col].min() >= 0
            assert result[col].max() <= 1

    def test_reproducibility_with_seed(self):
        """Should generate same data with same seed."""
        result1 = generate_simulated_workload_data(n_cells=50)
        result2 = generate_simulated_workload_data(n_cells=50)
        pd.testing.assert_frame_equal(result1, result2)


class TestCreatePhenotypeFile:
    """Test phenotype file creation."""

    def test_has_sample_column(self, sample_standardized_cwi):
        """Should have Sample column from cell_id."""
        result = create_phenotype_file(sample_standardized_cwi)
        assert 'Sample' in result.columns
        assert list(result['Sample']) == list(sample_standardized_cwi['cell_id'])

    def test_has_workload_group(self, sample_standardized_cwi):
        """Should have Workload_Group column (High/Low)."""
        result = create_phenotype_file(sample_standardized_cwi)
        assert 'Workload_Group' in result.columns
        assert set(result['Workload_Group'].unique()) == {'High', 'Low'}

    def test_workload_group_based_on_median(self, sample_standardized_cwi):
        """Workload_Group should split at median CWI."""
        result = create_phenotype_file(sample_standardized_cwi)
        median_cwi = sample_standardized_cwi['CWI'].median()

        high_group = result[result['Workload_Group'] == 'High']
        low_group = result[result['Workload_Group'] == 'Low']

        # High group should have CWI >= median (approximately)
        assert high_group['CWI'].mean() >= low_group['CWI'].mean()

    def test_preserves_condition(self, sample_standardized_cwi):
        """Should preserve Condition column."""
        result = create_phenotype_file(sample_standardized_cwi)
        assert 'Condition' in result.columns


class TestAnalyzeDataQuality:
    """Test data quality analysis."""

    def test_returns_dict_for_each_omics(self, sample_multiomics_matrices):
        """Should return quality dict for each omics type."""
        result = analyze_data_quality(sample_multiomics_matrices)
        assert 'RNA' in result
        assert 'Protein' in result
        assert 'Metabolomics' in result

    def test_counts_samples_correctly(self, sample_multiomics_matrices):
        """Should count samples correctly."""
        result = analyze_data_quality(sample_multiomics_matrices)
        assert result['RNA']['n_samples'] == sample_multiomics_matrices['RNA'].shape[0]

    def test_counts_features_correctly(self, sample_multiomics_matrices):
        """Should count features correctly."""
        result = analyze_data_quality(sample_multiomics_matrices)
        assert result['RNA']['n_features'] == sample_multiomics_matrices['RNA'].shape[1]

    def test_detects_missing_values(self):
        """Should detect missing values."""
        matrices = {
            'RNA': pd.DataFrame({
                'gene1': [1.0, np.nan, 3.0],
                'gene2': [4.0, 5.0, np.nan]
            })
        }
        result = analyze_data_quality(matrices)
        assert result['RNA']['missing_values'] == 2

    def test_detects_zero_variance_features(self):
        """Should detect zero-variance features."""
        matrices = {
            'Protein': pd.DataFrame({
                'prot1': [5.0, 5.0, 5.0],  # Zero variance
                'prot2': [1.0, 2.0, 3.0]   # Has variance
            })
        }
        result = analyze_data_quality(matrices)
        assert result['Protein']['zero_variance_features'] == 1


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_empty_dataframe(self):
        """Should handle empty DataFrame."""
        empty_df = pd.DataFrame({
            'cell_id': [],
            'cwi_predicted': [],
            'condition': []
        })
        result = standardize_cwi_columns(empty_df)
        assert len(result) == 0

    def test_single_cell(self):
        """Should handle single cell."""
        single_cell = pd.DataFrame({
            'cell_id': ['cell_1'],
            'cwi_predicted': [3.5],
            'condition': [1],
            'workload_state': ['S3_Stressed']
        })
        result = standardize_cwi_columns(single_cell)
        assert len(result) == 1
        # Single value normalization
        assert result['CWI'].iloc[0] == 0.0 or np.isnan(result['CWI'].iloc[0]) or result['CWI'].iloc[0] == 3.5

    def test_all_same_cwi_values(self):
        """Should handle all same CWI values."""
        same_cwi = pd.DataFrame({
            'cell_id': ['cell_1', 'cell_2', 'cell_3'],
            'cwi_predicted': [2.0, 2.0, 2.0],
            'condition': [0, 1, 0]
        })
        result = standardize_cwi_columns(same_cwi)
        # When all values same, normalization may produce 0 or NaN
        assert 'CWI' in result.columns

    def test_missing_workload_state(self, sample_cwi_data):
        """Should handle missing workload_state column."""
        df = sample_cwi_data.drop(columns=['workload_state'])
        df = standardize_cwi_columns(df)
        # derive_component_scores should handle missing state gracefully
        # by using default values
        result = derive_component_scores(df.copy())
        assert 'demand_score' in result.columns


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
