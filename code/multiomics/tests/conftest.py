"""
Pytest Configuration and Fixtures for Multiomics Tests
=======================================================

Note: This conftest is specific to multiomics tests.
"""

import pytest
import pandas as pd
import numpy as np


@pytest.fixture
def sample_cwi_data():
    """Sample CWI workload scores DataFrame."""
    np.random.seed(42)
    n = 50
    return pd.DataFrame({
        'cell_id': [f'cell_{i}' for i in range(n)],
        'cwi_predicted': np.random.uniform(0.5, 4.5, n),
        'cwi_literature': np.random.uniform(0.3, 3.0, n),
        'workload_state': np.random.choice(
            ['S1_Healthy', 'S2_Active', 'S3_Stressed', 'S4_Exhausted', 'S5_Failing'],
            n
        ),
        'condition': np.random.choice([0, 1], n, p=[0.6, 0.4])
    })


@pytest.fixture
def sample_standardized_cwi():
    """Sample standardized CWI data (after column standardization)."""
    np.random.seed(42)
    n = 50
    return pd.DataFrame({
        'cell_id': [f'cell_{i}' for i in range(n)],
        'CWI': np.random.uniform(0, 1, n),
        'condition': np.random.choice(['ND', 'T2D'], n, p=[0.6, 0.4]),
        'workload_state': np.random.choice(
            ['S1_Healthy', 'S2_Active', 'S3_Stressed', 'S4_Exhausted', 'S5_Failing'],
            n
        ),
        'demand_score': np.random.uniform(0, 1, n),
        'capacity_score': np.random.uniform(0, 1, n),
        'stress_score': np.random.uniform(0, 1, n),
        'dediff_score': np.random.uniform(0, 1, n)
    })


@pytest.fixture
def sample_multiomics_matrices(sample_standardized_cwi):
    """Sample multi-omics matrices."""
    n = len(sample_standardized_cwi)
    np.random.seed(42)

    return {
        'RNA': pd.DataFrame(
            np.random.normal(5, 2, (n, 20)),
            index=sample_standardized_cwi['cell_id'],
            columns=[f'Gene_{i}' for i in range(20)]
        ),
        'Protein': pd.DataFrame(
            np.random.normal(10, 3, (n, 15)),
            index=sample_standardized_cwi['cell_id'],
            columns=[f'Protein_{i}' for i in range(15)]
        ),
        'Metabolomics': pd.DataFrame(
            np.random.normal(8, 2, (n, 10)),
            index=sample_standardized_cwi['cell_id'],
            columns=[f'Metabolite_{i}' for i in range(10)]
        )
    }
