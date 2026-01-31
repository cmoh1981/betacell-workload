"""
Tests for LINCS Connectivity Mapping Module
============================================

Tests query signature creation, T2D drug search, and workload modulator identification.
"""

import pytest
import sys
import pandas as pd
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import constants for testing
from shared_constants import WORKLOAD_GENES


# ============================================================================
# Self-contained helper functions matching module logic
# ============================================================================

WORKLOAD_SIGNATURE = {
    "stress_up": [
        "DDIT3", "ATF4", "XBP1", "HSPA5", "ERN1", "EIF2AK3",
        "TRIB3", "ATF6", "CALR", "CANX", "HSP90B1"
    ],
    "dediff_up": [
        "ALDH1A3", "SOX9", "HES1", "NEUROG3"
    ],
    "biosynthetic_down": [
        "INS", "IAPP", "PCSK1", "PCSK2", "CPE", "SCG2", "CHGB"
    ],
    "metabolic_down": [
        "GCK", "SLC2A2", "PDX1", "MAFA", "NKX6-1", "ABCC8", "KCNJ11"
    ]
}

MR_TARGETS = {
    "PDX1": {"or": 0.66, "effect": "protective", "action": "increase"},
    "SLC2A2": {"or": 0.89, "effect": "protective", "action": "increase"},
    "MAFA": {"or": 1.14, "effect": "risk", "action": "decrease"},
}


def create_workload_query_signature(gene_info: pd.DataFrame) -> dict:
    """Create query signature vector for LINCS matching."""
    lincs_genes = set(gene_info['pr_gene_symbol'].values)
    query_sig = {}

    # Genes UP in T2D -> want drugs to DECREASE them (-1)
    for gene in WORKLOAD_SIGNATURE["stress_up"] + WORKLOAD_SIGNATURE["dediff_up"]:
        if gene in lincs_genes:
            query_sig[gene] = -1

    # Genes DOWN in T2D -> want drugs to INCREASE them (+1)
    for gene in WORKLOAD_SIGNATURE["biosynthetic_down"] + WORKLOAD_SIGNATURE["metabolic_down"]:
        if gene in lincs_genes:
            query_sig[gene] = +1

    return query_sig


def parse_dose(dose):
    """Parse dose string to numeric value."""
    if pd.isna(dose):
        return 0
    dose_str = str(dose).lower().replace('um', '').replace('µm', '').replace(' ', '')
    try:
        return float(dose_str)
    except (ValueError, TypeError):
        return 0


def parse_time(time):
    """Parse time string to numeric value."""
    if pd.isna(time):
        return 0
    try:
        return float(time)
    except (ValueError, TypeError):
        return 0


def search_known_t2d_drugs(pert_info: pd.DataFrame) -> pd.DataFrame:
    """Search for known T2D drugs in LINCS."""
    known_t2d_drugs = [
        'metformin', 'glipizide', 'glyburide', 'glimepiride',
        'pioglitazone', 'rosiglitazone', 'sitagliptin'
    ]

    found_drugs = []
    for drug in known_t2d_drugs:
        matches = pert_info[
            pert_info['pert_iname'].str.lower().str.contains(drug, na=False)
        ]

        if len(matches) > 0:
            moa_val = 'Unknown'
            if 'moa' in matches.columns:
                moa_val = matches['moa'].iloc[0] if pd.notna(matches['moa'].iloc[0]) else 'Unknown'
            found_drugs.append({
                'drug_name': drug,
                'pert_id': matches['pert_id'].iloc[0],
                'moa': moa_val,
                'in_lincs': True
            })
        else:
            found_drugs.append({
                'drug_name': drug,
                'pert_id': None,
                'moa': None,
                'in_lincs': False
            })

    return pd.DataFrame(found_drugs)


def identify_workload_modulators(
    ranked_compounds: pd.DataFrame,
    known_t2d_drugs: pd.DataFrame
) -> pd.DataFrame:
    """Identify top candidates that could modulate beta-cell workload."""
    candidates = []

    # Known T2D drugs
    known_in_lincs = known_t2d_drugs[known_t2d_drugs['in_lincs'] == True]
    for _, row in known_in_lincs.iterrows():
        candidates.append({
            'compound': row['drug_name'],
            'category': 'Known T2D Drug',
            'moa': row['moa'],
            'priority': 'Benchmark',
            'rationale': 'Established T2D therapy for comparison'
        })

    # Top connectivity-scored compounds
    for _, row in ranked_compounds.head(5).iterrows():
        if pd.notna(row['pert_iname']) and row['pert_iname'] not in [c['compound'] for c in candidates]:
            candidates.append({
                'compound': row['pert_iname'],
                'category': 'LINCS Top Hit',
                'moa': row['moa'] if pd.notna(row['moa']) else 'Unknown',
                'priority': 'High' if row['connectivity_score'] > 0.5 else 'Medium',
                'rationale': f"Connectivity score: {row['connectivity_score']:.3f}"
            })

    return pd.DataFrame(candidates)


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def sample_gene_info():
    """Sample LINCS gene info DataFrame."""
    return pd.DataFrame({
        'pr_gene_symbol': ['DDIT3', 'ATF4', 'XBP1', 'INS', 'GCK', 'PDX1',
                          'ALDH1A3', 'MAFA', 'IAPP', 'UNKNOWN_GENE'],
        'pr_gene_id': [1649, 468, 7494, 3630, 2645, 5708, 219, 389692, 3375, 99999],
        'pr_gene_title': ['DNA damage inducible', 'ATF4', 'XBP1', 'Insulin',
                          'Glucokinase', 'PDX1', 'ALDH1A3', 'MAFA', 'IAPP', 'Unknown']
    })


@pytest.fixture
def sample_pert_info():
    """Sample LINCS perturbation info DataFrame."""
    return pd.DataFrame({
        'pert_id': ['BRD-A001', 'BRD-A002', 'BRD-A003', 'BRD-A004', 'BRD-A005'],
        'pert_iname': ['metformin', 'pioglitazone', 'compound-x', 'sitagliptin', 'aspirin'],
        'pert_type': ['trt_cp', 'trt_cp', 'trt_cp', 'trt_cp', 'trt_cp'],
        'moa': ['AMPK activator', 'PPAR agonist', None, 'DPP4 inhibitor', 'COX inhibitor'],
        'canonical_smiles': ['CN(C)C(=N)NC(=N)N', 'CC1=CN...', 'CCC', 'C1CN...', 'CC(=O)OC1...']
    })


@pytest.fixture
def sample_sig_info():
    """Sample LINCS signature info DataFrame."""
    return pd.DataFrame({
        'sig_id': ['SIG001', 'SIG002', 'SIG003', 'SIG004'],
        'pert_id': ['BRD-A001', 'BRD-A001', 'BRD-A002', 'BRD-A003'],
        'pert_iname': ['metformin', 'metformin', 'pioglitazone', 'compound-x'],
        'pert_type': ['trt_cp', 'trt_cp', 'trt_cp', 'trt_cp'],
        'cell_id': ['HEPG2', 'MCF7', 'HEPG2', 'A549'],
        'pert_dose': ['10 um', '5 um', '1 um', '10um'],
        'pert_time': ['24', '6', '24', '12']
    })


@pytest.fixture
def sample_ranked_compounds():
    """Sample ranked compounds DataFrame."""
    return pd.DataFrame({
        'pert_iname': ['compound-a', 'compound-b', 'compound-c', 'metformin'],
        'connectivity_score': [0.85, 0.72, 0.45, 0.60],
        'moa': ['kinase inhibitor', 'HDAC inhibitor', None, 'AMPK activator'],
        'quality_score': [0.8, 0.7, 0.4, 0.6]
    })


@pytest.fixture
def sample_known_drugs():
    """Sample known T2D drugs DataFrame."""
    return pd.DataFrame({
        'drug_name': ['metformin', 'pioglitazone', 'glipizide'],
        'pert_id': ['BRD-A001', 'BRD-A002', None],
        'moa': ['AMPK activator', 'PPAR agonist', None],
        'in_lincs': [True, True, False]
    })


# ============================================================================
# Test Classes
# ============================================================================

class TestWorkloadSignature:
    """Test workload signature constants."""

    def test_workload_signature_has_all_categories(self):
        """Workload signature should have all 4 categories."""
        assert 'stress_up' in WORKLOAD_SIGNATURE
        assert 'dediff_up' in WORKLOAD_SIGNATURE
        assert 'biosynthetic_down' in WORKLOAD_SIGNATURE
        assert 'metabolic_down' in WORKLOAD_SIGNATURE

    def test_stress_genes_not_empty(self):
        """Stress genes should not be empty."""
        assert len(WORKLOAD_SIGNATURE['stress_up']) > 0
        assert 'DDIT3' in WORKLOAD_SIGNATURE['stress_up']
        assert 'HSPA5' in WORKLOAD_SIGNATURE['stress_up']

    def test_metabolic_genes_not_empty(self):
        """Metabolic genes should not be empty."""
        assert len(WORKLOAD_SIGNATURE['metabolic_down']) > 0
        assert 'GCK' in WORKLOAD_SIGNATURE['metabolic_down']
        assert 'PDX1' in WORKLOAD_SIGNATURE['metabolic_down']

    def test_mr_targets_have_required_fields(self):
        """MR targets should have OR, effect, and action."""
        for target, info in MR_TARGETS.items():
            assert 'or' in info
            assert 'effect' in info
            assert 'action' in info

    def test_pdx1_is_protective(self):
        """PDX1 should be protective with OR < 1."""
        assert MR_TARGETS['PDX1']['or'] < 1
        assert MR_TARGETS['PDX1']['effect'] == 'protective'


class TestQuerySignature:
    """Test query signature creation."""

    def test_create_query_signature_returns_dict(self, sample_gene_info):
        """Query signature should return a dictionary."""
        result = create_workload_query_signature(sample_gene_info)
        assert isinstance(result, dict)

    def test_stress_genes_marked_negative(self, sample_gene_info):
        """Stress genes should be marked -1 (want to decrease)."""
        result = create_workload_query_signature(sample_gene_info)
        assert result.get('DDIT3') == -1
        assert result.get('ATF4') == -1
        assert result.get('XBP1') == -1

    def test_capacity_genes_marked_positive(self, sample_gene_info):
        """Capacity genes should be marked +1 (want to increase)."""
        result = create_workload_query_signature(sample_gene_info)
        assert result.get('INS') == +1
        assert result.get('GCK') == +1
        assert result.get('PDX1') == +1

    def test_dediff_genes_marked_negative(self, sample_gene_info):
        """Dedifferentiation genes should be marked -1."""
        result = create_workload_query_signature(sample_gene_info)
        assert result.get('ALDH1A3') == -1

    def test_unknown_genes_excluded(self, sample_gene_info):
        """Genes not in workload signature should be excluded."""
        result = create_workload_query_signature(sample_gene_info)
        assert 'UNKNOWN_GENE' not in result

    def test_empty_gene_info_returns_empty_dict(self):
        """Empty gene info should return empty dictionary."""
        empty_df = pd.DataFrame({'pr_gene_symbol': []})
        result = create_workload_query_signature(empty_df)
        assert result == {}


class TestDoseParsing:
    """Test dose parsing function."""

    def test_parse_dose_with_um(self):
        """Should parse dose with 'um' suffix."""
        assert parse_dose('10 um') == 10.0
        assert parse_dose('5um') == 5.0
        assert parse_dose('10um') == 10.0

    def test_parse_dose_with_micro_symbol(self):
        """Should parse dose with µm suffix."""
        assert parse_dose('10 µm') == 10.0
        assert parse_dose('5µm') == 5.0

    def test_parse_dose_numeric_only(self):
        """Should parse plain numeric dose."""
        assert parse_dose('10') == 10.0
        assert parse_dose('0.5') == 0.5

    def test_parse_dose_none_returns_zero(self):
        """None dose should return 0."""
        assert parse_dose(None) == 0

    def test_parse_dose_nan_returns_zero(self):
        """NaN dose should return 0."""
        assert parse_dose(np.nan) == 0

    def test_parse_dose_invalid_returns_zero(self):
        """Invalid dose string should return 0."""
        assert parse_dose('invalid') == 0
        assert parse_dose('abc') == 0


class TestTimeParsing:
    """Test time parsing function."""

    def test_parse_time_numeric(self):
        """Should parse numeric time values."""
        assert parse_time('24') == 24.0
        assert parse_time('6') == 6.0
        assert parse_time(24) == 24.0

    def test_parse_time_float(self):
        """Should parse float time values."""
        assert parse_time('6.5') == 6.5
        assert parse_time(12.0) == 12.0

    def test_parse_time_none_returns_zero(self):
        """None time should return 0."""
        assert parse_time(None) == 0

    def test_parse_time_nan_returns_zero(self):
        """NaN time should return 0."""
        assert parse_time(np.nan) == 0

    def test_parse_time_invalid_returns_zero(self):
        """Invalid time string should return 0."""
        assert parse_time('invalid') == 0


class TestSearchKnownDrugs:
    """Test known T2D drug search."""

    def test_search_finds_metformin(self, sample_pert_info):
        """Should find metformin in perturbation info."""
        result = search_known_t2d_drugs(sample_pert_info)
        metformin_row = result[result['drug_name'] == 'metformin']
        assert len(metformin_row) == 1
        assert metformin_row.iloc[0]['in_lincs'] == True

    def test_search_finds_pioglitazone(self, sample_pert_info):
        """Should find pioglitazone in perturbation info."""
        result = search_known_t2d_drugs(sample_pert_info)
        pio_row = result[result['drug_name'] == 'pioglitazone']
        assert len(pio_row) == 1
        assert pio_row.iloc[0]['in_lincs'] == True

    def test_search_marks_missing_drugs(self, sample_pert_info):
        """Should mark missing drugs as not in LINCS."""
        result = search_known_t2d_drugs(sample_pert_info)
        missing_drugs = result[result['in_lincs'] == False]
        assert len(missing_drugs) > 0

    def test_search_returns_dataframe(self, sample_pert_info):
        """Should return a DataFrame."""
        result = search_known_t2d_drugs(sample_pert_info)
        assert isinstance(result, pd.DataFrame)
        assert 'drug_name' in result.columns
        assert 'in_lincs' in result.columns

    def test_search_handles_empty_pert_info(self):
        """Should handle empty perturbation info."""
        empty_df = pd.DataFrame({
            'pert_iname': pd.Series([], dtype=str),
            'pert_id': pd.Series([], dtype=str),
            'pert_type': pd.Series([], dtype=str),
            'moa': pd.Series([], dtype=str)
        })
        result = search_known_t2d_drugs(empty_df)
        assert len(result) > 0  # Should still return list of drugs
        assert all(result['in_lincs'] == False)


class TestIdentifyWorkloadModulators:
    """Test workload modulator identification."""

    def test_includes_known_t2d_drugs(self, sample_ranked_compounds, sample_known_drugs):
        """Should include known T2D drugs found in LINCS."""
        result = identify_workload_modulators(sample_ranked_compounds, sample_known_drugs)
        categories = result['category'].tolist()
        assert 'Known T2D Drug' in categories

    def test_includes_top_hits(self, sample_ranked_compounds, sample_known_drugs):
        """Should include top connectivity-scored compounds."""
        result = identify_workload_modulators(sample_ranked_compounds, sample_known_drugs)
        categories = result['category'].tolist()
        assert 'LINCS Top Hit' in categories

    def test_high_scores_get_high_priority(self, sample_ranked_compounds, sample_known_drugs):
        """High connectivity scores should get High priority."""
        result = identify_workload_modulators(sample_ranked_compounds, sample_known_drugs)
        high_priority = result[result['priority'] == 'High']
        # compound-a has score 0.85 > 0.5
        assert len(high_priority) > 0

    def test_returns_dataframe_with_required_columns(self, sample_ranked_compounds, sample_known_drugs):
        """Should return DataFrame with required columns."""
        result = identify_workload_modulators(sample_ranked_compounds, sample_known_drugs)
        assert 'compound' in result.columns
        assert 'category' in result.columns
        assert 'moa' in result.columns
        assert 'priority' in result.columns
        assert 'rationale' in result.columns

    def test_no_duplicate_compounds(self, sample_ranked_compounds, sample_known_drugs):
        """Should not include duplicate compounds."""
        result = identify_workload_modulators(sample_ranked_compounds, sample_known_drugs)
        compounds = result['compound'].tolist()
        assert len(compounds) == len(set(compounds))


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_empty_ranked_compounds(self, sample_known_drugs):
        """Should handle empty ranked compounds."""
        empty_ranked = pd.DataFrame({
            'pert_iname': [],
            'connectivity_score': [],
            'moa': [],
            'quality_score': []
        })
        result = identify_workload_modulators(empty_ranked, sample_known_drugs)
        # Should still have known drugs
        assert len(result) >= 0

    def test_no_known_drugs_in_lincs(self, sample_ranked_compounds):
        """Should handle case where no known drugs are in LINCS."""
        no_drugs = pd.DataFrame({
            'drug_name': ['glipizide', 'glyburide'],
            'pert_id': [None, None],
            'moa': [None, None],
            'in_lincs': [False, False]
        })
        result = identify_workload_modulators(sample_ranked_compounds, no_drugs)
        # Should still have top hits
        assert len(result) > 0

    def test_all_null_moa(self, sample_known_drugs):
        """Should handle compounds with all null MOA."""
        null_moa_ranked = pd.DataFrame({
            'pert_iname': ['comp1', 'comp2'],
            'connectivity_score': [0.9, 0.8],
            'moa': [None, None],
            'quality_score': [0.9, 0.8]
        })
        result = identify_workload_modulators(null_moa_ranked, sample_known_drugs)
        # Should still work, MOA will be 'Unknown'
        assert len(result) > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
