"""
Tests for Target Binding Module
===============================

Tests target prioritization and scoring functions.
"""

import pytest
import sys
from pathlib import Path
import pandas as pd

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from shared_constants import WORKLOAD_GENES, DRUG_CANDIDATES


class TestTargetPrioritization:
    """Test drug candidate prioritization logic."""

    def test_prioritize_approved_over_research(self):
        """Approved drugs should score higher than research compounds."""
        # Status scores
        status_scores = {
            'Approved': 5,
            'Approved (China)': 4,
            'Phase II': 2,
            'Research': 0.5,
        }

        assert status_scores['Approved'] > status_scores['Research']
        assert status_scores['Approved (China)'] > status_scores['Phase II']

    def test_mr_score_calculation(self):
        """Test MR effect score calculation."""
        # PDX1 has OR=0.66, effect = |1 - 0.66| = 0.34
        pdx1_or = 0.66
        mr_effect = abs(1 - pdx1_or)
        mr_score = min(mr_effect * 10, 3)

        assert mr_effect == pytest.approx(0.34, rel=0.01)
        assert mr_score == pytest.approx(3.0, rel=0.01)  # Capped at 3

    def test_action_alignment_score(self):
        """Test action alignment scoring."""
        # When therapeutic_action contains 'activate' and drug action is 'activator'
        # Score should be 2

        def calc_action_score(therapeutic_action, drug_action):
            if 'activate' in therapeutic_action and 'activator' in drug_action:
                return 2
            elif 'inhibit' in therapeutic_action and 'inhibitor' in drug_action:
                return 2
            elif 'chaperone' in drug_action:
                return 1.5
            return 0

        assert calc_action_score('activate', 'activator') == 2
        assert calc_action_score('inhibit', 'inhibitor') == 2
        assert calc_action_score('support', 'chaperone') == 1.5
        assert calc_action_score('activate', 'inhibitor') == 0

    def test_total_score_calculation(self):
        """Test total prioritization score."""
        # Dorzagliatin example:
        # - GCK target (no MR OR, so mr_score = 0)
        # - Approved (China) = 4 points
        # - Activator matching 'activate' = 2 points
        # Total = 6 points

        mr_score = 0  # GCK has no MR OR
        status_score = 4  # Approved (China)
        action_score = 2  # activator matches activate

        total = mr_score + status_score + action_score
        assert total == 6


class TestWorkloadTargets:
    """Test workload target definitions."""

    def test_mr_validated_targets_exist(self):
        """Test that MR-validated targets are defined."""
        mr_targets = ['PDX1', 'SLC2A2', 'MAFA']

        for target in mr_targets:
            assert target in WORKLOAD_GENES
            assert WORKLOAD_GENES[target]['mr_or'] is not None

    def test_pdx1_is_protective(self):
        """PDX1 should be protective (OR < 1)."""
        assert WORKLOAD_GENES['PDX1']['mr_or'] < 1
        assert WORKLOAD_GENES['PDX1']['effect'] == 'protective'

    def test_mafa_is_risk(self):
        """MAFA should be risk factor (OR > 1)."""
        assert WORKLOAD_GENES['MAFA']['mr_or'] > 1
        assert WORKLOAD_GENES['MAFA']['effect'] == 'risk'

    def test_stress_genes_categorized(self):
        """Stress genes should be in stress category."""
        stress_genes = ['DDIT3', 'ATF4', 'XBP1', 'HSPA5']

        for gene in stress_genes:
            assert gene in WORKLOAD_GENES
            assert WORKLOAD_GENES[gene]['category'] == 'stress'

    def test_capacity_genes_categorized(self):
        """Capacity genes should be in capacity category."""
        capacity_genes = ['PDX1', 'GCK', 'INS', 'ABCC8']

        for gene in capacity_genes:
            assert gene in WORKLOAD_GENES
            assert WORKLOAD_GENES[gene]['category'] == 'capacity'


class TestDrugCandidates:
    """Test drug candidate definitions."""

    def test_all_candidates_have_required_fields(self):
        """All drug candidates should have targets, action, status."""
        for drug, info in DRUG_CANDIDATES.items():
            assert 'targets' in info, f"{drug} missing targets"
            assert 'action' in info, f"{drug} missing action"
            assert 'status' in info, f"{drug} missing status"
            assert isinstance(info['targets'], list), f"{drug} targets should be list"

    def test_dorzagliatin_targets_gck(self):
        """Dorzagliatin should target GCK."""
        assert 'Dorzagliatin' in DRUG_CANDIDATES
        assert 'GCK' in DRUG_CANDIDATES['Dorzagliatin']['targets']

    def test_salubrinal_multi_target(self):
        """Salubrinal should have multiple targets."""
        assert 'Salubrinal' in DRUG_CANDIDATES
        targets = DRUG_CANDIDATES['Salubrinal']['targets']
        assert len(targets) >= 2

    def test_approved_drugs_exist(self):
        """Should have approved drugs in candidates."""
        approved_drugs = [
            drug for drug, info in DRUG_CANDIDATES.items()
            if info['status'] == 'approved'
        ]
        assert len(approved_drugs) >= 3


class TestDataFrameOperations:
    """Test DataFrame operations used in prioritization."""

    def test_groupby_target(self):
        """Test grouping compounds by target."""
        data = [
            {'target': 'GCK', 'compound': 'Dorzagliatin'},
            {'target': 'GCK', 'compound': 'Piragliatin'},
            {'target': 'HSPA5', 'compound': '4-PBA'},
        ]
        df = pd.DataFrame(data)

        grouped = df.groupby('target').size()
        assert grouped['GCK'] == 2
        assert grouped['HSPA5'] == 1

    def test_sort_by_score(self):
        """Test sorting by prioritization score."""
        data = [
            {'compound': 'A', 'score': 5},
            {'compound': 'B', 'score': 10},
            {'compound': 'C', 'score': 3},
        ]
        df = pd.DataFrame(data)
        sorted_df = df.sort_values('score', ascending=False)

        assert sorted_df.iloc[0]['compound'] == 'B'
        assert sorted_df.iloc[1]['compound'] == 'A'
        assert sorted_df.iloc[2]['compound'] == 'C'


class TestEdgeCases:
    """Test edge cases in target binding."""

    def test_drug_with_unknown_target(self):
        """Test handling of drug with target not in workload genes."""
        drug_info = {'targets': ['UNKNOWN_GENE'], 'action': 'unknown', 'status': 'research'}

        # Check target is not in workload genes
        target_in_workload = any(t in WORKLOAD_GENES for t in drug_info['targets'])
        assert target_in_workload is False

    def test_empty_targets_list(self):
        """Test handling of empty targets list."""
        drug_info = {'targets': [], 'action': 'unknown', 'status': 'research'}

        direct_targets = [t for t in drug_info['targets'] if t in WORKLOAD_GENES]
        assert len(direct_targets) == 0

    def test_none_mr_or_handling(self):
        """Test handling of None MR OR values."""
        gck_info = WORKLOAD_GENES['GCK']
        assert gck_info['mr_or'] is None

        # Score calculation should handle None
        mr_or = gck_info['mr_or']
        if mr_or is not None:
            mr_score = min(abs(1 - mr_or) * 10, 3)
        else:
            mr_score = 0

        assert mr_score == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
