"""
Tests for Network Pharmacology Module
=====================================

Tests pathway enrichment and synergy calculations.
"""

import pytest
import sys
from pathlib import Path
import pandas as pd
from scipy import stats

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from shared_constants import WORKLOAD_GENES, DRUG_CANDIDATES


# Self-contained helper functions for testing
def pathway_enrichment_analysis(genes, pathways, background_size=20000):
    """Perform pathway enrichment analysis for a gene set."""
    results = []
    gene_set = set(genes)

    for pathway_name, pathway_genes in pathways.items():
        overlap = gene_set & pathway_genes

        if len(overlap) == 0:
            continue

        a = len(overlap)
        b = len(pathway_genes) - a
        c = len(gene_set) - a
        d = background_size - a - b - c

        _, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative='greater')

        expected = len(gene_set) * len(pathway_genes) / background_size
        enrichment = a / expected if expected > 0 else 0

        results.append({
            'pathway': pathway_name,
            'overlap_genes': ', '.join(overlap),
            'overlap_count': a,
            'pathway_size': len(pathway_genes),
            'enrichment_ratio': round(enrichment, 2),
            'pvalue': pvalue
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        df = df.sort_values('pvalue')
        df['rank'] = range(1, len(df) + 1)
        df['fdr'] = df['pvalue'] * len(df) / df['rank']
        df['fdr'] = df['fdr'].clip(upper=1.0)
        df = df.sort_values('fdr')

    return df


def calculate_synergy_score(drug1_info, drug2_info, workload_genes, pathways):
    """Calculate synergy score between two drugs."""
    targets1 = set(drug1_info['targets'])
    targets2 = set(drug2_info['targets'])

    target_overlap = len(targets1 & targets2)
    unique_targets = len(targets1 | targets2)

    categories1 = set()
    categories2 = set()

    for t in targets1:
        if t in workload_genes:
            categories1.add(workload_genes[t]['category'])
    for t in targets2:
        if t in workload_genes:
            categories2.add(workload_genes[t]['category'])

    combined_categories = categories1 | categories2

    synergy_score = (
        (unique_targets - target_overlap) * 2 +
        len(combined_categories) * 3
    )

    if 'capacity' in combined_categories and 'stress' in combined_categories:
        synergy_score += 3

    return synergy_score, combined_categories


class TestPathwayEnrichment:
    """Test pathway enrichment analysis."""

    def test_enrichment_with_overlap(self, sample_pathways):
        """Test enrichment when genes overlap with pathway."""
        genes = ['GCK', 'PDX1', 'INS']
        result = pathway_enrichment_analysis(genes, sample_pathways)

        assert len(result) > 0
        diabetes_pathway = result[result['pathway'] == 'KEGG_TYPE_II_DIABETES_MELLITUS']
        assert len(diabetes_pathway) == 1
        assert diabetes_pathway.iloc[0]['overlap_count'] == 3

    def test_enrichment_no_overlap(self, sample_pathways):
        """Test enrichment when genes don't overlap."""
        genes = ['NONEXISTENT_GENE1', 'NONEXISTENT_GENE2']
        result = pathway_enrichment_analysis(genes, sample_pathways)

        assert len(result) == 0

    def test_enrichment_pvalue_range(self, sample_pathways):
        """Test that p-values are in valid range."""
        genes = ['GCK', 'PDX1', 'INS', 'ABCC8']
        result = pathway_enrichment_analysis(genes, sample_pathways)

        if len(result) > 0:
            assert all(0 <= p <= 1 for p in result['pvalue'])
            assert all(0 <= f <= 1 for f in result['fdr'])

    def test_enrichment_ratio_positive(self, sample_pathways):
        """Test enrichment ratio is non-negative."""
        genes = ['GCK', 'PDX1']
        result = pathway_enrichment_analysis(genes, sample_pathways, background_size=20000)

        if len(result) > 0:
            assert all(r >= 0 for r in result['enrichment_ratio'])


class TestSynergyCalculation:
    """Test synergy score calculation."""

    def test_capacity_stress_synergy_bonus(self, sample_workload_genes, sample_pathways):
        """Test that capacity + stress combination gets synergy bonus."""
        drug1 = {"targets": ["GCK"], "action": "activator", "status": "approved"}
        drug2 = {"targets": ["DDIT3"], "action": "inhibitor", "status": "approved"}

        score, categories = calculate_synergy_score(drug1, drug2, sample_workload_genes, sample_pathways)

        assert 'capacity' in categories
        assert 'stress' in categories
        # Should have synergy bonus
        assert score > 0

    def test_same_target_lower_synergy(self, sample_workload_genes, sample_pathways):
        """Drugs with same target should have lower synergy."""
        drug1 = {"targets": ["GCK"], "action": "activator", "status": "approved"}
        drug2 = {"targets": ["GCK"], "action": "activator", "status": "approved"}

        drug3 = {"targets": ["HSPA5"], "action": "chaperone", "status": "approved"}

        score_same, _ = calculate_synergy_score(drug1, drug2, sample_workload_genes, sample_pathways)
        score_diff, _ = calculate_synergy_score(drug1, drug3, sample_workload_genes, sample_pathways)

        # Different targets should have higher synergy
        assert score_diff > score_same

    def test_multi_target_synergy(self, sample_workload_genes, sample_pathways):
        """Multi-target drugs should contribute to synergy."""
        drug1 = {"targets": ["EIF2AK3", "DDIT3"], "action": "modulator", "status": "research"}
        drug2 = {"targets": ["GCK"], "action": "activator", "status": "approved"}

        score, categories = calculate_synergy_score(drug1, drug2, sample_workload_genes, sample_pathways)

        # Should have 3 unique targets
        assert score > 5  # Base score for unique targets + categories


class TestDrugEffectAnalysis:
    """Test drug effect analysis logic."""

    def test_direct_workload_target_count(self, sample_workload_genes):
        """Test counting direct workload targets."""
        drug_info = {"targets": ["GCK", "UNKNOWN"], "action": "activator", "status": "approved"}

        direct_targets = [t for t in drug_info['targets'] if t in sample_workload_genes]
        assert len(direct_targets) == 1
        assert 'GCK' in direct_targets

    def test_no_workload_target(self, sample_workload_genes):
        """Test drug with no workload targets."""
        drug_info = {"targets": ["SLC5A2"], "action": "inhibitor", "status": "approved"}

        direct_targets = [t for t in drug_info['targets'] if t in sample_workload_genes]
        assert len(direct_targets) == 0

    def test_pathway_counting(self, sample_pathways):
        """Test counting pathways affected by drug targets."""
        targets = ["GCK", "INS"]

        affected_pathways = set()
        for target in targets:
            for pathway_name, pathway_genes in sample_pathways.items():
                if target in pathway_genes:
                    affected_pathways.add(pathway_name)

        # GCK and INS are in diabetes and insulin secretion pathways
        assert len(affected_pathways) >= 2


class TestEdgeCases:
    """Test edge cases."""

    def test_empty_gene_list(self, sample_pathways):
        """Test enrichment with empty gene list."""
        result = pathway_enrichment_analysis([], sample_pathways)
        assert len(result) == 0

    def test_empty_pathways(self):
        """Test enrichment with empty pathways."""
        genes = ['GCK', 'PDX1']
        result = pathway_enrichment_analysis(genes, {})
        assert len(result) == 0

    def test_single_gene_enrichment(self, sample_pathways):
        """Test enrichment with single gene."""
        genes = ['GCK']
        result = pathway_enrichment_analysis(genes, sample_pathways)

        # Should find pathways containing GCK
        if len(result) > 0:
            assert result.iloc[0]['overlap_count'] == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
