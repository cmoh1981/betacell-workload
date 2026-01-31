"""
Pytest Configuration and Fixtures for Drug Discovery Tests
==========================================================
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path for imports
DRUG_DISCOVERY_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(DRUG_DISCOVERY_DIR))


@pytest.fixture
def sample_descriptors_good():
    """Sample descriptors for a drug-like molecule."""
    return {
        'MW': 350,
        'LogP': 2.5,
        'HBD': 2,
        'HBA': 5,
        'TPSA': 80,
        'RotatableBonds': 5,
        'NumRings': 2,
        'NumAromaticRings': 1,
        'FractionCSP3': 0.3,
        'NumHeavyAtoms': 25,
        'QED': 0.7
    }


@pytest.fixture
def sample_descriptors_bad():
    """Sample descriptors for a non-drug-like molecule."""
    return {
        'MW': 650,
        'LogP': 7,
        'HBD': 8,
        'HBA': 15,
        'TPSA': 200,
        'RotatableBonds': 15,
        'NumRings': 1,
        'NumAromaticRings': 0,
        'FractionCSP3': 0.1,
        'NumHeavyAtoms': 50,
        'QED': 0.2
    }


@pytest.fixture
def metformin_smiles():
    """SMILES string for metformin."""
    return "CN(C)C(=N)NC(=N)N"


@pytest.fixture
def dorzagliatin_smiles():
    """SMILES string for dorzagliatin (GCK activator)."""
    return "CC1=CC(=CC=C1)C2=CC(=NC(=N2)NC3=CC=C(C=C3)S(=O)(=O)C)C4=CC=CC=C4O"


@pytest.fixture
def sample_workload_genes():
    """Sample workload genes for testing."""
    return {
        "PDX1": {"category": "capacity", "mr_or": 0.66, "effect": "protective"},
        "GCK": {"category": "capacity", "mr_or": None, "effect": "key_metabolic"},
        "DDIT3": {"category": "stress", "mr_or": None, "effect": "apoptotic"},
        "HSPA5": {"category": "stress", "mr_or": None, "effect": "chaperone"},
    }


@pytest.fixture
def sample_drug_candidates():
    """Sample drug candidates for testing."""
    return {
        "Dorzagliatin": {"targets": ["GCK"], "action": "activator", "status": "approved"},
        "4-PBA": {"targets": ["HSPA5"], "action": "chaperone", "status": "approved"},
        "Salubrinal": {"targets": ["EIF2AK3", "DDIT3"], "action": "modulator", "status": "research"},
    }


@pytest.fixture
def sample_pathways():
    """Sample pathway gene sets for testing."""
    return {
        "KEGG_TYPE_II_DIABETES_MELLITUS": {"INS", "GCK", "PDX1", "ABCC8", "KCNJ11", "SLC2A2"},
        "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": {"HSPA5", "DDIT3", "ATF4", "XBP1", "ERN1"},
        "REACTOME_INSULIN_SECRETION": {"INS", "GCK", "ABCC8", "KCNJ11", "PCSK1"},
    }
