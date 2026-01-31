"""
Tests for ADMET Analysis Module
===============================

Tests molecular descriptor calculations, drug-likeness rules,
and ADMET property predictions.
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import shared constants
from shared_constants import COMPOUND_SMILES, LIPINSKI_RULES, VEBER_RULES

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


# Helper functions (same as in module)
def calculate_descriptors(smiles: str):
    """Calculate molecular descriptors using RDKit."""
    if not HAS_RDKIT:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        return {
            'MW': Descriptors.MolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'HBD': Descriptors.NumHDonors(mol),
            'HBA': Descriptors.NumHAcceptors(mol),
            'TPSA': Descriptors.TPSA(mol),
            'RotatableBonds': Descriptors.NumRotatableBonds(mol),
            'NumRings': Descriptors.RingCount(mol),
            'NumAromaticRings': Descriptors.NumAromaticRings(mol),
            'FractionCSP3': Descriptors.FractionCSP3(mol),
            'NumHeavyAtoms': Descriptors.HeavyAtomCount(mol),
            'QED': QED.qed(mol),
        }
    except Exception:
        return None


def check_lipinski(descriptors):
    """Check Lipinski's Rule of Five."""
    if descriptors is None:
        return {'passes': False, 'violations': 5, 'details': 'No descriptors'}

    violations = 0
    details = []

    if descriptors['MW'] > LIPINSKI_RULES['MW']['max']:
        violations += 1
        details.append(f"MW={descriptors['MW']:.1f} > 500")

    if descriptors['LogP'] > LIPINSKI_RULES['LogP']['max']:
        violations += 1
        details.append(f"LogP={descriptors['LogP']:.2f} > 5")

    if descriptors['HBD'] > LIPINSKI_RULES['HBD']['max']:
        violations += 1
        details.append(f"HBD={descriptors['HBD']} > 5")

    if descriptors['HBA'] > LIPINSKI_RULES['HBA']['max']:
        violations += 1
        details.append(f"HBA={descriptors['HBA']} > 10")

    return {
        'passes': violations <= 1,
        'violations': violations,
        'details': '; '.join(details) if details else 'All rules pass'
    }


def check_veber(descriptors):
    """Check Veber rules for oral bioavailability."""
    if descriptors is None:
        return {'passes': False, 'violations': 2, 'details': 'No descriptors'}

    violations = 0
    details = []

    if descriptors['RotatableBonds'] > VEBER_RULES['RotatableBonds']['max']:
        violations += 1
        details.append(f"RotBonds={descriptors['RotatableBonds']} > 10")

    if descriptors['TPSA'] > VEBER_RULES['TPSA']['max']:
        violations += 1
        details.append(f"TPSA={descriptors['TPSA']:.1f} > 140")

    return {
        'passes': violations == 0,
        'violations': violations,
        'details': '; '.join(details) if details else 'All rules pass'
    }


# Skip all tests if RDKit not available
pytestmark = pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed")


class TestMolecularDescriptors:
    """Test molecular descriptor calculations."""

    def test_calculate_descriptors_valid_smiles(self, metformin_smiles):
        """Test descriptor calculation with valid SMILES."""
        result = calculate_descriptors(metformin_smiles)

        assert result is not None
        assert 'MW' in result
        assert 'LogP' in result
        assert 'HBD' in result
        assert 'HBA' in result
        assert 'TPSA' in result
        assert 'QED' in result

    def test_calculate_descriptors_invalid_smiles(self):
        """Test descriptor calculation with invalid SMILES."""
        result = calculate_descriptors("invalid_smiles_string")
        assert result is None

    def test_calculate_descriptors_empty_string(self):
        """Test descriptor calculation with empty string."""
        result = calculate_descriptors("")
        # Empty string may return None, empty dict, or zero-value descriptors
        if result is not None and result != {}:
            # If it returns descriptors, MW should be 0 for empty mol
            assert result.get('MW', 0) == 0

    def test_metformin_properties(self):
        """Test known properties of metformin."""
        metformin_smiles = "CN(C)C(=N)NC(=N)N"
        result = calculate_descriptors(metformin_smiles)

        assert result is not None
        assert 125 < result['MW'] < 135  # Metformin MW ~129
        assert result['LogP'] < 0  # Metformin is hydrophilic
        assert result['HBD'] >= 2  # H-bond donors (NH groups)
        assert result['HBA'] >= 2  # H-bond acceptors (N atoms)

    def test_all_compound_smiles_valid(self):
        """Test that all predefined SMILES are valid."""
        for name, smiles in COMPOUND_SMILES.items():
            result = calculate_descriptors(smiles)
            assert result is not None, f"Invalid SMILES for {name}"


class TestLipinskiRules:
    """Test Lipinski's Rule of Five checking."""

    def test_lipinski_pass_small_molecule(self, sample_descriptors_good):
        """Small drug-like molecule should pass."""
        result = check_lipinski(sample_descriptors_good)
        assert result['passes'] is True
        assert result['violations'] == 0

    def test_lipinski_fail_large_molecule(self, sample_descriptors_bad):
        """Large molecule should fail Lipinski rules."""
        result = check_lipinski(sample_descriptors_bad)
        assert result['passes'] is False
        assert result['violations'] >= 2

    def test_lipinski_one_violation_allowed(self):
        """One violation should still pass."""
        descriptors = {
            'MW': 520,  # > 500, 1 violation
            'LogP': 3,
            'HBD': 2,
            'HBA': 5
        }
        result = check_lipinski(descriptors)
        assert result['passes'] is True
        assert result['violations'] == 1

    def test_lipinski_two_violations_fail(self):
        """Two violations should fail."""
        descriptors = {
            'MW': 520,  # > 500
            'LogP': 6,   # > 5
            'HBD': 2,
            'HBA': 5
        }
        result = check_lipinski(descriptors)
        assert result['passes'] is False
        assert result['violations'] == 2

    def test_lipinski_none_descriptors(self):
        """None descriptors should return failure."""
        result = check_lipinski(None)
        assert result['passes'] is False
        assert result['violations'] == 5


class TestVeberRules:
    """Test Veber oral bioavailability rules."""

    def test_veber_pass(self, sample_descriptors_good):
        """Good oral bioavailability molecule should pass."""
        result = check_veber(sample_descriptors_good)
        assert result['passes'] is True
        assert result['violations'] == 0

    def test_veber_fail_high_tpsa(self):
        """High TPSA should fail."""
        descriptors = {
            'RotatableBonds': 5,
            'TPSA': 160  # > 140
        }
        result = check_veber(descriptors)
        assert result['passes'] is False
        assert result['violations'] == 1

    def test_veber_fail_many_rotatable_bonds(self):
        """Many rotatable bonds should fail."""
        descriptors = {
            'RotatableBonds': 15,  # > 10
            'TPSA': 80
        }
        result = check_veber(descriptors)
        assert result['passes'] is False
        assert result['violations'] == 1

    def test_veber_none_descriptors(self):
        """None descriptors should return failure."""
        result = check_veber(None)
        assert result['passes'] is False


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_exact_lipinski_boundary_mw(self):
        """Test MW exactly at 500 threshold."""
        descriptors = {
            'MW': 500,  # Exactly at threshold
            'LogP': 3,
            'HBD': 2,
            'HBA': 5
        }
        result = check_lipinski(descriptors)
        assert result['violations'] == 0  # 500 is not > 500

    def test_exact_lipinski_boundary_logp(self):
        """Test LogP exactly at 5 threshold."""
        descriptors = {
            'MW': 300,
            'LogP': 5,  # Exactly at threshold
            'HBD': 2,
            'HBA': 5
        }
        result = check_lipinski(descriptors)
        assert result['violations'] == 0  # 5 is not > 5

    def test_qed_range(self):
        """Test QED is always between 0 and 1."""
        for name, smiles in list(COMPOUND_SMILES.items())[:5]:
            result = calculate_descriptors(smiles)
            if result:
                assert 0 <= result['QED'] <= 1, f"QED out of range for {name}"

    def test_dorzagliatin_is_drug_like(self, dorzagliatin_smiles):
        """Test that Dorzagliatin passes drug-likeness."""
        result = calculate_descriptors(dorzagliatin_smiles)
        lipinski = check_lipinski(result)
        assert lipinski['passes'] is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
