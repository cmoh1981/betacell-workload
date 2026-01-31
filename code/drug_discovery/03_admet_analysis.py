"""
ADMET Analysis for Drug Candidates
==================================

Predicts Absorption, Distribution, Metabolism, Excretion, Toxicity
for workload modulator candidates using RDKit and public models.

Author: Beta-Cell Workload Analysis Pipeline
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
import warnings

warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
RESULTS_DIR = WORKLOAD_DIR / "results" / "drug_discovery"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, AllChem, QED
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("RDKit not installed. Install with: pip install rdkit")


# Drug-likeness rules
LIPINSKI_RULES = {
    'MW': {'max': 500, 'description': 'Molecular Weight <= 500'},
    'LogP': {'max': 5, 'description': 'LogP <= 5'},
    'HBD': {'max': 5, 'description': 'H-bond donors <= 5'},
    'HBA': {'max': 10, 'description': 'H-bond acceptors <= 10'}
}

VEBER_RULES = {
    'RotatableBonds': {'max': 10, 'description': 'Rotatable bonds <= 10'},
    'TPSA': {'max': 140, 'description': 'TPSA <= 140 Å²'}
}


# Known compounds with SMILES for ADMET analysis
COMPOUND_SMILES = {
    # GCK Activators
    'Dorzagliatin': 'CC1=CC(=CC=C1)C2=CC(=NC(=N2)NC3=CC=C(C=C3)S(=O)(=O)C)C4=CC=CC=C4O',
    'Piragliatin': 'CC1=C(C=CC(=C1)C2=CC(=NC(=N2)N)C3=CC=CC=C3)S(=O)(=O)C',

    # ER Stress modulators
    '4-PBA': 'CCCCC1=CC=CC=C1C(=O)O',
    'TUDCA': 'CC(CCC(=O)NCCS(=O)(=O)O)C1CCC2C1(CCC3C2C(CC4C3(CCC(C4)O)C)O)C',
    'Salubrinal': 'CCOC(=O)C(=NNC(=S)NC1=CC=CC=C1Cl)C2=CC=CC=C2',

    # Sulfonylureas
    'Glibenclamide': 'COC1=CC=C(C=C1)CCNC(=O)NS(=O)(=O)C2=CC=C(C=C2)CCNC(=O)C3=CC=C(C=C3)Cl',
    'Glipizide': 'CC1=NC=C(C=C1)C(=O)NCCC2=CC=C(C=C2)S(=O)(=O)NC(=O)NC3CCCCC3',

    # KATP modulators
    'Diazoxide': 'CC1=NC2=CC=C(C=C2S(=O)(=O)N1)Cl',

    # Metformin
    'Metformin': 'CN(C)C(=N)NC(=N)N',

    # SGLT2 inhibitors
    'Dapagliflozin': 'CCOC1=CC=C(C=C1)CC2=C(C=CC(=C2)C3C(C(C(C(O3)CO)O)O)O)Cl',
    'Empagliflozin': 'OCC1OC(C(O)C(O)C1O)C2=CC(=C(OC3=CC=C(C=C3)Cl)C=C2)CC4=CC=CC=C4',
}


def calculate_descriptors(smiles: str) -> Optional[Dict]:
    """Calculate molecular descriptors using RDKit."""
    if not HAS_RDKIT:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        descriptors = {
            # Basic properties
            'MW': Descriptors.MolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'HBD': Descriptors.NumHDonors(mol),
            'HBA': Descriptors.NumHAcceptors(mol),
            'TPSA': Descriptors.TPSA(mol),
            'RotatableBonds': Descriptors.NumRotatableBonds(mol),

            # Additional properties
            'NumRings': Descriptors.RingCount(mol),
            'NumAromaticRings': Descriptors.NumAromaticRings(mol),
            'FractionCSP3': Descriptors.FractionCSP3(mol),
            'NumHeavyAtoms': Descriptors.HeavyAtomCount(mol),

            # Drug-likeness
            'QED': QED.qed(mol),  # Quantitative Estimate of Drug-likeness
        }

        return descriptors

    except Exception as e:
        print(f"Error calculating descriptors: {e}")
        return None


def check_lipinski(descriptors: Dict) -> Dict:
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
        'passes': violations <= 1,  # Allows 1 violation
        'violations': violations,
        'details': '; '.join(details) if details else 'All rules pass'
    }


def check_veber(descriptors: Dict) -> Dict:
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


def predict_admet_properties(smiles: str, compound_name: str) -> Dict:
    """
    Predict ADMET properties for a compound.

    Uses rule-based predictions. For production, would integrate
    with deep learning models (DeepChem, ADMETlab, etc.)
    """
    result = {
        'compound': compound_name,
        'smiles': smiles,
    }

    # Calculate descriptors
    descriptors = calculate_descriptors(smiles)

    if descriptors is None:
        result['status'] = 'failed'
        result['error'] = 'Could not parse SMILES or calculate descriptors'
        return result

    result.update(descriptors)

    # Check drug-likeness rules
    lipinski = check_lipinski(descriptors)
    veber = check_veber(descriptors)

    result['lipinski_passes'] = lipinski['passes']
    result['lipinski_violations'] = lipinski['violations']
    result['lipinski_details'] = lipinski['details']

    result['veber_passes'] = veber['passes']
    result['veber_violations'] = veber['violations']
    result['veber_details'] = veber['details']

    # Simple ADMET predictions based on descriptors
    # These are rough estimates - real predictions need ML models

    # Absorption (based on LogP, TPSA, MW)
    if descriptors['MW'] < 500 and descriptors['LogP'] < 5 and descriptors['TPSA'] < 140:
        result['absorption_prediction'] = 'Good'
    elif descriptors['TPSA'] > 140:
        result['absorption_prediction'] = 'Poor (high TPSA)'
    else:
        result['absorption_prediction'] = 'Moderate'

    # BBB penetration (simplified)
    if descriptors['MW'] < 400 and descriptors['TPSA'] < 90 and descriptors['HBD'] <= 3:
        result['bbb_prediction'] = 'Likely'
    else:
        result['bbb_prediction'] = 'Unlikely'

    # Metabolic stability (rough estimate based on structure)
    if descriptors['FractionCSP3'] > 0.25:
        result['metabolic_stability'] = 'Moderate-Good'
    else:
        result['metabolic_stability'] = 'May be unstable'

    # Overall drug-likeness
    result['drug_likeness'] = 'Good' if lipinski['passes'] and veber['passes'] else 'Moderate' if lipinski['passes'] else 'Poor'

    result['status'] = 'success'

    return result


def analyze_all_candidates() -> pd.DataFrame:
    """Analyze ADMET properties for all candidate compounds."""
    print("=" * 60)
    print("ADMET ANALYSIS FOR DRUG CANDIDATES")
    print("=" * 60)

    if not HAS_RDKIT:
        print("\nWARNING: RDKit not installed. Using limited analysis.")
        # Return basic info without calculations
        basic_results = []
        for name, smiles in COMPOUND_SMILES.items():
            basic_results.append({
                'compound': name,
                'smiles': smiles,
                'status': 'rdkit_not_installed'
            })
        return pd.DataFrame(basic_results)

    results = []

    for compound_name, smiles in COMPOUND_SMILES.items():
        print(f"\nAnalyzing: {compound_name}")
        result = predict_admet_properties(smiles, compound_name)
        results.append(result)

        if result['status'] == 'success':
            print(f"  MW: {result['MW']:.1f}, LogP: {result['LogP']:.2f}, QED: {result['QED']:.3f}")
            print(f"  Lipinski: {'PASS' if result['lipinski_passes'] else 'FAIL'} ({result['lipinski_violations']} violations)")
            print(f"  Drug-likeness: {result['drug_likeness']}")
        else:
            print(f"  Error: {result.get('error', 'Unknown')}")

    return pd.DataFrame(results)


def create_admet_report(results_df: pd.DataFrame) -> str:
    """Create summary ADMET report."""
    lines = [
        "=" * 70,
        "ADMET ANALYSIS REPORT - WORKLOAD MODULATOR CANDIDATES",
        "=" * 70,
        "",
    ]

    if 'drug_likeness' not in results_df.columns:
        lines.append("RDKit not installed - limited analysis available")
        return '\n'.join(lines)

    # Summary statistics
    total = len(results_df)
    good_dl = len(results_df[results_df['drug_likeness'] == 'Good'])
    lipinski_pass = results_df['lipinski_passes'].sum()

    lines.extend([
        "SUMMARY",
        "-" * 40,
        f"Total compounds analyzed: {total}",
        f"Good drug-likeness: {good_dl} ({good_dl/total*100:.1f}%)",
        f"Lipinski compliant: {lipinski_pass} ({lipinski_pass/total*100:.1f}%)",
        "",
        "COMPOUND DETAILS",
        "-" * 40,
    ])

    # Sort by QED score
    if 'QED' in results_df.columns:
        sorted_df = results_df.sort_values('QED', ascending=False)
    else:
        sorted_df = results_df

    for _, row in sorted_df.iterrows():
        if row.get('status') != 'success':
            continue

        lines.extend([
            f"\n{row['compound']}",
            f"  Drug-likeness: {row['drug_likeness']}",
            f"  QED Score: {row['QED']:.3f}",
            f"  MW: {row['MW']:.1f}, LogP: {row['LogP']:.2f}",
            f"  Absorption: {row['absorption_prediction']}",
            f"  BBB: {row['bbb_prediction']}",
        ])

    # Recommendations
    lines.extend([
        "",
        "=" * 70,
        "RECOMMENDATIONS",
        "=" * 70,
        "",
        "Top candidates by drug-likeness (QED score):",
    ])

    if 'QED' in results_df.columns:
        top_qed = sorted_df.head(5)
        for i, (_, row) in enumerate(top_qed.iterrows(), 1):
            if row.get('status') == 'success':
                lines.append(f"  {i}. {row['compound']} (QED={row['QED']:.3f})")

    return '\n'.join(lines)


def main():
    """Main ADMET analysis pipeline."""
    print("=" * 60)
    print("ADMET PREDICTION PIPELINE")
    print("=" * 60)

    # Analyze all candidates
    results_df = analyze_all_candidates()

    # Create report
    report = create_admet_report(results_df)

    # Save results
    print("\n" + "="*60)
    print("Saving Results")
    print("="*60)

    results_df.to_csv(RESULTS_DIR / "admet_analysis.csv", index=False)

    with open(RESULTS_DIR / "admet_report.txt", 'w') as f:
        f.write(report)

    print(f"\nResults saved to {RESULTS_DIR}")
    print("\n" + report)

    return results_df


if __name__ == "__main__":
    results = main()
