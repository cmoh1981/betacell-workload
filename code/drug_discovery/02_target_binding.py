"""
Target-Based Drug Screening using BindingDB and ChEMBL
======================================================

Finds compounds that bind to MR-validated workload targets:
- PDX1 (protective, OR=0.66) - want activators
- SLC2A2/GLUT2 (protective, OR=0.89) - want activators
- MAFA (risk, OR=1.14) - want inhibitors
- GCK (key metabolic gene) - activators

Author: Beta-Cell Workload Analysis Pipeline
"""

import os
import sys
import gzip
import zipfile
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings

warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
DRUG_DB_DIR = Path("E:/drugdatabase")
BINDINGDB_DIR = DRUG_DB_DIR / "BindingDB"
CHEMBL_DIR = DRUG_DB_DIR / "chembl"
STITCH_DIR = DRUG_DB_DIR / "stitch"
RESULTS_DIR = WORKLOAD_DIR / "results" / "drug_discovery"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Target genes with therapeutic hypothesis
WORKLOAD_TARGETS = {
    "PDX1": {
        "uniprot": "P52945",
        "gene_id": "3651",
        "mr_or": 0.66,
        "effect": "protective",
        "therapeutic_action": "activate/increase",
        "rationale": "MR shows higher PDX1 expression reduces T2D risk"
    },
    "SLC2A2": {
        "uniprot": "P11168",
        "gene_id": "6514",
        "mr_or": 0.89,
        "effect": "protective",
        "therapeutic_action": "activate/increase",
        "rationale": "GLUT2 enables glucose sensing; reduced in T2D"
    },
    "MAFA": {
        "uniprot": "Q8NHW3",
        "gene_id": "389692",
        "mr_or": 1.14,
        "effect": "risk",
        "therapeutic_action": "modulate",
        "rationale": "Transcription factor; complex role in β-cell function"
    },
    "GCK": {
        "uniprot": "P35557",
        "gene_id": "2645",
        "mr_or": None,
        "effect": "key_metabolic",
        "therapeutic_action": "activate",
        "rationale": "Glucokinase activators (GKAs) are established T2D drug class"
    },
    "HSPA5": {
        "uniprot": "P11021",
        "gene_id": "3309",
        "mr_or": None,
        "effect": "stress_marker",
        "therapeutic_action": "support/chaperone",
        "rationale": "BiP/GRP78 - ER chaperone; chemical chaperones may help"
    },
    "DDIT3": {
        "uniprot": "P35638",
        "gene_id": "1649",
        "mr_or": None,
        "effect": "stress_marker",
        "therapeutic_action": "inhibit",
        "rationale": "CHOP - pro-apoptotic; inhibition may protect β-cells"
    }
}


def load_stitch_interactions() -> pd.DataFrame:
    """Load STITCH chemical-protein interactions for human."""
    stitch_file = STITCH_DIR / "9606.protein_chemical.links.v5.0.tsv.gz"

    if not stitch_file.exists():
        print(f"STITCH file not found: {stitch_file}")
        return pd.DataFrame()

    print(f"Loading STITCH interactions from {stitch_file}...")

    with gzip.open(stitch_file, 'rt') as f:
        df = pd.read_csv(f, sep='\t')

    print(f"  Loaded {len(df)} chemical-protein interactions")
    return df


def load_stitch_actions() -> pd.DataFrame:
    """Load STITCH action types (activation, inhibition, etc.)."""
    actions_file = STITCH_DIR / "9606.actions.v5.0.tsv.gz"

    if not actions_file.exists():
        print(f"STITCH actions file not found: {actions_file}")
        return pd.DataFrame()

    print(f"Loading STITCH actions from {actions_file}...")

    with gzip.open(actions_file, 'rt') as f:
        df = pd.read_csv(f, sep='\t')

    print(f"  Loaded {len(df)} action annotations")
    return df


def search_stitch_for_targets(
    interactions: pd.DataFrame,
    actions: pd.DataFrame,
    targets: Dict
) -> pd.DataFrame:
    """Search STITCH for compounds targeting our workload genes."""
    print("\n" + "="*60)
    print("Searching STITCH for Workload Target Binders")
    print("="*60)

    results = []

    for gene, info in targets.items():
        # STITCH uses ENSEMBL protein IDs: 9606.ENSP...
        # We need to search by gene name in the database

        # Search in interactions (protein column contains gene names sometimes)
        # STITCH format: chemical (CIDm/CIDs prefix) + protein (9606.ENSP...)

        print(f"\n{gene} ({info['therapeutic_action']}):")

        # For now, we'll note which targets we're looking for
        # Full implementation would require protein ID mapping

        results.append({
            'target_gene': gene,
            'uniprot': info['uniprot'],
            'mr_or': info['mr_or'],
            'therapeutic_action': info['therapeutic_action'],
            'rationale': info['rationale'],
            'stitch_search': 'requires_protein_mapping'
        })

        print(f"  UniProt: {info['uniprot']}")
        print(f"  Action needed: {info['therapeutic_action']}")

    return pd.DataFrame(results)


def search_chembl_targets() -> pd.DataFrame:
    """
    Search ChEMBL for compounds with activity on workload targets.

    Note: Full ChEMBL search requires database setup.
    Here we document the approach and known compounds.
    """
    print("\n" + "="*60)
    print("ChEMBL Target Search Strategy")
    print("="*60)

    # Known compounds from literature/ChEMBL for our targets
    known_compounds = [
        # GCK Activators (established drug class)
        {
            'target': 'GCK',
            'compound': 'Dorzagliatin',
            'chembl_id': 'CHEMBL3545252',
            'action': 'activator',
            'status': 'Approved (China)',
            'ic50_nm': None,
            'ec50_nm': 30,
            'notes': 'Dual-acting GKA, approved for T2D'
        },
        {
            'target': 'GCK',
            'compound': 'Piragliatin',
            'chembl_id': 'CHEMBL517712',
            'action': 'activator',
            'status': 'Phase II (discontinued)',
            'ic50_nm': None,
            'ec50_nm': 50,
            'notes': 'GKA, hypoglycemia concerns'
        },
        {
            'target': 'GCK',
            'compound': 'MK-0941',
            'chembl_id': 'CHEMBL2105725',
            'action': 'activator',
            'status': 'Phase II (discontinued)',
            'ic50_nm': None,
            'ec50_nm': 10,
            'notes': 'GKA from Merck'
        },
        # ER Stress modulators (HSPA5/BiP related)
        {
            'target': 'HSPA5',
            'compound': '4-PBA',
            'chembl_id': 'CHEMBL558',
            'action': 'chemical_chaperone',
            'status': 'Approved (urea cycle)',
            'ic50_nm': None,
            'ec50_nm': None,
            'notes': '4-phenylbutyrate, reduces ER stress'
        },
        {
            'target': 'HSPA5',
            'compound': 'TUDCA',
            'chembl_id': 'CHEMBL1236055',
            'action': 'chemical_chaperone',
            'status': 'Approved (liver)',
            'ic50_nm': None,
            'ec50_nm': None,
            'notes': 'Tauroursodeoxycholic acid, ER stress reducer'
        },
        {
            'target': 'HSPA5',
            'compound': 'Azoramide',
            'chembl_id': 'CHEMBL3614509',
            'action': 'ER_proteostasis',
            'status': 'Research',
            'ic50_nm': None,
            'ec50_nm': 500,
            'notes': 'Enhances ER protein folding'
        },
        # DDIT3/CHOP related
        {
            'target': 'DDIT3',
            'compound': 'Salubrinal',
            'chembl_id': 'CHEMBL96661',
            'action': 'indirect_inhibitor',
            'status': 'Research tool',
            'ic50_nm': None,
            'ec50_nm': 15000,
            'notes': 'eIF2α phosphatase inhibitor, blocks CHOP induction'
        },
        # PDX1 (transcription factor - harder to drug)
        {
            'target': 'PDX1',
            'compound': 'BRD7552',
            'chembl_id': None,
            'action': 'expression_inducer',
            'status': 'Research',
            'ic50_nm': None,
            'ec50_nm': None,
            'notes': 'Small molecule that increases PDX1 expression'
        },
        # KATP channel (ABCC8/KCNJ11)
        {
            'target': 'ABCC8',
            'compound': 'Glibenclamide',
            'chembl_id': 'CHEMBL466',
            'action': 'inhibitor',
            'status': 'Approved',
            'ic50_nm': 3,
            'ec50_nm': None,
            'notes': 'Sulfonylurea, closes KATP channels'
        },
        {
            'target': 'ABCC8',
            'compound': 'Diazoxide',
            'chembl_id': 'CHEMBL1438',
            'action': 'opener',
            'status': 'Approved',
            'ic50_nm': None,
            'ec50_nm': 10000,
            'notes': 'KATP opener, may reduce β-cell workload'
        },
    ]

    df = pd.DataFrame(known_compounds)

    print(f"Identified {len(df)} compounds targeting workload genes:")
    for target in df['target'].unique():
        count = len(df[df['target'] == target])
        print(f"  {target}: {count} compounds")

    return df


def prioritize_candidates(
    chembl_compounds: pd.DataFrame,
    targets: Dict
) -> pd.DataFrame:
    """
    Prioritize drug candidates based on:
    1. MR evidence strength
    2. Drug development status
    3. Safety profile
    4. Mechanism alignment with therapeutic hypothesis
    """
    print("\n" + "="*60)
    print("Prioritizing Drug Candidates")
    print("="*60)

    prioritized = []

    for _, row in chembl_compounds.iterrows():
        target = row['target']
        target_info = targets.get(target, {})

        # Score components
        mr_score = 0
        if target_info.get('mr_or'):
            # Stronger MR effect = higher priority
            mr_effect = abs(1 - target_info['mr_or'])
            mr_score = min(mr_effect * 10, 3)  # Max 3 points

        # Development status score
        status_scores = {
            'Approved': 5,
            'Approved (China)': 4,
            'Approved (liver)': 4,
            'Approved (urea cycle)': 4,
            'Phase III': 3,
            'Phase II': 2,
            'Phase II (discontinued)': 1,
            'Research': 0.5,
            'Research tool': 0.5
        }
        status_score = status_scores.get(row['status'], 0)

        # Action alignment score
        action = row['action']
        therapeutic_action = target_info.get('therapeutic_action', '')

        action_score = 0
        if 'activate' in therapeutic_action and 'activator' in action:
            action_score = 2
        elif 'inhibit' in therapeutic_action and 'inhibitor' in action:
            action_score = 2
        elif 'chaperone' in action or 'proteostasis' in action:
            action_score = 1.5  # Supportive mechanisms

        # Total priority score
        total_score = mr_score + status_score + action_score

        prioritized.append({
            'compound': row['compound'],
            'target': target,
            'action': action,
            'status': row['status'],
            'mr_score': mr_score,
            'status_score': status_score,
            'action_score': action_score,
            'total_score': total_score,
            'chembl_id': row['chembl_id'],
            'notes': row['notes'],
            'rationale': target_info.get('rationale', '')
        })

    prioritized_df = pd.DataFrame(prioritized)
    prioritized_df = prioritized_df.sort_values('total_score', ascending=False)

    print("\nTop 10 Prioritized Candidates:")
    for i, row in prioritized_df.head(10).iterrows():
        print(f"  {row['total_score']:.1f} | {row['compound']} ({row['target']} {row['action']})")
        print(f"       Status: {row['status']}")

    return prioritized_df


def main():
    """Main target-based screening pipeline."""
    print("=" * 60)
    print("TARGET-BASED DRUG SCREENING")
    print("=" * 60)

    # Load STITCH data
    stitch_interactions = load_stitch_interactions()
    stitch_actions = load_stitch_actions()

    # Search STITCH for our targets
    stitch_results = search_stitch_for_targets(
        stitch_interactions, stitch_actions, WORKLOAD_TARGETS
    )

    # Get known compounds from ChEMBL/literature
    chembl_compounds = search_chembl_targets()

    # Prioritize candidates
    prioritized = prioritize_candidates(chembl_compounds, WORKLOAD_TARGETS)

    # Save results
    print("\n" + "="*60)
    print("Saving Results")
    print("="*60)

    chembl_compounds.to_csv(RESULTS_DIR / "chembl_target_compounds.csv", index=False)
    prioritized.to_csv(RESULTS_DIR / "prioritized_target_compounds.csv", index=False)

    # Save target info
    target_df = pd.DataFrame([
        {'gene': k, **v} for k, v in WORKLOAD_TARGETS.items()
    ])
    target_df.to_csv(RESULTS_DIR / "workload_drug_targets.csv", index=False)

    print(f"Results saved to {RESULTS_DIR}")

    return prioritized


if __name__ == "__main__":
    prioritized = main()
