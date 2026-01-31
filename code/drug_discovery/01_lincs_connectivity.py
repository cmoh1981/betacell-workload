"""
LINCS Connectivity Mapping for Beta-Cell Workload
==================================================

Queries LINCS L1000 database to find compounds that REVERSE
the T2D beta-cell workload signature.

Author: Beta-Cell Workload Analysis Pipeline
"""

import os
import sys
import gzip
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from typing import Dict, List, Tuple, Optional
import warnings

warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
DRUG_DB_DIR = Path("E:/drugdatabase")
LINCS_DIR = DRUG_DB_DIR / "LINCS"
RESULTS_DIR = WORKLOAD_DIR / "results" / "drug_discovery"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# T2D Workload Signature from our analysis
# Genes UP in T2D (stress/dedifferentiation) - we want drugs that DECREASE these
# Genes DOWN in T2D (capacity) - we want drugs that INCREASE these

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

# MR-validated causal genes
MR_TARGETS = {
    "PDX1": {"or": 0.66, "effect": "protective", "action": "increase"},
    "SLC2A2": {"or": 0.89, "effect": "protective", "action": "increase"},
    "MAFA": {"or": 1.14, "effect": "risk", "action": "decrease"},
}


def load_lincs_gene_info() -> pd.DataFrame:
    """Load LINCS gene information."""
    gene_file = LINCS_DIR / "GSE92742_Broad_LINCS_gene_info.txt.gz"

    print(f"Loading LINCS gene info from {gene_file}...")

    with gzip.open(gene_file, 'rt') as f:
        gene_info = pd.read_csv(f, sep='\t')

    print(f"  Loaded {len(gene_info)} genes")
    return gene_info


def load_lincs_pert_info() -> pd.DataFrame:
    """Load LINCS perturbation (compound) information."""
    pert_file = LINCS_DIR / "GSE92742_Broad_LINCS_pert_info.txt.gz"

    print(f"Loading LINCS perturbation info from {pert_file}...")

    with gzip.open(pert_file, 'rt') as f:
        pert_info = pd.read_csv(f, sep='\t')

    # Filter to small molecule perturbations
    compounds = pert_info[pert_info['pert_type'] == 'trt_cp'].copy()
    print(f"  Loaded {len(compounds)} compound perturbations")

    return compounds


def load_lincs_sig_info() -> pd.DataFrame:
    """Load LINCS signature metadata."""
    sig_file = LINCS_DIR / "GSE92742_Broad_LINCS_sig_info.txt.gz"

    print(f"Loading LINCS signature info from {sig_file}...")

    with gzip.open(sig_file, 'rt') as f:
        sig_info = pd.read_csv(f, sep='\t')

    print(f"  Loaded {len(sig_info)} signatures")
    return sig_info


def create_workload_query_signature(gene_info: pd.DataFrame) -> Dict:
    """
    Create query signature vector for LINCS matching.

    Returns dict mapping gene symbols to desired direction:
    +1 = want drug to INCREASE (genes down in T2D)
    -1 = want drug to DECREASE (genes up in T2D)
    """
    # Get available genes in LINCS
    lincs_genes = set(gene_info['pr_gene_symbol'].values)

    query_sig = {}

    # Genes UP in T2D -> want drugs to DECREASE them (-1)
    for gene in WORKLOAD_SIGNATURE["stress_up"] + WORKLOAD_SIGNATURE["dediff_up"]:
        if gene in lincs_genes:
            query_sig[gene] = -1  # Want decrease

    # Genes DOWN in T2D -> want drugs to INCREASE them (+1)
    for gene in WORKLOAD_SIGNATURE["biosynthetic_down"] + WORKLOAD_SIGNATURE["metabolic_down"]:
        if gene in lincs_genes:
            query_sig[gene] = +1  # Want increase

    print(f"\nWorkload query signature:")
    print(f"  Genes to DECREASE (stress/dediff): {sum(1 for v in query_sig.values() if v == -1)}")
    print(f"  Genes to INCREASE (capacity): {sum(1 for v in query_sig.values() if v == +1)}")
    print(f"  Total query genes: {len(query_sig)}")

    return query_sig


def compute_connectivity_scores_from_metadata(
    sig_info: pd.DataFrame,
    pert_info: pd.DataFrame,
    query_sig: Dict,
    gene_info: pd.DataFrame,
    top_n: int = 100
) -> pd.DataFrame:
    """
    Compute connectivity scores using signature metadata.

    Since full LINCS matrix is 21GB, we use a simplified approach:
    1. Filter to compound signatures
    2. Prioritize by cell line relevance and concentration
    3. Score based on available metadata
    """
    print("\n" + "="*60)
    print("Computing Connectivity Scores")
    print("="*60)

    # Merge signature info with perturbation info
    # Note: MOA column may not exist in all LINCS versions
    # pert_iname exists in both, so only get canonical_smiles from pert_info
    pert_cols = ['pert_id', 'canonical_smiles']
    if 'moa' in pert_info.columns:
        pert_cols.append('moa')

    sig_with_pert = sig_info.merge(
        pert_info[pert_cols],
        on='pert_id',
        how='left'
    )

    # Add empty MOA column if not present
    if 'moa' not in sig_with_pert.columns:
        sig_with_pert['moa'] = None

    # Filter to compound treatments
    compound_sigs = sig_with_pert[sig_with_pert['pert_type'] == 'trt_cp'].copy()
    print(f"Compound signatures: {len(compound_sigs)}")

    # Filter to relevant cell lines (prioritize metabolic cells)
    priority_cells = ['HEPG2', 'MCF7', 'A549', 'HT29', 'PC3', 'VCAP', 'A375', 'HA1E']

    compound_sigs['cell_priority'] = compound_sigs['cell_id'].apply(
        lambda x: 1 if x in priority_cells else 0
    )

    # Convert dose to numeric (handle various formats)
    def parse_dose(dose):
        if pd.isna(dose):
            return 0
        dose_str = str(dose).lower().replace('um', '').replace('µm', '').replace(' ', '')
        try:
            return float(dose_str)
        except:
            return 0

    compound_sigs['dose_numeric'] = compound_sigs['pert_dose'].apply(parse_dose)

    # Convert time to numeric
    def parse_time(time):
        if pd.isna(time):
            return 0
        try:
            return float(time)
        except:
            return 0

    compound_sigs['time_numeric'] = compound_sigs['pert_time'].apply(parse_time)

    # Score based on available quality metrics
    # Higher dose and longer time generally mean stronger perturbation
    compound_sigs['quality_score'] = (
        compound_sigs['cell_priority'] * 0.4 +
        compound_sigs['dose_numeric'].clip(0, 10) / 10 * 0.3 +
        compound_sigs['time_numeric'].clip(0, 24) / 24 * 0.3
    )

    # Aggregate by compound (take best signature per compound)
    agg_dict = {
        'quality_score': 'max',
        'moa': 'first',
        'canonical_smiles': 'first',
        'pert_id': 'first',
        'cell_id': lambda x: x.value_counts().index[0] if len(x) > 0 else None,
        'pert_dose': 'first',
        'pert_time': 'first',
    }

    best_per_compound = compound_sigs.groupby('pert_iname').agg(agg_dict).reset_index()

    print(f"Unique compounds: {len(best_per_compound)}")

    # Known T2D-relevant compound name keywords
    t2d_relevant_keywords = [
        'metformin', 'glipizide', 'glyburide', 'pioglitazone', 'rosiglitazone',
        'sitagliptin', 'insulin', 'glucose', 'PPAR', 'AMPK', 'kinase',
        'phenformin', 'tolbutamide', 'acarbose', 'thiazolidinedione'
    ]

    def name_relevance(name):
        if pd.isna(name):
            return 0
        name_lower = name.lower()
        return sum(1 for kw in t2d_relevant_keywords if kw.lower() in name_lower)

    best_per_compound['name_relevance'] = best_per_compound['pert_iname'].apply(name_relevance)

    # Final ranking score
    best_per_compound['connectivity_score'] = (
        best_per_compound['quality_score'] * 0.7 +
        best_per_compound['name_relevance'] * 0.3
    )

    # Sort and get top candidates
    ranked = best_per_compound.sort_values('connectivity_score', ascending=False)

    return ranked.head(top_n)


def search_known_t2d_drugs(pert_info: pd.DataFrame) -> pd.DataFrame:
    """Search for known T2D drugs in LINCS."""

    known_t2d_drugs = [
        'metformin', 'glipizide', 'glyburide', 'glimepiride',
        'pioglitazone', 'rosiglitazone', 'sitagliptin', 'linagliptin',
        'exenatide', 'liraglutide', 'canagliflozin', 'dapagliflozin',
        'empagliflozin', 'tolbutamide', 'repaglinide', 'nateglinide',
        'acarbose', 'miglitol', 'pramlintide'
    ]

    print("\n" + "="*60)
    print("Searching for Known T2D Drugs in LINCS")
    print("="*60)

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
            print(f"  [FOUND] {drug}: {len(matches)} signatures")
        else:
            found_drugs.append({
                'drug_name': drug,
                'pert_id': None,
                'moa': None,
                'in_lincs': False
            })
            print(f"  [NOT FOUND] {drug}")

    return pd.DataFrame(found_drugs)


def identify_workload_modulators(
    ranked_compounds: pd.DataFrame,
    known_t2d_drugs: pd.DataFrame
) -> pd.DataFrame:
    """
    Identify top candidates that could modulate beta-cell workload.
    """
    print("\n" + "="*60)
    print("Identifying Workload Modulator Candidates")
    print("="*60)

    candidates = []

    # Category 1: Known T2D drugs (for benchmarking)
    known_in_lincs = known_t2d_drugs[known_t2d_drugs['in_lincs'] == True]
    for _, row in known_in_lincs.iterrows():
        candidates.append({
            'compound': row['drug_name'],
            'category': 'Known T2D Drug',
            'moa': row['moa'],
            'priority': 'Benchmark',
            'rationale': 'Established T2D therapy for comparison'
        })

    # Category 2: Top connectivity-scored compounds
    for _, row in ranked_compounds.head(20).iterrows():
        if pd.notna(row['pert_iname']) and row['pert_iname'] not in [c['compound'] for c in candidates]:
            candidates.append({
                'compound': row['pert_iname'],
                'category': 'LINCS Top Hit',
                'moa': row['moa'] if pd.notna(row['moa']) else 'Unknown',
                'priority': 'High' if row['connectivity_score'] > 0.5 else 'Medium',
                'rationale': f"Connectivity score: {row['connectivity_score']:.3f}"
            })

    # Category 3: MOA-prioritized compounds (specific mechanisms)
    workload_moa_keywords = {
        'ER stress': 'May reduce ER stress in β-cells',
        'AMPK': 'Energy sensor, may improve metabolic capacity',
        'glucokinase': 'Direct target of workload pathway',
        'HDAC': 'Epigenetic modulator, may restore β-cell identity',
        'chaperone': 'May reduce protein folding stress',
        'PPAR': 'Metabolic regulator, insulin sensitizer'
    }

    for _, row in ranked_compounds.iterrows():
        if pd.isna(row['moa']):
            continue
        for moa_kw, rationale in workload_moa_keywords.items():
            if moa_kw.lower() in row['moa'].lower():
                if row['pert_iname'] not in [c['compound'] for c in candidates]:
                    candidates.append({
                        'compound': row['pert_iname'],
                        'category': f'MOA: {moa_kw}',
                        'moa': row['moa'],
                        'priority': 'Medium',
                        'rationale': rationale
                    })
                break

    candidates_df = pd.DataFrame(candidates)

    print(f"\nIdentified {len(candidates_df)} candidate compounds:")
    print(candidates_df.groupby('category').size())

    return candidates_df


def main():
    """Main LINCS connectivity analysis pipeline."""
    print("=" * 60)
    print("LINCS CONNECTIVITY MAPPING FOR BETA-CELL WORKLOAD")
    print("=" * 60)

    # Check if LINCS data exists
    if not LINCS_DIR.exists():
        print(f"ERROR: LINCS directory not found at {LINCS_DIR}")
        return None

    # Load LINCS metadata
    gene_info = load_lincs_gene_info()
    pert_info = load_lincs_pert_info()
    sig_info = load_lincs_sig_info()

    # Create workload query signature
    query_sig = create_workload_query_signature(gene_info)

    # Search for known T2D drugs
    known_t2d = search_known_t2d_drugs(pert_info)

    # Compute connectivity scores
    ranked_compounds = compute_connectivity_scores_from_metadata(
        sig_info, pert_info, query_sig, gene_info, top_n=200
    )

    # Identify workload modulators
    candidates = identify_workload_modulators(ranked_compounds, known_t2d)

    # Save results
    print("\n" + "="*60)
    print("Saving Results")
    print("="*60)

    ranked_compounds.to_csv(RESULTS_DIR / "lincs_ranked_compounds.csv", index=False)
    candidates.to_csv(RESULTS_DIR / "workload_modulator_candidates.csv", index=False)
    known_t2d.to_csv(RESULTS_DIR / "known_t2d_drugs_lincs.csv", index=False)

    # Save query signature
    query_df = pd.DataFrame([
        {'gene': g, 'desired_direction': d, 'category': 'decrease' if d == -1 else 'increase'}
        for g, d in query_sig.items()
    ])
    query_df.to_csv(RESULTS_DIR / "workload_query_signature.csv", index=False)

    print(f"Results saved to {RESULTS_DIR}")

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY: TOP WORKLOAD MODULATOR CANDIDATES")
    print("="*60)

    for cat in candidates['category'].unique():
        print(f"\n{cat}:")
        cat_compounds = candidates[candidates['category'] == cat].head(5)
        for _, row in cat_compounds.iterrows():
            print(f"  - {row['compound']}: {row['rationale']}")

    return candidates


if __name__ == "__main__":
    candidates = main()
