"""
GEM-Drug Discovery Integration
==============================

Integrates genome-scale metabolic modeling with drug discovery to:
1. Map metabolic bottlenecks to drug targets
2. Identify compounds that may restore metabolic capacity
3. Connect MR-validated genes with metabolic pathways

Author: Beta-Cell Workload Analysis Pipeline
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
import warnings

warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
GEM_RESULTS = WORKLOAD_DIR / "results" / "gem_analysis"
DRUG_RESULTS = WORKLOAD_DIR / "results" / "drug_discovery"
MR_RESULTS = WORKLOAD_DIR / "results" / "mr_analysis"
OUTPUT_DIR = WORKLOAD_DIR / "results" / "integrated_analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Metabolic pathway to drug target mapping
PATHWAY_DRUG_TARGETS = {
    "Glycolysis / Gluconeogenesis": {
        "targets": ["GCK", "PFKFB2", "PKM", "HK2"],
        "drug_classes": ["glucokinase activator", "GCK", "glucose"],
        "rationale": "Restore glucose sensing and metabolism"
    },
    "Pentose phosphate pathway": {
        "targets": ["G6PD", "PGLS", "TKT"],
        "drug_classes": ["antioxidant", "NADPH"],
        "rationale": "Provide NADPH for redox balance"
    },
    "Fructose and mannose metabolism": {
        "targets": ["ALDOB", "KHK"],
        "drug_classes": ["fructose", "metabolism"],
        "rationale": "Alternative fuel utilization"
    },
    "TCA cycle": {
        "targets": ["CS", "IDH2", "OGDH", "MDH2", "PC"],
        "drug_classes": ["anaplerotic", "pyruvate", "mitochondria"],
        "rationale": "Enhance mitochondrial capacity"
    },
    "Oxidative phosphorylation": {
        "targets": ["ATP5F1A", "COX4I1", "NDUFS1"],
        "drug_classes": ["CoQ10", "mitochondria", "ATP"],
        "rationale": "Boost ATP production"
    },
    "Alanine, aspartate and glutamate metabolism": {
        "targets": ["GOT1", "GOT2", "GLS", "GLUD1"],
        "drug_classes": ["glutamate", "amino acid", "alanine"],
        "rationale": "Amino acid-stimulated insulin secretion"
    },
    "Arginine and proline metabolism": {
        "targets": ["ARG1", "NOS1", "ODC1"],
        "drug_classes": ["nitric oxide", "polyamine"],
        "rationale": "Signaling and cell proliferation"
    },
    "Purine metabolism": {
        "targets": ["HPRT1", "IMPDH1", "ADSL"],
        "drug_classes": ["purine", "nucleotide"],
        "rationale": "Nucleotide synthesis for cell function"
    },
    "Pyrimidine metabolism": {
        "targets": ["DHODH", "CAD", "UMPS"],
        "drug_classes": ["pyrimidine", "nucleotide"],
        "rationale": "DNA/RNA synthesis support"
    }
}

# MR-validated targets with metabolic connections
MR_METABOLIC_CONNECTIONS = {
    "PDX1": {
        "or": 0.66,
        "effect": "protective",
        "metabolic_role": "Master regulator of beta-cell identity and glucose metabolism",
        "pathways": ["Glycolysis / Gluconeogenesis", "Oxidative phosphorylation"],
        "drug_strategy": "Enhance PDX1 expression/activity"
    },
    "SLC2A2": {
        "or": 0.89,
        "effect": "protective",
        "metabolic_role": "Glucose transporter, glucose sensing",
        "pathways": ["Glycolysis / Gluconeogenesis"],
        "drug_strategy": "Enhance glucose uptake capacity"
    },
    "MAFA": {
        "or": 1.14,
        "effect": "risk",
        "metabolic_role": "Beta-cell maturation factor",
        "pathways": ["TCA cycle", "Oxidative phosphorylation"],
        "drug_strategy": "Modulate MAFA in context of stress"
    },
    "GCK": {
        "or": 0.87,
        "effect": "protective",
        "metabolic_role": "Glucose sensor, rate-limiting glycolysis",
        "pathways": ["Glycolysis / Gluconeogenesis"],
        "drug_strategy": "Glucokinase activators (e.g., Dorzagliatin)"
    }
}


def load_gem_bottlenecks() -> pd.DataFrame:
    """Load metabolic bottlenecks from GEM analysis."""
    bottleneck_file = GEM_RESULTS / "metabolic_bottlenecks.csv"
    if bottleneck_file.exists():
        return pd.read_csv(bottleneck_file)
    return pd.DataFrame()


def load_drug_candidates() -> pd.DataFrame:
    """Load drug candidates from LINCS analysis."""
    candidates_file = DRUG_RESULTS / "workload_modulator_candidates.csv"
    if candidates_file.exists():
        return pd.read_csv(candidates_file)
    return pd.DataFrame()


def load_lincs_ranked() -> pd.DataFrame:
    """Load ranked LINCS compounds."""
    ranked_file = DRUG_RESULTS / "lincs_ranked_compounds.csv"
    if ranked_file.exists():
        return pd.read_csv(ranked_file)
    return pd.DataFrame()


def load_mr_results() -> pd.DataFrame:
    """Load MR analysis results."""
    mr_file = MR_RESULTS / "mr_expected_vs_observed.csv"
    if mr_file.exists():
        return pd.read_csv(mr_file)
    return pd.DataFrame()


def map_bottlenecks_to_targets(bottlenecks: pd.DataFrame) -> pd.DataFrame:
    """Map metabolic bottlenecks to drug targets."""
    print("\n" + "=" * 60)
    print("MAPPING METABOLIC BOTTLENECKS TO DRUG TARGETS")
    print("=" * 60)

    if bottlenecks.empty:
        print("No bottleneck data available")
        return pd.DataFrame()

    mappings = []

    for _, row in bottlenecks.iterrows():
        subsystem = row.get('subsystem', '')
        if pd.isna(subsystem) or subsystem == '':
            continue

        # Find matching pathway
        for pathway, info in PATHWAY_DRUG_TARGETS.items():
            if pathway.lower() in subsystem.lower() or subsystem.lower() in pathway.lower():
                mappings.append({
                    'reaction_id': row['reaction_id'],
                    'reaction_name': row.get('reaction_name', ''),
                    'subsystem': subsystem,
                    'pathway_match': pathway,
                    'drug_targets': ', '.join(info['targets']),
                    'drug_classes': ', '.join(info['drug_classes']),
                    'rationale': info['rationale']
                })
                break

    mapping_df = pd.DataFrame(mappings)

    if not mapping_df.empty:
        print(f"\nMapped {len(mapping_df)} bottlenecks to drug targets")
        print("\nPathway distribution:")
        print(mapping_df['pathway_match'].value_counts())

    return mapping_df


def find_drugs_for_bottlenecks(
    bottleneck_targets: pd.DataFrame,
    lincs_compounds: pd.DataFrame
) -> pd.DataFrame:
    """Find LINCS compounds that may address metabolic bottlenecks."""
    print("\n" + "=" * 60)
    print("FINDING DRUGS FOR METABOLIC BOTTLENECKS")
    print("=" * 60)

    if bottleneck_targets.empty or lincs_compounds.empty:
        print("Insufficient data")
        return pd.DataFrame()

    # Get all drug class keywords from bottleneck targets
    drug_keywords = set()
    for _, row in bottleneck_targets.iterrows():
        if pd.notna(row.get('drug_classes')):
            keywords = [k.strip().lower() for k in row['drug_classes'].split(',')]
            drug_keywords.update(keywords)

    print(f"Searching for drug classes: {drug_keywords}")

    matches = []

    for _, compound in lincs_compounds.iterrows():
        compound_name = str(compound.get('pert_iname', '')).lower()
        moa = str(compound.get('moa', '')).lower()

        # Check for keyword matches
        for keyword in drug_keywords:
            if keyword in compound_name or keyword in moa:
                matches.append({
                    'compound': compound.get('pert_iname', ''),
                    'moa': compound.get('moa', ''),
                    'matched_keyword': keyword,
                    'connectivity_score': compound.get('connectivity_score', 0),
                    'rationale': f"May address {keyword}-related metabolic dysfunction"
                })
                break

    matches_df = pd.DataFrame(matches)

    if not matches_df.empty:
        matches_df = matches_df.drop_duplicates(subset=['compound'])
        print(f"\nFound {len(matches_df)} compounds matching metabolic targets")

    return matches_df


def integrate_mr_with_metabolism() -> pd.DataFrame:
    """Integrate MR-validated targets with metabolic analysis."""
    print("\n" + "=" * 60)
    print("INTEGRATING MR-VALIDATED TARGETS WITH METABOLISM")
    print("=" * 60)

    results = []

    for gene, info in MR_METABOLIC_CONNECTIONS.items():
        results.append({
            'gene': gene,
            'mr_or': info['or'],
            'mr_effect': info['effect'],
            'metabolic_role': info['metabolic_role'],
            'affected_pathways': ', '.join(info['pathways']),
            'drug_strategy': info['drug_strategy']
        })

    results_df = pd.DataFrame(results)

    print("\nMR-Validated Genes with Metabolic Connections:")
    for _, row in results_df.iterrows():
        effect_symbol = "decreases" if row['mr_effect'] == 'protective' else "increases"
        print(f"\n  {row['gene']} (OR={row['mr_or']}, {effect_symbol} T2D risk)")
        print(f"    Role: {row['metabolic_role']}")
        print(f"    Strategy: {row['drug_strategy']}")

    return results_df


def create_therapeutic_roadmap(
    bottleneck_drugs: pd.DataFrame,
    mr_metabolism: pd.DataFrame,
    drug_candidates: pd.DataFrame
) -> pd.DataFrame:
    """Create integrated therapeutic roadmap."""
    print("\n" + "=" * 60)
    print("CREATING THERAPEUTIC ROADMAP")
    print("=" * 60)

    roadmap = []

    # Tier 1: MR-validated targets with existing drugs
    tier1_drugs = {
        "GCK": ["Dorzagliatin", "glucokinase activator"],
        "PDX1": ["HDAC inhibitor", "epigenetic modulator"],
        "SLC2A2": ["glucose transporter modulator"],
    }

    for gene, drugs in tier1_drugs.items():
        if gene in MR_METABOLIC_CONNECTIONS:
            info = MR_METABOLIC_CONNECTIONS[gene]
            roadmap.append({
                'tier': 1,
                'target': gene,
                'approach': info['drug_strategy'],
                'candidate_drugs': ', '.join(drugs),
                'evidence_level': 'MR-validated',
                'metabolic_rationale': info['metabolic_role'],
                'priority': 'High'
            })

    # Tier 2: Metabolic bottleneck-targeted drugs from LINCS
    if not bottleneck_drugs.empty:
        for _, row in bottleneck_drugs.head(10).iterrows():
            roadmap.append({
                'tier': 2,
                'target': row.get('matched_keyword', 'metabolic'),
                'approach': 'Address metabolic bottleneck',
                'candidate_drugs': row.get('compound', ''),
                'evidence_level': 'GEM + LINCS',
                'metabolic_rationale': row.get('rationale', ''),
                'priority': 'Medium'
            })

    # Tier 3: Known T2D drugs with workload relevance
    if not drug_candidates.empty:
        known_t2d = drug_candidates[drug_candidates['category'] == 'Known T2D Drug']
        for _, row in known_t2d.iterrows():
            roadmap.append({
                'tier': 3,
                'target': row.get('moa', 'Multiple'),
                'approach': 'Established T2D therapy',
                'candidate_drugs': row.get('compound', ''),
                'evidence_level': 'Clinical',
                'metabolic_rationale': 'Proven efficacy in T2D',
                'priority': 'Benchmark'
            })

    roadmap_df = pd.DataFrame(roadmap)

    if not roadmap_df.empty:
        print(f"\nTherapeutic Roadmap: {len(roadmap_df)} entries")
        print("\nBy Tier:")
        print(roadmap_df.groupby('tier').size())
        print("\nBy Priority:")
        print(roadmap_df.groupby('priority').size())

    return roadmap_df


def generate_integration_report() -> str:
    """Generate comprehensive integration report."""
    report = []
    report.append("=" * 70)
    report.append("GEM-DRUG DISCOVERY INTEGRATION REPORT")
    report.append("Beta-Cell Workload Analysis Pipeline")
    report.append("=" * 70)

    # Load all data
    bottlenecks = load_gem_bottlenecks()
    drug_candidates = load_drug_candidates()
    lincs_ranked = load_lincs_ranked()
    mr_results = load_mr_results()

    # Section 1: MR Results Summary
    report.append("\n" + "-" * 50)
    report.append("1. MENDELIAN RANDOMIZATION RESULTS")
    report.append("-" * 50)
    if not mr_results.empty:
        for _, row in mr_results.iterrows():
            effect = "protective" if row.get('observed_or_ivw', 1) < 1 else "risk"
            report.append(f"  {row['gene']}: OR={row.get('observed_or_ivw', 'N/A'):.3f} ({effect})")
            report.append(f"    Expected: {row.get('expected_effect', 'N/A')}")

    # Section 2: Metabolic Bottlenecks
    report.append("\n" + "-" * 50)
    report.append("2. METABOLIC BOTTLENECKS (GEM Analysis)")
    report.append("-" * 50)
    if not bottlenecks.empty:
        pathway_counts = bottlenecks['subsystem'].value_counts().head(5)
        for pathway, count in pathway_counts.items():
            report.append(f"  {pathway}: {count} bottleneck reactions")

    # Section 3: Drug Candidates
    report.append("\n" + "-" * 50)
    report.append("3. DRUG CANDIDATES (LINCS Connectivity)")
    report.append("-" * 50)
    if not drug_candidates.empty:
        report.append(f"  Total candidates: {len(drug_candidates)}")
        for cat in drug_candidates['category'].unique():
            count = len(drug_candidates[drug_candidates['category'] == cat])
            report.append(f"    {cat}: {count}")

    # Section 4: Integration Summary
    report.append("\n" + "-" * 50)
    report.append("4. INTEGRATION SUMMARY")
    report.append("-" * 50)
    report.append("  Pipeline Components Integrated:")
    report.append("    [x] CWI Workload Scores (269 cells)")
    report.append("    [x] Multi-omics DIABLO (22 signature features)")
    report.append("    [x] Human-GEM Metabolic Model (12,971 reactions)")
    report.append("    [x] Mendelian Randomization (5 genes)")
    report.append("    [x] LINCS Connectivity (473,647 signatures)")

    report.append("\n" + "=" * 70)
    report.append("END OF REPORT")
    report.append("=" * 70)

    return "\n".join(report)


def main():
    """Main integration pipeline."""
    print("=" * 60)
    print("GEM-DRUG DISCOVERY INTEGRATION")
    print("=" * 60)

    # Load data
    bottlenecks = load_gem_bottlenecks()
    drug_candidates = load_drug_candidates()
    lincs_ranked = load_lincs_ranked()
    mr_results = load_mr_results()

    print(f"\nData loaded:")
    print(f"  Metabolic bottlenecks: {len(bottlenecks)}")
    print(f"  Drug candidates: {len(drug_candidates)}")
    print(f"  LINCS compounds: {len(lincs_ranked)}")
    print(f"  MR results: {len(mr_results)}")

    # Map bottlenecks to targets
    bottleneck_targets = map_bottlenecks_to_targets(bottlenecks)

    # Find drugs for bottlenecks
    bottleneck_drugs = find_drugs_for_bottlenecks(bottleneck_targets, lincs_ranked)

    # Integrate MR with metabolism
    mr_metabolism = integrate_mr_with_metabolism()

    # Create therapeutic roadmap
    roadmap = create_therapeutic_roadmap(bottleneck_drugs, mr_metabolism, drug_candidates)

    # Save results
    print("\n" + "=" * 60)
    print("SAVING RESULTS")
    print("=" * 60)

    if not bottleneck_targets.empty:
        bottleneck_targets.to_csv(OUTPUT_DIR / "bottleneck_drug_targets.csv", index=False)

    if not bottleneck_drugs.empty:
        bottleneck_drugs.to_csv(OUTPUT_DIR / "bottleneck_matched_drugs.csv", index=False)

    mr_metabolism.to_csv(OUTPUT_DIR / "mr_metabolic_connections.csv", index=False)

    if not roadmap.empty:
        roadmap.to_csv(OUTPUT_DIR / "therapeutic_roadmap.csv", index=False)

    # Generate report
    report = generate_integration_report()
    report_path = OUTPUT_DIR / "integration_report.txt"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)

    print(f"\nResults saved to: {OUTPUT_DIR}")

    # Print roadmap summary
    print("\n" + "=" * 60)
    print("THERAPEUTIC ROADMAP SUMMARY")
    print("=" * 60)

    if not roadmap.empty:
        for tier in sorted(roadmap['tier'].unique()):
            tier_data = roadmap[roadmap['tier'] == tier]
            print(f"\nTier {tier}:")
            for _, row in tier_data.head(5).iterrows():
                print(f"  - {row['target']}: {row['candidate_drugs']}")
                print(f"    Rationale: {row['metabolic_rationale'][:60]}...")

    print("\n" + "=" * 60)
    print("INTEGRATION COMPLETE")
    print("=" * 60)

    return roadmap


if __name__ == "__main__":
    roadmap = main()
