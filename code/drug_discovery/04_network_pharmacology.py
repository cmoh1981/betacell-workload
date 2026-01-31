"""
Network Pharmacology Analysis for Beta-Cell Workload
=====================================================

Integrates drug-target-pathway networks to identify:
1. Multi-target drug effects on workload pathways
2. Synergistic drug combinations
3. Network-based drug prioritization

Uses: STRING (PPI), STITCH (drug-protein), MSigDB (pathways)

Author: Beta-Cell Workload Analysis Pipeline
"""

import os
import sys
import gzip
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
import warnings

warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
DRUG_DB_DIR = Path("E:/drugdatabase")
STRING_DIR = DRUG_DB_DIR / "string"
STITCH_DIR = DRUG_DB_DIR / "stitch"
MSIGDB_DIR = DRUG_DB_DIR / "msigdb_data"
RESULTS_DIR = WORKLOAD_DIR / "results" / "drug_discovery"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Try to import networkx
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    print("NetworkX not installed. Install with: pip install networkx")

# Workload-related genes from our analysis
WORKLOAD_GENES = {
    # MR-validated targets
    "PDX1": {"category": "capacity", "mr_or": 0.66, "effect": "protective"},
    "SLC2A2": {"category": "capacity", "mr_or": 0.89, "effect": "protective"},
    "MAFA": {"category": "capacity", "mr_or": 1.14, "effect": "risk"},

    # Metabolic/functional genes
    "GCK": {"category": "capacity", "mr_or": None, "effect": "key_metabolic"},
    "INS": {"category": "capacity", "mr_or": None, "effect": "biosynthetic"},
    "IAPP": {"category": "capacity", "mr_or": None, "effect": "biosynthetic"},
    "NKX6-1": {"category": "capacity", "mr_or": None, "effect": "identity"},
    "ABCC8": {"category": "capacity", "mr_or": None, "effect": "secretion"},
    "KCNJ11": {"category": "capacity", "mr_or": None, "effect": "secretion"},

    # Stress markers
    "DDIT3": {"category": "stress", "mr_or": None, "effect": "apoptotic"},
    "ATF4": {"category": "stress", "mr_or": None, "effect": "UPR"},
    "XBP1": {"category": "stress", "mr_or": None, "effect": "UPR"},
    "HSPA5": {"category": "stress", "mr_or": None, "effect": "chaperone"},
    "ERN1": {"category": "stress", "mr_or": None, "effect": "UPR"},
    "EIF2AK3": {"category": "stress", "mr_or": None, "effect": "UPR"},

    # Dedifferentiation markers
    "ALDH1A3": {"category": "dediff", "mr_or": None, "effect": "dediff_marker"},
    "SOX9": {"category": "dediff", "mr_or": None, "effect": "progenitor"},
}

# Drug candidates from our target-based screening
DRUG_CANDIDATES = {
    "Dorzagliatin": {"targets": ["GCK"], "action": "activator", "status": "approved"},
    "Piragliatin": {"targets": ["GCK"], "action": "activator", "status": "phase2"},
    "4-PBA": {"targets": ["HSPA5"], "action": "chaperone", "status": "approved"},
    "TUDCA": {"targets": ["HSPA5"], "action": "chaperone", "status": "approved"},
    "Salubrinal": {"targets": ["EIF2AK3", "DDIT3"], "action": "modulator", "status": "research"},
    "Glibenclamide": {"targets": ["ABCC8"], "action": "inhibitor", "status": "approved"},
    "Diazoxide": {"targets": ["ABCC8"], "action": "opener", "status": "approved"},
    "BRD7552": {"targets": ["PDX1"], "action": "inducer", "status": "research"},
    "Metformin": {"targets": ["PRKAA1", "PRKAA2"], "action": "activator", "status": "approved"},
    "Dapagliflozin": {"targets": ["SLC5A2"], "action": "inhibitor", "status": "approved"},
    "Empagliflozin": {"targets": ["SLC5A2"], "action": "inhibitor", "status": "approved"},
}

# Beta-cell relevant pathways
RELEVANT_PATHWAYS = [
    "KEGG_TYPE_II_DIABETES_MELLITUS",
    "KEGG_INSULIN_SIGNALING_PATHWAY",
    "KEGG_MATURITY_ONSET_DIABETES_OF_THE_YOUNG",
    "REACTOME_REGULATION_OF_INSULIN_SECRETION",
    "REACTOME_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_PANCREAS_BETA_CELLS",
    "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
    "GOBP_INSULIN_SECRETION",
    "GOBP_GLUCOSE_HOMEOSTASIS",
]


def load_string_aliases() -> Dict[str, str]:
    """Load STRING protein aliases to map gene symbols to STRING IDs."""
    alias_file = STRING_DIR / "9606.protein.aliases.v12.0.txt.gz"

    if not alias_file.exists():
        print(f"STRING aliases file not found: {alias_file}")
        return {}

    print(f"Loading STRING aliases from {alias_file}...")

    gene_to_string = {}
    with gzip.open(alias_file, 'rt', encoding='utf-8') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                string_id = parts[0]  # 9606.ENSP...
                alias = parts[1]
                source = parts[2]

                # Prioritize gene symbols
                if source in ['BioMart_HUGO', 'Ensembl_HGNC', 'BLAST_KEGG_NAME']:
                    gene_to_string[alias] = string_id

    print(f"  Loaded {len(gene_to_string)} gene-to-STRING mappings")
    return gene_to_string


def load_string_interactions(gene_to_string: Dict, genes_of_interest: Set[str],
                            min_score: int = 400) -> pd.DataFrame:
    """Load STRING protein-protein interactions for genes of interest."""
    links_file = STRING_DIR / "9606.protein.links.v12.0.txt.gz"

    if not links_file.exists():
        print(f"STRING links file not found: {links_file}")
        return pd.DataFrame()

    print(f"Loading STRING interactions from {links_file}...")

    # Get STRING IDs for our genes
    string_ids = set()
    for gene in genes_of_interest:
        if gene in gene_to_string:
            string_ids.add(gene_to_string[gene])

    print(f"  Found STRING IDs for {len(string_ids)}/{len(genes_of_interest)} genes")

    # Load interactions involving our genes
    interactions = []
    with gzip.open(links_file, 'rt', encoding='utf-8') as f:
        header = next(f).strip().split()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                protein1, protein2, score = parts[0], parts[1], int(parts[2])

                if score >= min_score:
                    if protein1 in string_ids or protein2 in string_ids:
                        interactions.append({
                            'protein1': protein1,
                            'protein2': protein2,
                            'combined_score': score
                        })

    df = pd.DataFrame(interactions)
    print(f"  Loaded {len(df)} interactions (score >= {min_score})")
    return df


def load_msigdb_pathways(pathway_file: str = "c2.all.v7.5.1.symbols.gmt") -> Dict[str, Set[str]]:
    """Load MSigDB pathway gene sets."""
    gmt_file = MSIGDB_DIR / pathway_file

    if not gmt_file.exists():
        print(f"MSigDB file not found: {gmt_file}")
        return {}

    print(f"Loading MSigDB pathways from {gmt_file}...")

    pathways = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pathway_name = parts[0]
                genes = set(parts[2:])
                pathways[pathway_name] = genes

    print(f"  Loaded {len(pathways)} pathways")
    return pathways


def load_hallmark_pathways() -> Dict[str, Set[str]]:
    """Load MSigDB Hallmark gene sets."""
    return load_msigdb_pathways("h.all.v7.5.1.symbols.gmt")


def build_workload_network(interactions: pd.DataFrame,
                          gene_to_string: Dict,
                          pathways: Dict[str, Set[str]]) -> Optional[object]:
    """Build a network of workload genes and their interactions."""
    if not HAS_NETWORKX:
        print("NetworkX not available, skipping network construction")
        return None

    print("\n" + "="*60)
    print("Building Workload Gene Network")
    print("="*60)

    # Create reverse mapping (STRING ID to gene symbol)
    string_to_gene = {v: k for k, v in gene_to_string.items()}

    # Create network
    G = nx.Graph()

    # Add workload genes as nodes
    for gene, info in WORKLOAD_GENES.items():
        G.add_node(gene, **info, node_type='workload_gene')

    # Add edges from STRING interactions
    edges_added = 0
    for _, row in interactions.iterrows():
        gene1 = string_to_gene.get(row['protein1'])
        gene2 = string_to_gene.get(row['protein2'])

        if gene1 and gene2:
            G.add_edge(gene1, gene2, weight=row['combined_score']/1000)
            edges_added += 1

    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Annotate nodes with pathway membership
    for gene in G.nodes():
        gene_pathways = []
        for pathway_name, pathway_genes in pathways.items():
            if gene in pathway_genes:
                gene_pathways.append(pathway_name)
        G.nodes[gene]['pathways'] = gene_pathways
        G.nodes[gene]['n_pathways'] = len(gene_pathways)

    return G


def analyze_drug_network_effects(drug_candidates: Dict,
                                 workload_genes: Dict,
                                 pathways: Dict[str, Set[str]],
                                 G: Optional[object] = None) -> pd.DataFrame:
    """Analyze how drugs affect the workload gene network."""
    print("\n" + "="*60)
    print("Analyzing Drug Network Effects")
    print("="*60)

    results = []

    for drug_name, drug_info in drug_candidates.items():
        targets = drug_info['targets']
        action = drug_info['action']

        # Count workload genes affected
        direct_workload_targets = [t for t in targets if t in workload_genes]

        # Find pathways affected
        affected_pathways = set()
        for target in targets:
            for pathway_name, pathway_genes in pathways.items():
                if target in pathway_genes:
                    affected_pathways.add(pathway_name)

        # Count relevant pathways (beta-cell related)
        relevant_affected = [p for p in affected_pathways
                           if any(kw in p.upper() for kw in
                                 ['DIABETES', 'INSULIN', 'PANCREA', 'GLUCOSE',
                                  'UNFOLDED', 'ENDOPLASMIC', 'BETA'])]

        # Network analysis if available
        network_centrality = 0
        neighbors_count = 0
        if G is not None and HAS_NETWORKX:
            for target in targets:
                if target in G.nodes():
                    # Get neighbors (potential indirect effects)
                    neighbors = set(G.neighbors(target))
                    neighbors_count += len(neighbors)

                    # Calculate centrality
                    try:
                        centrality = nx.degree_centrality(G)
                        network_centrality += centrality.get(target, 0)
                    except:
                        pass

        # Calculate network pharmacology score
        # Higher is better for workload modulation
        np_score = (
            len(direct_workload_targets) * 3 +  # Direct workload targets
            len(relevant_affected) * 2 +         # Relevant pathways
            min(neighbors_count / 10, 2) +       # Network reach
            network_centrality * 5               # Network importance
        )

        results.append({
            'drug': drug_name,
            'targets': ', '.join(targets),
            'action': action,
            'status': drug_info['status'],
            'direct_workload_targets': len(direct_workload_targets),
            'total_pathways_affected': len(affected_pathways),
            'relevant_pathways': len(relevant_affected),
            'network_neighbors': neighbors_count,
            'network_centrality': round(network_centrality, 4),
            'np_score': round(np_score, 2),
            'relevant_pathway_names': '; '.join(relevant_affected[:5])
        })

        print(f"\n{drug_name}:")
        print(f"  Targets: {', '.join(targets)}")
        print(f"  Direct workload targets: {len(direct_workload_targets)}")
        print(f"  Relevant pathways: {len(relevant_affected)}")
        print(f"  NP Score: {np_score:.2f}")

    df = pd.DataFrame(results)
    df = df.sort_values('np_score', ascending=False)

    return df


def identify_synergistic_combinations(drug_effects: pd.DataFrame,
                                     workload_genes: Dict,
                                     pathways: Dict) -> pd.DataFrame:
    """Identify potentially synergistic drug combinations."""
    print("\n" + "="*60)
    print("Identifying Synergistic Drug Combinations")
    print("="*60)

    combinations = []
    drugs = drug_effects['drug'].tolist()

    for i, drug1 in enumerate(drugs):
        for drug2 in drugs[i+1:]:
            drug1_info = DRUG_CANDIDATES[drug1]
            drug2_info = DRUG_CANDIDATES[drug2]

            targets1 = set(drug1_info['targets'])
            targets2 = set(drug2_info['targets'])

            # Check for complementary mechanisms
            # Different targets = potentially synergistic
            target_overlap = len(targets1 & targets2)
            unique_targets = len(targets1 | targets2)

            # Check category coverage
            categories1 = set()
            categories2 = set()

            for t in targets1:
                if t in workload_genes:
                    categories1.add(workload_genes[t]['category'])
            for t in targets2:
                if t in workload_genes:
                    categories2.add(workload_genes[t]['category'])

            combined_categories = categories1 | categories2

            # Pathway coverage
            pathways1 = set()
            pathways2 = set()
            for t in targets1:
                for pname, pgenes in pathways.items():
                    if t in pgenes:
                        pathways1.add(pname)
            for t in targets2:
                for pname, pgenes in pathways.items():
                    if t in pgenes:
                        pathways2.add(pname)

            unique_pathways = len(pathways1 | pathways2)
            pathway_overlap = len(pathways1 & pathways2)

            # Synergy score
            # Reward: different targets, multiple categories, unique pathways
            # Penalize: same targets, same category only
            synergy_score = (
                (unique_targets - target_overlap) * 2 +  # Unique targets
                len(combined_categories) * 3 +            # Category coverage
                unique_pathways * 0.1 -                   # Pathway breadth
                pathway_overlap * 0.05                    # Slight penalty for overlap
            )

            # Mechanism complementarity
            mechanisms = []
            if 'capacity' in combined_categories and 'stress' in combined_categories:
                mechanisms.append("Capacity + Stress reduction")
                synergy_score += 3
            if drug1_info['action'] != drug2_info['action']:
                mechanisms.append(f"{drug1_info['action']} + {drug2_info['action']}")
                synergy_score += 1

            # Status bonus (approved + approved = clinically feasible)
            if drug1_info['status'] == 'approved' and drug2_info['status'] == 'approved':
                synergy_score += 2
                mechanisms.append("Both approved")

            combinations.append({
                'drug1': drug1,
                'drug2': drug2,
                'targets1': ', '.join(targets1),
                'targets2': ', '.join(targets2),
                'unique_targets': unique_targets,
                'categories_covered': ', '.join(combined_categories),
                'n_categories': len(combined_categories),
                'unique_pathways': unique_pathways,
                'synergy_score': round(synergy_score, 2),
                'mechanisms': '; '.join(mechanisms) if mechanisms else 'Standard combination'
            })

    df = pd.DataFrame(combinations)
    df = df.sort_values('synergy_score', ascending=False)

    print(f"\nTop 10 Synergistic Combinations:")
    for _, row in df.head(10).iterrows():
        print(f"  {row['synergy_score']:.1f} | {row['drug1']} + {row['drug2']}")
        print(f"       Categories: {row['categories_covered']}")
        print(f"       {row['mechanisms']}")

    return df


def pathway_enrichment_analysis(genes: List[str],
                               pathways: Dict[str, Set[str]],
                               background_size: int = 20000) -> pd.DataFrame:
    """Perform pathway enrichment analysis for a gene set."""
    from scipy import stats

    results = []
    gene_set = set(genes)

    for pathway_name, pathway_genes in pathways.items():
        # Overlap
        overlap = gene_set & pathway_genes

        if len(overlap) == 0:
            continue

        # Fisher's exact test
        # a = overlap, b = pathway - overlap
        # c = query - overlap, d = background - all
        a = len(overlap)
        b = len(pathway_genes) - a
        c = len(gene_set) - a
        d = background_size - a - b - c

        _, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative='greater')

        # Enrichment ratio
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
        # Multiple testing correction (Benjamini-Hochberg)
        df = df.sort_values('pvalue')
        df['rank'] = range(1, len(df) + 1)
        df['fdr'] = df['pvalue'] * len(df) / df['rank']
        df['fdr'] = df['fdr'].clip(upper=1.0)
        df = df.sort_values('fdr')

    return df


def create_network_report(drug_effects: pd.DataFrame,
                         synergies: pd.DataFrame,
                         workload_enrichment: pd.DataFrame) -> str:
    """Create summary network pharmacology report."""
    lines = [
        "=" * 70,
        "NETWORK PHARMACOLOGY REPORT - WORKLOAD MODULATORS",
        "=" * 70,
        "",
        "SUMMARY",
        "-" * 40,
        f"Drugs analyzed: {len(drug_effects)}",
        f"Drug combinations evaluated: {len(synergies)}",
        f"Enriched pathways: {len(workload_enrichment[workload_enrichment['fdr'] < 0.05]) if len(workload_enrichment) > 0 else 0}",
        "",
        "TOP DRUGS BY NETWORK PHARMACOLOGY SCORE",
        "-" * 40,
    ]

    for _, row in drug_effects.head(5).iterrows():
        lines.extend([
            f"\n{row['drug']} (NP Score: {row['np_score']})",
            f"  Targets: {row['targets']}",
            f"  Action: {row['action']}",
            f"  Status: {row['status']}",
            f"  Direct workload targets: {row['direct_workload_targets']}",
            f"  Relevant pathways: {row['relevant_pathways']}",
        ])

    lines.extend([
        "",
        "=" * 70,
        "TOP SYNERGISTIC COMBINATIONS",
        "=" * 70,
    ])

    for _, row in synergies.head(5).iterrows():
        lines.extend([
            f"\n{row['drug1']} + {row['drug2']} (Synergy: {row['synergy_score']})",
            f"  Categories: {row['categories_covered']}",
            f"  Mechanism: {row['mechanisms']}",
            f"  Unique targets: {row['unique_targets']}",
        ])

    if len(workload_enrichment) > 0:
        sig_pathways = workload_enrichment[workload_enrichment['fdr'] < 0.05]
        if len(sig_pathways) > 0:
            lines.extend([
                "",
                "=" * 70,
                "ENRICHED PATHWAYS (FDR < 0.05)",
                "=" * 70,
            ])

            for _, row in sig_pathways.head(10).iterrows():
                lines.append(f"  {row['pathway']}: {row['overlap_count']} genes (FDR={row['fdr']:.2e})")

    lines.extend([
        "",
        "=" * 70,
        "THERAPEUTIC RECOMMENDATIONS",
        "=" * 70,
        "",
        "Based on network analysis, recommended treatment strategies:",
        "",
        "1. MONOTHERAPY (highest NP scores):",
    ])

    for _, row in drug_effects.head(3).iterrows():
        lines.append(f"   - {row['drug']}: {row['action']} of {row['targets']}")

    lines.extend([
        "",
        "2. COMBINATION THERAPY (highest synergy):",
    ])

    for _, row in synergies.head(3).iterrows():
        lines.append(f"   - {row['drug1']} + {row['drug2']}: {row['mechanisms']}")

    lines.extend([
        "",
        "3. MULTI-TARGET STRATEGY:",
        "   Combine drugs targeting different workload components:",
        "   - Capacity enhancers (GCK activators)",
        "   - Stress reducers (ER chaperones)",
        "   - Demand reducers (SGLT2 inhibitors)",
    ])

    return '\n'.join(lines)


def main():
    """Main network pharmacology analysis pipeline."""
    print("=" * 60)
    print("NETWORK PHARMACOLOGY ANALYSIS")
    print("=" * 60)

    # Load STRING data
    gene_to_string = load_string_aliases()

    # Get all genes of interest
    all_genes = set(WORKLOAD_GENES.keys())
    for drug_info in DRUG_CANDIDATES.values():
        all_genes.update(drug_info['targets'])

    # Load STRING interactions
    interactions = load_string_interactions(gene_to_string, all_genes)

    # Load pathways
    pathways = load_msigdb_pathways("c2.all.v7.5.1.symbols.gmt")
    hallmark = load_hallmark_pathways()
    pathways.update(hallmark)

    # Also load GO biological process
    go_bp = load_msigdb_pathways("c5.all.v7.5.1.symbols.gmt")
    pathways.update(go_bp)

    print(f"\nTotal pathways loaded: {len(pathways)}")

    # Build network
    G = build_workload_network(interactions, gene_to_string, pathways)

    # Analyze drug network effects
    drug_effects = analyze_drug_network_effects(
        DRUG_CANDIDATES, WORKLOAD_GENES, pathways, G
    )

    # Identify synergistic combinations
    synergies = identify_synergistic_combinations(
        drug_effects, WORKLOAD_GENES, pathways
    )

    # Pathway enrichment for workload genes
    print("\n" + "="*60)
    print("Workload Gene Pathway Enrichment")
    print("="*60)

    workload_enrichment = pathway_enrichment_analysis(
        list(WORKLOAD_GENES.keys()), pathways
    )

    sig_pathways = workload_enrichment[workload_enrichment['fdr'] < 0.05]
    print(f"\nSignificant pathways (FDR < 0.05): {len(sig_pathways)}")
    for _, row in sig_pathways.head(10).iterrows():
        print(f"  {row['pathway']}: {row['overlap_count']} genes (FDR={row['fdr']:.2e})")

    # Create report
    report = create_network_report(drug_effects, synergies, workload_enrichment)

    # Save results
    print("\n" + "="*60)
    print("Saving Results")
    print("="*60)

    drug_effects.to_csv(RESULTS_DIR / "network_drug_effects.csv", index=False)
    synergies.to_csv(RESULTS_DIR / "synergistic_combinations.csv", index=False)
    workload_enrichment.to_csv(RESULTS_DIR / "workload_pathway_enrichment.csv", index=False)

    with open(RESULTS_DIR / "network_pharmacology_report.txt", 'w') as f:
        f.write(report)

    print(f"\nResults saved to {RESULTS_DIR}")
    print("\n" + report)

    return drug_effects, synergies


if __name__ == "__main__":
    drug_effects, synergies = main()
