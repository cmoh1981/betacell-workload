"""
Download Public Multi-Omics Data for Beta-Cell Analysis
========================================================

Downloads:
1. Human Protein Atlas - Beta-cell proteomics
2. HMDB - Metabolomics reference data
3. GTEx - eQTL data for pancreas
4. DIAGRAM - T2D GWAS summary statistics

Author: Beta-Cell Workload Analysis Pipeline
"""

import os
import sys
import requests
import gzip
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional
import warnings

warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = WORKLOAD_DIR / "data" / "public_databases"
RESULTS_DIR = WORKLOAD_DIR / "results" / "mr_analysis"

DATA_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Our workload target genes
WORKLOAD_GENES = [
    # MR-validated
    "PDX1", "SLC2A2", "MAFA",
    # Capacity genes
    "GCK", "INS", "NKX6-1", "IAPP", "ABCC8", "KCNJ11",
    # Stress genes
    "DDIT3", "ATF4", "XBP1", "HSPA5", "ERN1", "EIF2AK3",
    # Dedifferentiation
    "ALDH1A3", "SOX9"
]


def download_file(url: str, output_path: Path, description: str = "") -> bool:
    """Download a file with progress indication."""
    try:
        print(f"  Downloading {description}...")
        response = requests.get(url, stream=True, timeout=300)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))

        with open(output_path, 'wb') as f:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    pct = downloaded / total_size * 100
                    print(f"\r    Progress: {pct:.1f}%", end="", flush=True)

        print(f"\r    Downloaded: {output_path.name}")
        return True

    except Exception as e:
        print(f"  ERROR downloading {url}: {e}")
        return False


def download_human_protein_atlas():
    """
    Download Human Protein Atlas pancreas-specific proteomics data.

    Source: https://www.proteinatlas.org/about/download
    """
    print("\n" + "="*60)
    print("DOWNLOADING HUMAN PROTEIN ATLAS DATA")
    print("="*60)

    hpa_dir = DATA_DIR / "human_protein_atlas"
    hpa_dir.mkdir(exist_ok=True)

    # HPA tissue-specific RNA expression
    urls = {
        "rna_tissue_consensus": "https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip",
        "subcellular_location": "https://www.proteinatlas.org/download/subcellular_location.tsv.zip",
        "normal_tissue": "https://www.proteinatlas.org/download/normal_tissue.tsv.zip",
    }

    for name, url in urls.items():
        output_file = hpa_dir / f"{name}.tsv.zip"
        if not output_file.exists():
            download_file(url, output_file, name)
        else:
            print(f"  Already exists: {name}")

    # Extract and filter for pancreas/beta-cell genes
    print("\n  Processing HPA data for pancreas genes...")

    return hpa_dir


def get_gtex_eqtl_data():
    """
    Get GTEx eQTL data for pancreas tissue.

    GTEx V8: https://gtexportal.org/home/downloads/adult-gtex
    Note: Large files, will provide instructions if not available locally.
    """
    print("\n" + "="*60)
    print("GTEx eQTL DATA FOR PANCREAS")
    print("="*60)

    gtex_dir = DATA_DIR / "gtex"
    gtex_dir.mkdir(exist_ok=True)

    # GTEx pancreas eQTL file (V8)
    # Note: GTEx requires agreement to data use policy
    eqtl_file = gtex_dir / "Pancreas.v8.signif_variant_gene_pairs.txt.gz"

    if not eqtl_file.exists():
        print("""
  GTEx pancreas eQTL data not found locally.

  To download GTEx V8 pancreas eQTLs:
  1. Go to: https://gtexportal.org/home/downloads/adult-gtex
  2. Sign data use agreement
  3. Download: GTEx_Analysis_v8_eQTL/Pancreas.v8.signif_variant_gene_pairs.txt.gz
  4. Place in: {gtex_dir}

  Alternative: Use eQTL Catalogue API (we'll use this as backup)
        """.format(gtex_dir=gtex_dir))

        # Use eQTL Catalogue as backup
        return fetch_eqtl_catalogue_data(gtex_dir)

    return eqtl_file


def fetch_eqtl_catalogue_data(output_dir: Path) -> pd.DataFrame:
    """
    Fetch eQTL data from eQTL Catalogue API for our target genes.

    eQTL Catalogue: https://www.ebi.ac.uk/eqtl/api-docs/
    """
    print("\n  Fetching eQTL data from eQTL Catalogue API...")

    # eQTL Catalogue REST API
    base_url = "https://www.ebi.ac.uk/eqtl/api/v2"

    # Pancreas datasets in eQTL Catalogue
    pancreas_datasets = [
        "QTD000584",  # GTEx V8 Pancreas
    ]

    all_eqtls = []

    for gene in WORKLOAD_GENES[:5]:  # Start with key genes
        print(f"    Fetching eQTLs for {gene}...")

        try:
            # Query by gene symbol
            url = f"{base_url}/genes/{gene}/associations"
            params = {
                "size": 1000,
                "qtl_group": "Pancreas"
            }

            response = requests.get(url, params=params, timeout=60)

            if response.status_code == 200:
                data = response.json()
                if '_embedded' in data and 'associations' in data['_embedded']:
                    associations = data['_embedded']['associations']
                    for assoc in associations:
                        all_eqtls.append({
                            'gene': gene,
                            'variant': assoc.get('variant_id', ''),
                            'rsid': assoc.get('rsid', ''),
                            'beta': assoc.get('beta', 0),
                            'se': assoc.get('se', 0),
                            'pvalue': assoc.get('pvalue', 1),
                            'tissue': assoc.get('qtl_group', ''),
                        })
                    print(f"      Found {len(associations)} eQTLs")
            else:
                print(f"      No data found (status: {response.status_code})")

        except Exception as e:
            print(f"      Error: {e}")

    if all_eqtls:
        eqtl_df = pd.DataFrame(all_eqtls)
        output_file = output_dir / "pancreas_eqtls_workload_genes.csv"
        eqtl_df.to_csv(output_file, index=False)
        print(f"\n  Saved {len(eqtl_df)} eQTLs to {output_file.name}")
        return eqtl_df

    return pd.DataFrame()


def get_diagram_gwas():
    """
    Get DIAGRAM consortium T2D GWAS summary statistics.

    Source: https://diagram-consortium.org/downloads.html
    Using: Mahajan et al. 2018 (largest European T2D GWAS)
    """
    print("\n" + "="*60)
    print("DIAGRAM T2D GWAS SUMMARY STATISTICS")
    print("="*60)

    diagram_dir = DATA_DIR / "diagram"
    diagram_dir.mkdir(exist_ok=True)

    # DIAGRAM T2D GWAS (Mahajan 2018)
    # Available at: http://diagram-consortium.org/downloads.html
    gwas_file = diagram_dir / "Mahajan.NatGenet2018b.T2D.European.txt.gz"

    if not gwas_file.exists():
        print("""
  DIAGRAM T2D GWAS data not found locally.

  To download:
  1. Go to: https://diagram-consortium.org/downloads.html
  2. Download: Mahajan et al. 2018 European T2D GWAS
     File: Mahajan.NatGenet2018b.T2D.European.txt.gz
  3. Place in: {diagram_dir}

  Alternative: We'll use OpenGWAS API to fetch relevant SNPs
        """.format(diagram_dir=diagram_dir))

        # Use OpenGWAS as backup
        return fetch_opengwas_t2d(diagram_dir)

    print(f"  Found DIAGRAM GWAS: {gwas_file.name}")
    return gwas_file


def fetch_opengwas_t2d(output_dir: Path) -> pd.DataFrame:
    """
    Fetch T2D GWAS data from OpenGWAS (IEU) for our target gene regions.

    OpenGWAS: https://gwas.mrcieu.ac.uk/
    """
    print("\n  Fetching T2D GWAS data from OpenGWAS...")

    # OpenGWAS API
    base_url = "https://gwas-api.mrcieu.ac.uk"

    # T2D GWAS IDs in OpenGWAS (DIAGRAM-based)
    t2d_gwas_ids = [
        "ebi-a-GCST006867",  # DIAGRAM T2D (Mahajan 2018)
        "ieu-b-24",          # T2D (UKB + DIAGRAM)
    ]

    # Gene regions for our targets (hg38)
    gene_regions = {
        "PDX1": {"chr": "13", "start": 27920000, "end": 27950000},
        "SLC2A2": {"chr": "3", "start": 171016000, "end": 171080000},
        "MAFA": {"chr": "8", "start": 143290000, "end": 143320000},
        "GCK": {"chr": "7", "start": 44183000, "end": 44240000},
        "INS": {"chr": "11", "start": 2158000, "end": 2162000},
    }

    all_associations = []

    for gene, region in gene_regions.items():
        print(f"    Fetching T2D associations near {gene}...")

        try:
            # Query by chromosome region
            url = f"{base_url}/associations"
            params = {
                "id": t2d_gwas_ids[0],
                "chr": region["chr"],
                "start": region["start"],
                "end": region["end"],
            }

            response = requests.post(
                f"{base_url}/associations",
                json={"id": [t2d_gwas_ids[0]],
                      "variant": [],
                      "proxies": 0},
                timeout=60
            )

            # Try alternative approach - tophits
            url = f"{base_url}/tophits"
            response = requests.post(
                url,
                json={"id": [t2d_gwas_ids[0]], "clump": 0},
                timeout=60
            )

            if response.status_code == 200:
                data = response.json()
                for item in data:
                    all_associations.append({
                        'gene_region': gene,
                        'rsid': item.get('rsid', ''),
                        'chr': item.get('chr', ''),
                        'position': item.get('position', 0),
                        'beta': item.get('beta', 0),
                        'se': item.get('se', 0),
                        'pvalue': item.get('p', 1),
                        'effect_allele': item.get('ea', ''),
                        'other_allele': item.get('nea', ''),
                    })
                print(f"      Found associations from GWAS")
                break  # Only need one successful fetch

        except Exception as e:
            print(f"      Error: {e}")

    if all_associations:
        gwas_df = pd.DataFrame(all_associations)
        output_file = output_dir / "t2d_gwas_workload_regions.csv"
        gwas_df.to_csv(output_file, index=False)
        print(f"\n  Saved {len(gwas_df)} associations to {output_file.name}")
        return gwas_df

    return pd.DataFrame()


def download_hmdb_metabolomics():
    """
    Download HMDB metabolomics reference data.

    Source: https://hmdb.ca/downloads
    """
    print("\n" + "="*60)
    print("HMDB METABOLOMICS REFERENCE DATA")
    print("="*60)

    hmdb_dir = DATA_DIR / "hmdb"
    hmdb_dir.mkdir(exist_ok=True)

    # HMDB metabolites related to beta-cell function
    beta_cell_metabolites = [
        # Glucose metabolism
        "Glucose", "Glucose-6-phosphate", "Fructose-6-phosphate",
        "Pyruvate", "Lactate", "ATP", "ADP", "AMP",
        # Amino acids
        "Glutamate", "Glutamine", "Alanine", "Leucine",
        # Lipids
        "Palmitate", "Oleate", "Ceramide",
        # Insulin-related
        "Proinsulin", "C-peptide",
        # TCA cycle
        "Citrate", "Succinate", "Malate",
    ]

    # Create reference file
    metabolite_df = pd.DataFrame({
        'metabolite': beta_cell_metabolites,
        'category': [
            'glucose_metabolism'] * 8 + ['amino_acid'] * 4 +
            ['lipid'] * 3 + ['insulin'] * 2 + ['tca_cycle'] * 3
    })

    output_file = hmdb_dir / "beta_cell_metabolites_reference.csv"
    metabolite_df.to_csv(output_file, index=False)
    print(f"  Created metabolite reference: {output_file.name}")

    return hmdb_dir


def create_data_summary():
    """Create summary of downloaded data."""
    print("\n" + "="*60)
    print("DATA DOWNLOAD SUMMARY")
    print("="*60)

    summary = {
        'component': [],
        'status': [],
        'location': [],
        'notes': []
    }

    # Check each data source
    data_sources = [
        ("Human Protein Atlas", DATA_DIR / "human_protein_atlas", "Pancreas proteomics"),
        ("GTEx eQTLs", DATA_DIR / "gtex", "Pancreas eQTL data"),
        ("DIAGRAM GWAS", DATA_DIR / "diagram", "T2D GWAS summary stats"),
        ("HMDB Reference", DATA_DIR / "hmdb", "Metabolomics reference"),
    ]

    for name, path, notes in data_sources:
        exists = path.exists() and any(path.iterdir()) if path.exists() else False
        summary['component'].append(name)
        summary['status'].append('Available' if exists else 'Missing')
        summary['location'].append(str(path))
        summary['notes'].append(notes)

    summary_df = pd.DataFrame(summary)
    print("\n", summary_df.to_string(index=False))

    # Save summary
    summary_df.to_csv(RESULTS_DIR / "data_availability_summary.csv", index=False)

    return summary_df


def main():
    """Main data download pipeline."""
    print("="*60)
    print("PUBLIC DATABASE DOWNLOAD FOR MR ANALYSIS")
    print("="*60)
    print(f"\nTarget genes: {len(WORKLOAD_GENES)}")
    print(f"Data directory: {DATA_DIR}")

    # Download each data source
    download_human_protein_atlas()
    get_gtex_eqtl_data()
    get_diagram_gwas()
    download_hmdb_metabolomics()

    # Create summary
    summary = create_data_summary()

    print("\n" + "="*60)
    print("DOWNLOAD COMPLETE")
    print("="*60)
    print(f"\nNext steps:")
    print("  1. Run: python 02_mendelian_randomization.py")
    print("  2. Run: python 03_colocalization_analysis.py")

    return summary


if __name__ == "__main__":
    summary = main()
