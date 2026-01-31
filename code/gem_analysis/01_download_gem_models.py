"""
Download Genome-Scale Metabolic Models for Beta-Cell Analysis
=============================================================

Downloads:
1. Recon3D - Human metabolic reconstruction
2. Pancreas/Islet-specific GEM if available
3. iHsa model - Alternative human GEM

Author: Beta-Cell Workload Analysis Pipeline
"""

import os
import sys
import requests
import gzip
import shutil
from pathlib import Path
from typing import Optional
import warnings

warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = WORKLOAD_DIR / "data" / "gem_models"
RESULTS_DIR = WORKLOAD_DIR / "results" / "gem_analysis"

DATA_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# Model sources
GEM_MODELS = {
    "recon3d": {
        "name": "Recon3D",
        "description": "Human metabolic reconstruction (3D)",
        "url": "https://www.vmh.life/files/reconstructions/Recon/3D.01/Recon3D.xml.gz",
        "filename": "Recon3D.xml.gz",
        "reference": "Brunk et al. 2018, Nature Biotechnology"
    },
    "human1": {
        "name": "Human-GEM",
        "description": "Human genome-scale metabolic model",
        "url": "https://github.com/SysBioChalmers/Human-GEM/raw/main/model/Human-GEM.xml",
        "filename": "Human-GEM.xml",
        "reference": "Robinson et al. 2020, Science Signaling"
    },
}

# Beta-cell specific pathways to focus on
BETACELL_PATHWAYS = {
    "glucose_metabolism": [
        "Glycolysis/Gluconeogenesis",
        "Pentose phosphate pathway",
        "Pyruvate metabolism"
    ],
    "insulin_synthesis": [
        "Protein synthesis",
        "ER processing",
        "Vesicle transport"
    ],
    "mitochondrial": [
        "TCA cycle",
        "Oxidative phosphorylation",
        "Fatty acid oxidation"
    ],
    "lipid_metabolism": [
        "Fatty acid synthesis",
        "Sphingolipid metabolism",
        "Cholesterol metabolism"
    ],
    "amino_acid": [
        "Glutamate metabolism",
        "Alanine metabolism",
        "Leucine metabolism"
    ]
}

# Key beta-cell metabolic genes
BETACELL_METABOLIC_GENES = {
    "glucose_sensing": ["GCK", "SLC2A2", "G6PC2", "PFKFB2"],
    "pyruvate_metabolism": ["PDK1", "PDK2", "PDK3", "PDK4", "PC", "LDHA", "LDHB"],
    "tca_cycle": ["CS", "ACO2", "IDH2", "OGDH", "SUCLA2", "SUCLG2", "FH", "MDH2"],
    "oxphos": ["NDUFS1", "SDHA", "UQCRC1", "COX4I1", "ATP5F1A", "ATP5F1B"],
    "fatty_acid": ["CPT1A", "CPT2", "ACADVL", "ACADM", "HADHA", "HADHB"],
    "insulin_processing": ["PCSK1", "PCSK2", "CPE", "CHGB", "SLC30A8"],
    "upr_metabolic": ["XBP1", "ATF6", "EIF2AK3", "ATF4", "DDIT3"]
}


def download_file(url: str, output_path: Path, description: str = "") -> bool:
    """Download a file with progress indication."""
    try:
        print(f"  Downloading {description}...")
        response = requests.get(url, stream=True, timeout=600, allow_redirects=True)
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


def download_recon3d():
    """Download Recon3D human metabolic model."""
    print("\n" + "=" * 60)
    print("DOWNLOADING RECON3D MODEL")
    print("=" * 60)

    model_info = GEM_MODELS["recon3d"]
    output_gz = DATA_DIR / model_info["filename"]
    output_xml = DATA_DIR / "Recon3D.xml"

    if output_xml.exists():
        print(f"  Already exists: {output_xml.name}")
        return output_xml

    if not output_gz.exists():
        success = download_file(
            model_info["url"],
            output_gz,
            model_info["name"]
        )
        if not success:
            return None

    # Decompress
    print("  Decompressing...")
    try:
        with gzip.open(output_gz, 'rb') as f_in:
            with open(output_xml, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"  Decompressed to: {output_xml.name}")
        return output_xml
    except Exception as e:
        print(f"  ERROR decompressing: {e}")
        return None


def download_human_gem():
    """Download Human-GEM model."""
    print("\n" + "=" * 60)
    print("DOWNLOADING HUMAN-GEM MODEL")
    print("=" * 60)

    model_info = GEM_MODELS["human1"]
    output_path = DATA_DIR / model_info["filename"]

    if output_path.exists():
        print(f"  Already exists: {output_path.name}")
        return output_path

    success = download_file(
        model_info["url"],
        output_path,
        model_info["name"]
    )

    return output_path if success else None


def create_betacell_gene_mapping():
    """Create mapping of beta-cell genes to metabolic pathways."""
    print("\n" + "=" * 60)
    print("CREATING BETA-CELL GENE MAPPING")
    print("=" * 60)

    import pandas as pd

    # Flatten gene mapping
    rows = []
    for pathway, genes in BETACELL_METABOLIC_GENES.items():
        for gene in genes:
            rows.append({
                'gene_symbol': gene,
                'pathway': pathway,
                'category': 'metabolic'
            })

    gene_df = pd.DataFrame(rows)
    output_path = DATA_DIR / "betacell_metabolic_genes.csv"
    gene_df.to_csv(output_path, index=False)
    print(f"  Saved: {output_path.name} ({len(gene_df)} genes)")

    return gene_df


def verify_model_structure(model_path: Path) -> dict:
    """Verify model can be loaded and report basic statistics."""
    print(f"\n  Verifying model: {model_path.name}")

    try:
        from cobra.io import read_sbml_model

        model = read_sbml_model(str(model_path))

        stats = {
            'model_id': model.id,
            'n_reactions': len(model.reactions),
            'n_metabolites': len(model.metabolites),
            'n_genes': len(model.genes),
            'n_compartments': len(model.compartments)
        }

        print(f"    Model ID: {stats['model_id']}")
        print(f"    Reactions: {stats['n_reactions']}")
        print(f"    Metabolites: {stats['n_metabolites']}")
        print(f"    Genes: {stats['n_genes']}")
        print(f"    Compartments: {stats['n_compartments']}")

        return stats

    except ImportError:
        print("    COBRApy not installed. Install with: pip install cobra")
        return None
    except Exception as e:
        print(f"    ERROR loading model: {e}")
        return None


def create_download_summary():
    """Create summary of downloaded models."""
    print("\n" + "=" * 60)
    print("DOWNLOAD SUMMARY")
    print("=" * 60)

    import pandas as pd

    summary = {
        'model': [],
        'status': [],
        'file': [],
        'reference': []
    }

    for model_id, info in GEM_MODELS.items():
        # Check for both .xml and .xml.gz files
        xml_path = DATA_DIR / info["filename"].replace(".gz", "")
        gz_path = DATA_DIR / info["filename"]

        if xml_path.exists():
            status = "Available"
            filepath = str(xml_path)
        elif gz_path.exists():
            status = "Compressed"
            filepath = str(gz_path)
        else:
            status = "Not downloaded"
            filepath = ""

        summary['model'].append(info['name'])
        summary['status'].append(status)
        summary['file'].append(filepath)
        summary['reference'].append(info['reference'])

    df = pd.DataFrame(summary)
    print(df.to_string(index=False))

    # Save summary
    df.to_csv(RESULTS_DIR / "gem_models_summary.csv", index=False)

    return df


def main():
    """Main download pipeline."""
    print("=" * 60)
    print("GEM MODEL DOWNLOAD FOR BETA-CELL ANALYSIS")
    print("=" * 60)
    print(f"\nData directory: {DATA_DIR}")

    # Download models
    recon3d_path = download_recon3d()
    human_gem_path = download_human_gem()

    # Create gene mapping
    gene_mapping = create_betacell_gene_mapping()

    # Verify models if downloaded
    if recon3d_path and recon3d_path.exists():
        verify_model_structure(recon3d_path)

    if human_gem_path and human_gem_path.exists():
        verify_model_structure(human_gem_path)

    # Create summary
    summary = create_download_summary()

    print("\n" + "=" * 60)
    print("DOWNLOAD COMPLETE")
    print("=" * 60)
    print("\nNext steps:")
    print("  1. Run: python 02_betacell_gem_analysis.py")
    print("  2. Run: python 03_workload_metabolic_integration.py")

    return summary


if __name__ == "__main__":
    summary = main()
