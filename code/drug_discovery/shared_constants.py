"""
Shared Constants for Drug Discovery Pipeline
=============================================

Centralizes all constants used across modules to avoid duplication.
"""

from pathlib import Path

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
DRUG_DB_DIR = Path("E:/drugdatabase")
RESULTS_DIR = WORKLOAD_DIR / "results" / "drug_discovery"

# Drug-likeness rules
LIPINSKI_RULES = {
    'MW': {'max': 500, 'description': 'Molecular Weight <= 500'},
    'LogP': {'max': 5, 'description': 'LogP <= 5'},
    'HBD': {'max': 5, 'description': 'H-bond donors <= 5'},
    'HBA': {'max': 10, 'description': 'H-bond acceptors <= 10'}
}

VEBER_RULES = {
    'RotatableBonds': {'max': 10, 'description': 'Rotatable bonds <= 10'},
    'TPSA': {'max': 140, 'description': 'TPSA <= 140 A^2'}
}

# Workload-related genes from MR analysis
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

# Drug candidates from target-based screening
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

# Known compound SMILES
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

# Workload query signature
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
