# Beta-Cell T2D Workload Project

## Project Overview

This project focuses on advanced single-cell RNA-seq analysis and multi-omics data integration for identifying therapeutic targets in Type 2 Diabetes (T2D) beta-cell dysfunction.

## Project Focus Areas

1. **Single-Cell Analysis** - Deep analysis of scRNA-seq data from human pancreatic islets
2. **Literature Review** - Systematic review of beta-cell dysfunction mechanisms
3. **Data Integration** - Multi-dataset integration and cross-validation

## Directory Structure

```
workload/
├── data/
│   ├── raw/           # Raw input data files
│   └── processed/     # Processed AnnData objects
├── code/
│   ├── singlecell/    # Single-cell analysis scripts
│   ├── integration/   # Data integration pipelines
│   └── utils/         # Utility functions
├── results/
│   ├── figures/       # Generated figures
│   └── tables/        # Output tables and statistics
├── manuscript/        # Manuscript drafts
├── references/        # Literature and citations
└── docs/              # Documentation
```

## Data Sources

| Dataset | GEO/ArrayExpress | Cells | Description |
|---------|------------------|-------|-------------|
| Segerstolpe | E-MTAB-5061 | ~3,500 | T2D vs Normal islets |
| GSE221156 | GSE221156 | ~245,000 | Large islet atlas |
| Enge | GSE81547 | ~2,500 | Aging pancreas |

## Key Analysis Targets

- Beta-cell identity markers (INS, MAFA, PDX1, NKX6-1, UCN3)
- ER stress/UPR pathway genes (XBP1, ATF6, HSPA5)
- Dedifferentiation markers (ALDH1A3, GASTRIN)
- Metabolic regulators (GCK, PC, TXNIP)

## Getting Started

```bash
# Install dependencies
pip install -r requirements.txt

# Run initial analysis
python code/singlecell/01_load_data.py
```

## Requirements

- Python 3.9+
- scanpy >= 1.9.0
- anndata >= 0.9.0
- pandas >= 2.0.0
- numpy >= 1.24.0

---
*Project initialized: 2026-01-30*
