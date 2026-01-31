# Project Status - Beta-Cell T2D Workload

**Created:** 2026-01-30
**Status:** Architecture Complete

## Workload Framework Defined

Based on literature review (PMC5021190, Frontiers 2024), we define beta-cell workload through:

### Three-Pillar Framework
1. **Biosynthetic Demand** - Insulin production, ER folding, secretory machinery
2. **Metabolic Capacity** - Glucose sensing, mitochondrial function, β-cell identity
3. **Stress Response** - UPR activation, oxidative stress, inflammation

### Five Workload States
| State | CWI Range | Description |
|-------|-----------|-------------|
| S1_Resting | < 0.5 | Low demand, healthy |
| S2_Active | 0.5-1.0 | Normal secretory activity |
| S3_Stressed | 1.0-1.5 | Adaptive UPR, compensating |
| S4_Exhausted | 1.5-2.5 | Terminal UPR, decompensating |
| S5_Failing | > 2.5 | Dedifferentiation, identity loss |

## Project Structure

```
workload/
├── README.md                 ✅ Project overview
├── requirements.txt          ✅ Python dependencies
├── config.yaml               ✅ Configuration settings
├── PROJECT_STATUS.md         ✅ This file
├── data/
│   ├── raw/                  ✅ For raw input data
│   └── processed/            ✅ For processed h5ad files
├── code/
│   ├── singlecell/
│   │   ├── 01_load_data.py                ✅ Data loading
│   │   ├── 02_differential_expression.py  ✅ DE analysis
│   │   ├── 03_workload_scoring.py         ✅ CWI computation
│   │   └── 04_workload_visualization.py   ✅ Visualization
│   ├── integration/
│   │   └── 01_integrate_datasets.py       ✅ Multi-dataset integration
│   └── utils/
│       └── helpers.py                     ✅ Utility functions
├── results/
│   ├── figures/              ✅ Output figures
│   └── tables/               ✅ Output tables
├── manuscript/               ✅ Manuscript drafts
├── references/               ✅ Literature citations
└── docs/
    ├── ANALYSIS_PLAN.md                   ✅ Analysis roadmap
    └── BETA_CELL_WORKLOAD_ARCHITECTURE.md ✅ Workload framework
```

## Data Sources (from betacell project)

| Dataset | Status | Path |
|---------|--------|------|
| Segerstolpe | Available | ../betacell/data/processed/segerstolpe_processed.h5ad |
| GSE221156 | Available | ../betacell/data/processed/gse221156_atlas_processed.h5ad |
| Enge | Available | ../betacell/data/processed/enge_aging_processed.h5ad |

## Next Steps

1. Install dependencies: `pip install -r requirements.txt`
2. Run data loading: `python code/singlecell/01_load_data.py`
3. Verify data access and quality
4. Begin differential expression analysis

## Notes

- Project builds on previous betacell analysis (~34,500 cells analyzed)
- Focus: Single-cell analysis, literature review, data integration
- Key targets: XBP1, GCK, PC, TXNIP, PDX1
