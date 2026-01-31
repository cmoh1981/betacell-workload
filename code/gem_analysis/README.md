# GEM Analysis Module for Beta-Cell Workload

## Overview

This module integrates Genome-Scale Metabolic Modeling (GEM) with the Composite Workload Index (CWI) analysis to understand metabolic changes across beta-cell workload states.

## Pipeline Components

### 1. `01_download_gem_models.py`
Downloads genome-scale metabolic models:
- **Human-GEM**: 12,971 reactions, 8,455 metabolites, 2,887 genes
- Creates beta-cell specific gene mappings

### 2. `02_betacell_gem_analysis.py`
Performs constraint-based metabolic analysis:
- **Flux Balance Analysis (FBA)**: Predict optimal metabolic fluxes
- **Parsimonious FBA (pFBA)**: Minimize total flux at optimality
- **Flux Variability Analysis (FVA)**: Determine flux ranges
- **Gene Knockout Analysis**: Identify essential genes

### 3. `03_workload_metabolic_integration.py`
Integrates CWI with GEM:
- Maps workload states to metabolic phenotypes
- Identifies metabolic bottlenecks in stressed cells
- Predicts therapeutic targets

## Workload State → Metabolic Phenotype Mapping

| Workload State | Glucose Uptake | Insulin Demand | ATP Requirement | UPR Activity |
|----------------|----------------|----------------|-----------------|--------------|
| S1_Resting     | Low (5)        | Low (0.2)      | Normal (1.0)    | None (0.0)   |
| S2_Active      | High (15)      | High (0.8)     | Elevated (1.5)  | Low (0.1)    |
| S3_Stressed    | Very High (18) | Max (1.0)      | High (1.8)      | Moderate (0.5)|
| S4_Exhausted   | Reduced (12)   | Medium (0.6)   | Elevated (1.2)  | High (0.8)   |
| S5_Failing     | Low (8)        | Low (0.3)      | Reduced (0.8)   | Max (1.0)    |

## Key Pathways Analyzed

1. **Glucose Sensing**: GCK, SLC2A2, G6PC2
2. **Glycolysis**: PFK, ALDOA, PKM
3. **TCA Cycle**: CS, IDH2, OGDH, MDH2, PC
4. **Oxidative Phosphorylation**: NADH dehydrogenase, ATP synthase
5. **Insulin Synthesis**: PCSK1, PCSK2, CPE
6. **UPR/Stress Response**: XBP1, ATF4, DDIT3

## Results

Results are saved to `results/gem_analysis/`:
- `workload_metabolic_integration.csv`: CWI-GEM mapping
- `metabolic_bottlenecks.csv`: Constrained reactions in exhausted state
- `therapeutic_targets.csv`: Predicted drug targets
- `gene_knockout_results.csv`: Essential gene analysis
- `fva_top_variable.csv`: Most variable reactions

## Usage

```bash
# 1. Download GEM models
python 01_download_gem_models.py

# 2. Run GEM analysis
python 02_betacell_gem_analysis.py

# 3. Integrate with CWI workload scores
python 03_workload_metabolic_integration.py
```

## Dependencies

- COBRApy >= 0.30.0
- pandas
- numpy
- requests

Install with: `pip install cobra pandas numpy requests`

## Integration with Workload Analysis

This module connects with:
- **CWI Workload Scores** from `results/multiomics/workload_scores.csv`
- **DIABLO Multi-omics Signatures** from multiomics analysis
- **MR-validated Targets** (PDX1, SLC2A2, MAFA) for genetic validation

## References

- Human-GEM: Robinson et al. 2020, Science Signaling
- COBRApy: Ebrahim et al. 2013, BMC Systems Biology
- Beta-cell metabolism: Prentki & Nolan 2006, JCI
