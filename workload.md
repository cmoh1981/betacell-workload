# Beta-Cell Workload Hypothesis: A Systems Biology Approach to Type 2 Diabetes

## Executive Summary

This project presents a comprehensive systems biology framework for understanding Type 2 Diabetes (T2D) progression through the lens of **beta-cell workload**. By integrating multi-omics data, genome-scale metabolic modeling, Mendelian randomization, and drug discovery databases, we provide evidence that chronic metabolic demand exceeding beta-cell capacity drives progressive dysfunction leading to T2D.

---

## Table of Contents

1. [The Workload Hypothesis](#1-the-workload-hypothesis)
2. [Data Sources & Datasets](#2-data-sources--datasets)
3. [Tools & Methods](#3-tools--methods)
4. [Analysis Pipeline](#4-analysis-pipeline)
5. [Key Results](#5-key-results)
6. [Therapeutic Implications](#6-therapeutic-implications)
7. [File Structure](#7-file-structure)
8. [Reproducibility](#8-reproducibility)
9. [References](#9-references)

---

## 1. The Workload Hypothesis

### Core Concept

The **beta-cell workload hypothesis** proposes that T2D develops through a progressive cascade when chronic metabolic **demand exceeds capacity**:

```
NORMAL → COMPENSATION → STRESS → EXHAUSTION → FAILURE
   │           │            │          │           │
   │           │            │          │           └── Clinical T2D
   │           │            │          └── Dedifferentiation
   │           │            └── ER Stress/UPR
   │           └── Metabolic Bottlenecks
   └── Protective Genes Active
```

### Five Stages of Progression

| Stage | CWI Score | Key Events | Molecular Markers |
|-------|-----------|------------|-------------------|
| **Normal** | < 0.3 | Healthy function, demand matches capacity | GCK↑, SLC2A2↑, PDX1↑, MAFA↑ |
| **Compensation** | 0.3-0.6 | Hyperfunction, metabolic strain | INS↑, HSPA5↑, bottlenecks emerge |
| **Stress** | 0.6-0.8 | ER stress, UPR activation | DDIT3↑, XBP1↑, ATF4↑ |
| **Exhaustion** | > 0.8 | Dedifferentiation, identity loss | PDX1↓, MAFA↓, ALDH1A3↑, SOX9↑ |
| **Failure** | Clinical T2D | Beta-cell mass < 50% | Fasting hyperglycemia |

### Composite Workload Index (CWI)

The CWI integrates four dimensions of beta-cell stress:
- **Demand**: Insulin secretion requirements
- **Capacity**: Maximum functional output
- **Stress**: ER stress and oxidative markers
- **Dedifferentiation**: Loss of mature beta-cell identity

---

## 2. Data Sources & Datasets

### Primary Data

| Dataset | Source | Description | Size |
|---------|--------|-------------|------|
| Single-cell RNA-seq | In-house | Human islet transcriptomics | 269 cells |
| CWI Scores | Computed | Workload indices per cell | 269 scores |
| Multi-omics Signatures | DIABLO | Integrated feature selection | 22 features |

### External Databases

| Database | URL | Data Used | Access |
|----------|-----|-----------|--------|
| **Human-GEM** | https://github.com/SysBioChalmers/Human-GEM | Genome-scale metabolic model | Public |
| **DIAGRAM Consortium** | https://diagram-consortium.org | T2D GWAS summary statistics | Public |
| **LINCS L1000** | https://lincsproject.org | Drug perturbation signatures | Public |
| **DisGeNET** | https://disgenet.com | Gene-disease associations | API Key |
| **DisGeNET (GitHub)** | https://github.com/dhimmel/disgenet | Processed GDA data | Public |

### GitHub Repositories Used

```bash
# Human Genome-scale Metabolic Model
git clone https://github.com/SysBioChalmers/Human-GEM.git

# DisGeNET processed data (fallback)
# Source: https://github.com/dhimmel/disgenet
# File: https://raw.githubusercontent.com/dhimmel/disgenet/master/data/consolidated.tsv

# COBRApy (Constraint-Based Reconstruction and Analysis)
pip install cobra
```

---

## 3. Tools & Methods

### Computational Tools

| Tool | Version | Purpose | Reference |
|------|---------|---------|-----------|
| **COBRApy** | 0.30.0 | Metabolic modeling (FBA, FVA, pFBA) | [Ebrahim et al., 2013](https://doi.org/10.1186/1752-0509-7-74) |
| **Human-GEM** | 1.17.0 | Human metabolic reconstruction | [Robinson et al., 2020](https://doi.org/10.1126/scisignal.aaz1482) |
| **Pandas** | 2.x | Data manipulation | - |
| **Matplotlib** | 3.x | Visualization | - |
| **Requests** | 2.x | API interactions | - |

### Analysis Methods

#### 3.1 Flux Balance Analysis (FBA)
```python
from cobra.io import load_model
model = load_model("Human-GEM")
solution = model.optimize()
```

#### 3.2 Mendelian Randomization
- **Method**: Two-sample MR (IVW, Weighted Median, MR-Egger)
- **Exposure**: Gene expression instruments from eQTL
- **Outcome**: T2D risk from DIAGRAM GWAS

#### 3.3 LINCS Connectivity Mapping
- Query workload gene signatures against L1000 database
- Identify compounds that reverse stress signatures
- 473,647 perturbation signatures analyzed

#### 3.4 DisGeNET Integration
- Gene-disease association validation
- Novel target discovery
- 82,833 curated associations from processed data

---

## 4. Analysis Pipeline

### Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                     BETA-CELL WORKLOAD ANALYSIS PIPELINE                │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  [1] DATA INTEGRATION                                                   │
│       │                                                                 │
│       ├── Single-cell RNA-seq (269 cells)                              │
│       ├── CWI Score Computation                                         │
│       └── Multi-omics DIABLO (22 features)                             │
│                     │                                                   │
│  [2] METABOLIC MODELING ─────────────────────────────────────────────  │
│       │                                                                 │
│       ├── Human-GEM (12,971 reactions)                                 │
│       ├── FBA/pFBA/FVA Analysis                                        │
│       └── Bottleneck Identification (43 reactions)                     │
│                     │                                                   │
│  [3] CAUSAL INFERENCE ───────────────────────────────────────────────  │
│       │                                                                 │
│       ├── Mendelian Randomization                                       │
│       └── 5 Validated Targets (PDX1, SLC2A2, MAFA, GCK, INS)          │
│                     │                                                   │
│  [4] DRUG DISCOVERY ─────────────────────────────────────────────────  │
│       │                                                                 │
│       ├── LINCS Connectivity (473,647 signatures)                      │
│       ├── DisGeNET Validation (82,833 associations)                    │
│       └── 24 Drug Candidates Identified                                 │
│                     │                                                   │
│  [5] INTEGRATION & VISUALIZATION ────────────────────────────────────  │
│       │                                                                 │
│       ├── Therapeutic Roadmap                                           │
│       └── Publication-ready Figures                                     │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

### Code Modules

| Script | Purpose | Key Outputs |
|--------|---------|-------------|
| `01_download_gem_models.py` | Download Human-GEM | `Human-GEM.xml` |
| `02_betacell_gem_analysis.py` | Metabolic state analysis | Flux distributions |
| `03_workload_metabolic_integration.py` | CWI + GEM integration | 269 cell profiles |
| `04_gem_drug_integration.py` | Drug target mapping | Therapeutic roadmap |
| `05_disgenet_integration.py` | Disease association validation | MR target validation |
| `create_workload_schematic.py` | Figure generation | Publication figures |

---

## 5. Key Results

### 5.1 Mendelian Randomization Targets

All 5 MR-validated targets confirmed with DisGeNET evidence:

| Gene | MR Effect | Odds Ratio | DisGeNET Score | Validation |
|------|-----------|------------|----------------|------------|
| **SLC2A2** | Protective | 0.834 | 0.558 | ✓ VALIDATED |
| **GCK** | Protective | 0.872 | 0.523 | ✓ VALIDATED |
| **INS** | Complex | 1.012 | 0.443 | ✓ VALIDATED |
| **PDX1** | Protective | 0.857 | 0.344 | ✓ VALIDATED |
| **MAFA** | Risk | 1.135 | 0.101 | ✓ VALIDATED |

### 5.2 Metabolic Bottlenecks

GEM analysis identified 43 metabolic bottleneck reactions:

| Pathway | Bottleneck Reactions | Significance |
|---------|---------------------|--------------|
| Arginine/Proline metabolism | 24 | Amino acid sensing |
| Fructose/Mannose metabolism | 6 | Glycolytic flux |
| Pyrimidine metabolism | 4 | DNA repair capacity |
| Nucleotide metabolism | 4 | Proliferation |
| Alanine/Glutamate metabolism | 3 | Anaplerosis |

### 5.3 Novel Therapeutic Targets

Top 10 novel T2D-associated genes from DisGeNET:

| Rank | Gene | DisGeNET Score | Function |
|------|------|----------------|----------|
| 1 | **IRS1** | 0.907 | Insulin receptor substrate |
| 2 | **HNF4A** | 0.725 | Transcription factor |
| 3 | **AKT2** | 0.721 | Kinase signaling |
| 4 | **HNF1B** | 0.692 | MODY gene |
| 5 | **TCF7L2** | 0.648 | GWAS top hit |
| 6 | **MAPK8IP1** | 0.620 | JNK scaffold |
| 7 | **CAPN10** | 0.616 | Calpain protease |
| 8 | **IL6** | 0.593 | Inflammatory cytokine |
| 9 | **INSR** | 0.590 | Insulin receptor |
| 10 | **SLC2A4** | 0.551 | GLUT4 transporter |

### 5.4 Drug Candidates

24 compounds identified from LINCS connectivity:

**Known T2D Drugs (Benchmarks):**
- Metformin, Glipizide, Glimepiride, Pioglitazone
- Sitagliptin, Repaglinide, Acarbose, etc.

**Novel LINCS Hits (High Priority):**
- RHO-kinase inhibitors (Connectivity: 1.0)
- SRC-kinase inhibitors (Connectivity: 1.0)
- VEGF-receptor inhibitors (Connectivity: 1.0)
- Phenformin (Metformin analog)

---

## 6. Therapeutic Implications

### Tiered Therapeutic Roadmap

| Tier | Priority | Target | Approach | Evidence |
|------|----------|--------|----------|----------|
| **1** | High | GCK | Glucokinase activators | MR + GEM |
| **1** | High | PDX1 | HDAC inhibitors | MR + DisGeNET |
| **1** | High | SLC2A2 | Transporter modulators | MR + DisGeNET |
| **2** | Medium | IRS1 | Signaling enhancers | DisGeNET novel |
| **2** | Medium | HNF4A | TF activators | DisGeNET novel |
| **3** | Benchmark | Multiple | Existing T2D drugs | Clinical validation |

### Stage-Specific Interventions

```
STAGE           INTERVENTION TARGET           THERAPEUTIC APPROACH
─────────────────────────────────────────────────────────────────
Compensation    GCK (glucose sensor)         Glucokinase activators
                                              (e.g., Dorzagliatin)

Stress          ER stress pathways           Chemical chaperones
                XBP1, HSPA5                   HDAC inhibitors

Exhaustion      PDX1, MAFA restoration       Epigenetic modulators
                                              Identity preservation

Failure         Multi-target                  Combination therapy
                                              GLP-1 + SGLT2i + novel
```

---

## 7. File Structure

```
workload/
├── workload.md                          # This documentation
├── code/
│   ├── gem_analysis/
│   │   ├── 01_download_gem_models.py    # Human-GEM download
│   │   ├── 02_betacell_gem_analysis.py  # Metabolic modeling
│   │   ├── 03_workload_metabolic_integration.py
│   │   ├── 04_gem_drug_integration.py   # Drug discovery
│   │   ├── 05_disgenet_integration.py   # Disease associations
│   │   └── create_workload_schematic.py # Figure generation
│   └── mr_analysis/
│       └── 02_mendelian_randomization.py
├── data/
│   └── gem_models/
│       └── Human-GEM.xml                # 43MB metabolic model
├── results/
│   ├── gem_analysis/
│   │   ├── metabolic_bottlenecks.csv
│   │   └── workload_metabolic_integration.csv
│   ├── mr_analysis/
│   │   └── mr_expected_vs_observed.csv
│   ├── drug_discovery/
│   │   └── workload_modulator_candidates.csv
│   ├── disgenet/
│   │   ├── workload_gene_associations.csv
│   │   ├── t2d_gene_associations.csv
│   │   ├── mr_target_validation.csv
│   │   ├── novel_targets.csv
│   │   └── disgenet_report.txt
│   └── integrated_analysis/
│       ├── integration_report.txt
│       └── therapeutic_roadmap.csv
└── figures/
    ├── betacell_workload_model.png      # Main schematic
    └── betacell_workload_model.pdf      # Publication version
```

---

## 8. Reproducibility

### Environment Setup

```bash
# Create conda environment
conda create -n workload python=3.12
conda activate workload

# Install dependencies
pip install cobra pandas numpy matplotlib requests

# Verify COBRApy installation
python -c "import cobra; print(f'COBRApy version: {cobra.__version__}')"
```

### Running the Pipeline

```bash
# 1. Download Human-GEM model
python code/gem_analysis/01_download_gem_models.py

# 2. Run metabolic analysis
python code/gem_analysis/02_betacell_gem_analysis.py

# 3. Integrate with workload scores
python code/gem_analysis/03_workload_metabolic_integration.py

# 4. Drug discovery integration
python code/gem_analysis/04_gem_drug_integration.py

# 5. DisGeNET validation
python code/gem_analysis/05_disgenet_integration.py

# 6. Generate figures
python code/gem_analysis/create_workload_schematic.py
```

### API Keys Required

| Service | Key Location | Purpose |
|---------|--------------|---------|
| DisGeNET | User account | Gene-disease queries |
| OpenRouter | Environment variable | AI figure generation (optional) |

---

## 9. References

### Primary Methods

1. **COBRApy**: Ebrahim A, et al. (2013). COBRApy: COnstraints-Based Reconstruction and Analysis for Python. *BMC Systems Biology*, 7:74.

2. **Human-GEM**: Robinson JL, et al. (2020). An atlas of human metabolism. *Science Signaling*, 13(624):eaaz1482.

3. **DisGeNET**: Piñero J, et al. (2020). The DisGeNET knowledge platform for disease genomics: 2019 update. *Nucleic Acids Research*, 48(D1):D845-D855.

4. **LINCS**: Subramanian A, et al. (2017). A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. *Cell*, 171(6):1437-1452.

### Beta-Cell Biology

5. **Beta-cell dedifferentiation**: Talchai C, et al. (2012). Pancreatic β cell dedifferentiation as a mechanism of diabetic β cell failure. *Cell*, 150(6):1223-1234.

6. **ER stress in diabetes**: Eizirik DL, et al. (2008). The role for endoplasmic reticulum stress in diabetes mellitus. *Endocrine Reviews*, 29(1):42-61.

7. **Beta-cell identity**: Swisa A, et al. (2017). PAX6 maintains β cell identity by repressing genes of alternative islet cell types. *Journal of Clinical Investigation*, 127(1):230-243.

### Mendelian Randomization

8. **MR methods**: Burgess S, et al. (2015). Using published data in Mendelian randomization: a blueprint for efficient identification of causal risk factors. *European Journal of Epidemiology*, 30(7):543-552.

9. **DIAGRAM Consortium**: Mahajan A, et al. (2018). Fine-mapping type 2 diabetes loci to single-variant resolution using high-density imputation and islet-specific epigenome maps. *Nature Genetics*, 50(11):1505-1513.

---

## Acknowledgments

This work integrates publicly available resources:
- Human-GEM from SysBioChalmers (Chalmers University of Technology)
- DisGeNET processed data from Daniel Himmelstein (dhimmel/disgenet)
- LINCS L1000 from the NIH LINCS Program
- DIAGRAM Consortium GWAS data

---

## License

This analysis pipeline is provided for academic research purposes. Please cite the original data sources when using these results.

---

*Last updated: 2026-01-31*
*Contact: Beta-Cell Workload Analysis Pipeline*
