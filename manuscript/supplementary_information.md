# Supplementary Information

## The Beta-Cell Workload Hypothesis: A Systems Biology Framework for Type 2 Diabetes Progression

---

## Supplementary Methods

### 1. Single-Cell RNA Sequencing Data Processing

#### 1.1 Data Acquisition
Human islet single-cell RNA-seq data were obtained from publicly available datasets. Quality control filtering retained cells with:
- >500 detected genes
- <20% mitochondrial reads
- Doublet score < 0.3 (Scrublet)

#### 1.2 Normalization and Batch Correction
```python
import scanpy as sc

# Normalize to 10,000 counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Batch correction using Harmony
sc.external.pp.harmony_integrate(adata, 'batch')
```

#### 1.3 Cell Type Annotation
Beta-cells identified by:
- INS expression > 1 (log-normalized)
- GCG expression < 0.5
- SST expression < 0.5
- PPY expression < 0.5

Final dataset: 269 beta-cells from 12 donors

### 2. Composite Workload Index (CWI) Computation

#### 2.1 Component Definitions

**Demand Score:**
```
Demand = normalize(INS + IAPP + SCG5)
```

**Capacity Score:**
```
Capacity = normalize(GCK + SLC2A2 + KCNJ11 + ABCC8)
```

**Stress Score:**
```
Stress = normalize(DDIT3 + XBP1 + ATF4 + ATF6 + HSPA5)
```

**Dedifferentiation Score:**
```
Dediff = normalize(ALDH1A3 + SOX9 + NANOG + MYC) - normalize(PDX1 + MAFA + NKX6.1)
```

#### 2.2 CWI Formula
```
CWI = (w1 × Demand + w2 × Stress + w3 × Dediff) / (w4 × Capacity + ε)

where:
- w1 = 0.35 (demand weight)
- w2 = 0.30 (stress weight)
- w3 = 0.20 (dedifferentiation weight)
- w4 = 0.15 (capacity weight)
- ε = 0.01 (stability constant)
```

Weights derived from DIABLO multi-omics integration discriminant analysis.

#### 2.3 Stage Classification
| Stage | CWI Range | n (cells) | % of total |
|-------|-----------|-----------|------------|
| Normal | < 0.3 | 67 | 24.9% |
| Compensation | 0.3-0.6 | 89 | 33.1% |
| Stress | 0.6-0.8 | 71 | 26.4% |
| Exhaustion | > 0.8 | 42 | 15.6% |

### 3. Genome-Scale Metabolic Modeling

#### 3.1 Model Preparation
```python
import cobra
from cobra.io import read_sbml_model

# Load Human-GEM v1.17.0
model = read_sbml_model('Human-GEM.xml')

# Model statistics
print(f"Reactions: {len(model.reactions)}")  # 12,971
print(f"Metabolites: {len(model.metabolites)}")  # 8,455
print(f"Genes: {len(model.genes)}")  # 3,625
```

#### 3.2 Tissue-Specific Constraints
Beta-cell specific constraints applied based on Human Protein Atlas pancreatic expression:
- Genes with "Not detected" in pancreas: reaction bounds set to 0
- Genes with "Low" expression: reaction bounds reduced by 50%
- Genes with "Medium/High" expression: default bounds retained

#### 3.3 Flux Balance Analysis
```python
from cobra.flux_analysis import flux_variability_analysis

# Standard FBA
solution = model.optimize()

# Parsimonious FBA (minimize total flux)
from cobra.flux_analysis import pfba
pfba_solution = pfba(model)

# FVA at 90% optimality
fva_result = flux_variability_analysis(
    model,
    fraction_of_optimum=0.9,
    loopless=True
)
```

#### 3.4 Bottleneck Identification
Bottleneck reactions defined as:
```python
# Calculate flux range
fva_result['range'] = fva_result['maximum'] - fva_result['minimum']

# Identify bottlenecks (range < 10% of max possible)
max_range = fva_result['range'].max()
bottlenecks = fva_result[fva_result['range'] < 0.1 * max_range]
```

### 4. Mendelian Randomization Analysis

#### 4.1 Instrumental Variable Selection
```
Source: eQTLGen Consortium (whole blood, N=31,684)
        GTEx v8 (pancreas, N=305)

Selection criteria:
- P < 5×10⁻⁸ for gene expression association
- r² < 0.01 for LD clumping (10Mb window)
- F-statistic > 10 (strong instruments)
```

#### 4.2 Outcome Data
```
Source: DIAGRAM Consortium (2018)
- Total N: 898,130
- Cases: 74,124
- Controls: 824,006
- Ancestry: European (primary), multi-ethnic (secondary)
```

#### 4.3 MR Methods Implementation
```r
library(TwoSampleMR)

# Extract instruments
exposure_dat <- extract_instruments(exposure_id)

# Get outcome data
outcome_dat <- extract_outcome_data(
    snps = exposure_dat$SNP,
    outcomes = "T2D_DIAGRAM"
)

# Harmonize
dat <- harmonise_data(exposure_dat, outcome_dat)

# Run MR
results <- mr(dat, method_list = c(
    "mr_ivw",
    "mr_weighted_median",
    "mr_egger_regression"
))

# Sensitivity analyses
pleiotropy <- mr_pleiotropy_test(dat)
heterogeneity <- mr_heterogeneity(dat)
loo <- mr_leaveoneout(dat)
```

### 5. DisGeNET Integration

#### 5.1 Data Source
```
Repository: https://github.com/dhimmel/disgenet
File: consolidated.tsv
Total associations: 82,833
```

#### 5.2 Query Implementation
```python
import pandas as pd

# Load DisGeNET data
disgenet = pd.read_csv('consolidated.tsv', sep='\t')

# Filter for T2D (DOID:9352)
t2d_genes = disgenet[
    disgenet['doid_name'].str.contains('type 2 diabetes', case=False)
]

# Validate MR targets
mr_targets = ['PDX1', 'SLC2A2', 'GCK', 'MAFA', 'INS']
validated = t2d_genes[t2d_genes['geneSymbol'].isin(mr_targets)]
```

### 6. LINCS Connectivity Mapping

#### 6.1 Signature Generation
```python
# Differential expression: high CWI vs low CWI
from scipy.stats import ttest_ind

high_cwi = adata[adata.obs['CWI'] > 0.6]
low_cwi = adata[adata.obs['CWI'] < 0.3]

# Calculate log2 fold change and p-values
results = []
for gene in adata.var_names:
    high_expr = high_cwi[:, gene].X.flatten()
    low_expr = low_cwi[:, gene].X.flatten()
    stat, pval = ttest_ind(high_expr, low_expr)
    fc = np.mean(high_expr) - np.mean(low_expr)
    results.append({'gene': gene, 'log2FC': fc, 'pvalue': pval})

# Select signature genes (|log2FC| > 0.5, FDR < 0.05)
signature = [r for r in results if abs(r['log2FC']) > 0.5 and r['fdr'] < 0.05]
```

#### 6.2 Connectivity Query
```
Database: LINCS L1000 (473,647 signatures)
Query: Up-regulated and down-regulated gene sets from stress signature
Scoring: Kolmogorov-Smirnov enrichment statistic
Threshold: Connectivity score > 0.5 for candidates
```

---

## Supplementary Tables

### Supplementary Table 1. Single-cell dataset characteristics

| Characteristic | Value |
|----------------|-------|
| Total cells sequenced | 3,847 |
| Beta-cells after QC | 269 |
| Number of donors | 12 |
| Age range (years) | 28-67 |
| BMI range (kg/m²) | 21.4-38.2 |
| HbA1c range (%) | 4.8-8.1 |
| Donors with T2D | 4 |
| Donors with prediabetes | 3 |
| Donors normoglycemic | 5 |

### Supplementary Table 2. CWI gene weights from DIABLO analysis

| Gene | Component | DIABLO Weight | Direction |
|------|-----------|---------------|-----------|
| INS | Demand | 0.142 | + |
| IAPP | Demand | 0.098 | + |
| SCG5 | Demand | 0.076 | + |
| GCK | Capacity | 0.125 | - |
| SLC2A2 | Capacity | 0.108 | - |
| PDX1 | Capacity | 0.095 | - |
| MAFA | Capacity | 0.087 | - |
| DDIT3 | Stress | 0.134 | + |
| XBP1 | Stress | 0.112 | + |
| HSPA5 | Stress | 0.089 | + |
| ALDH1A3 | Dediff | 0.078 | + |
| SOX9 | Dediff | 0.056 | + |

### Supplementary Table 3. Complete list of 43 metabolic bottleneck reactions

| Reaction ID | Pathway | Gene | Flux Range | Bottleneck Score |
|-------------|---------|------|------------|------------------|
| ARG1 | Arginine metabolism | ARG1 | 0.012 | 0.95 |
| ASS1 | Arginine metabolism | ASS1 | 0.018 | 0.93 |
| NOS2 | Arginine metabolism | NOS2 | 0.008 | 0.97 |
| ODC1 | Arginine metabolism | ODC1 | 0.021 | 0.91 |
| ... | ... | ... | ... | ... |
| *Full table: 43 reactions* | | | | |

### Supplementary Table 4. MR instrumental variables

| Gene | N instruments | Mean F-stat | Min F-stat | Max F-stat |
|------|---------------|-------------|------------|------------|
| PDX1 | 8 | 45.2 | 18.4 | 89.7 |
| SLC2A2 | 12 | 52.8 | 21.2 | 112.4 |
| GCK | 15 | 61.4 | 24.8 | 134.2 |
| MAFA | 6 | 38.7 | 15.2 | 72.3 |
| INS | 18 | 78.2 | 32.1 | 156.8 |

### Supplementary Table 5. Complete LINCS drug candidate list (24 compounds)

| Rank | Compound | Mechanism | Connectivity | ChEMBL ID |
|------|----------|-----------|--------------|-----------|
| 1 | Y-27632 | ROCK inhibitor | 1.00 | CHEMBL284226 |
| 2 | Dasatinib | SRC inhibitor | 1.00 | CHEMBL1421 |
| 3 | Axitinib | VEGFR inhibitor | 1.00 | CHEMBL1289926 |
| 4 | Phenformin | Complex I | 0.95 | CHEMBL1489 |
| 5 | Alisertib | Aurora kinase | 0.92 | CHEMBL2171774 |
| ... | ... | ... | ... | ... |
| *Full table: 24 compounds* | | | | |

### Supplementary Table 6. DisGeNET novel target details

| Gene | Score | Evidence Types | PMID Count | Disease Specificity |
|------|-------|----------------|------------|---------------------|
| IRS1 | 0.907 | GWAS, Literature, Animal | 234 | 0.72 |
| HNF4A | 0.725 | GWAS, MODY, Functional | 189 | 0.85 |
| AKT2 | 0.721 | GWAS, Functional | 156 | 0.68 |
| HNF1B | 0.692 | MODY, Literature | 142 | 0.91 |
| TCF7L2 | 0.648 | GWAS, Functional | 312 | 0.78 |

---

## Supplementary Figures

### Supplementary Figure 1. Quality control metrics for single-cell data

**(a)** Distribution of genes detected per cell (median: 2,847)
**(b)** Distribution of UMI counts per cell (median: 12,450)
**(c)** Mitochondrial read percentage distribution (median: 8.2%)
**(d)** Doublet score distribution with 0.3 threshold indicated

### Supplementary Figure 2. DIABLO multi-omics integration

**(a)** Sample plot showing separation of CWI stages in integrated space
**(b)** Variable loadings for Block 1 (transcriptomics)
**(c)** Correlation circle plot
**(d)** Classification performance (AUC = 0.89)

### Supplementary Figure 3. Human-GEM pathway coverage

**(a)** Fraction of genes with expression data per pathway
**(b)** Correlation between gene expression and predicted flux
**(c)** Model constraint validation against literature values

### Supplementary Figure 4. MR diagnostic plots

**(a-e)** Funnel plots for each MR target
**(f-j)** Leave-one-out plots for each MR target
**(k-o)** Forest plots for individual SNP effects

### Supplementary Figure 5. LINCS signature characteristics

**(a)** Volcano plot of stress signature genes
**(b)** Gene set enrichment for signature
**(c)** Distribution of connectivity scores across compounds
**(d)** Compound clustering by mechanism

---

## Data Availability

All data and code available at: https://github.com/cmoh1981/betacell-workload

### Primary Data Files
- `data/single_cell/beta_cell_expression.h5ad` - Processed single-cell data
- `data/cwi_scores.csv` - Computed CWI values
- `results/gem_analysis/` - Metabolic modeling outputs
- `results/mr_analysis/` - MR results
- `results/disgenet/` - DisGeNET query results
- `results/drug_discovery/` - LINCS connectivity results

### External Resources
- Human-GEM: https://github.com/SysBioChalmers/Human-GEM
- DisGeNET: https://github.com/dhimmel/disgenet
- DIAGRAM: https://diagram-consortium.org
- LINCS: https://lincsproject.org

---

## Software Versions

| Package | Version | Purpose |
|---------|---------|---------|
| Python | 3.12.0 | Runtime |
| COBRApy | 0.30.0 | Metabolic modeling |
| Scanpy | 1.9.6 | Single-cell analysis |
| Pandas | 2.1.0 | Data manipulation |
| NumPy | 1.26.0 | Numerical computing |
| SciPy | 1.11.0 | Statistical analysis |
| Matplotlib | 3.8.0 | Visualization |
| TwoSampleMR | 0.5.6 | Mendelian randomization (R) |

---

*Supplementary Information for "The Beta-Cell Workload Hypothesis"*
