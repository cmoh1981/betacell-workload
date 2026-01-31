# Beta-Cell Workload Project - Data Inventory

**Date:** 2026-01-30
**Source Directory:** `../betacell/`

---

## Summary

| Category | Count | Description |
|----------|-------|-------------|
| **Processed h5ad** | 11 | Ready-to-use AnnData objects |
| **Raw Datasets** | 4 | GEO/ArrayExpress downloads |
| **DEG Results** | 13 | Cell-type specific differential expression |
| **Analysis Results** | 20+ | Multi-omics, trajectories, drug targets |

---

## 1. Processed Single-Cell Data (h5ad)

### Primary Datasets

| File | Location | Description | Cells (est.) |
|------|----------|-------------|--------------|
| `segerstolpe_processed.h5ad` | data/processed/ | E-MTAB-5061, T2D vs Normal | ~3,500 |
| `gse221156_atlas_processed.h5ad` | data/processed/ | Large human islet atlas | ~245,000 |
| `gse221156_atlas.h5ad` | .lobster_workspace/ | Atlas (alternate processing) | ~245,000 |

### Trajectory Analysis

| File | Location | Description |
|------|----------|-------------|
| `beta_trajectory.h5ad` | results/celltype_deep_analysis/ | Beta cell pseudotime |
| `endocrine_trajectory.h5ad` | results/celltype_deep_analysis/ | Endocrine lineage |
| `immune_trajectory.h5ad` | results/celltype_deep_analysis/ | Immune cell trajectory |

### Specialized Analyses

| File | Location | Description |
|------|----------|-------------|
| `tf_activities.h5ad` | results/tf_analysis/ | Transcription factor activities |
| `spatial_validation.h5ad` | results/spatial_integration/ | Spatial transcriptomics |
| `synthetic_pancreas.h5ad` | data/raw/ | Synthetic reference data |

---

## 2. Raw Data Sources

### E-MTAB-5061 (Segerstolpe et al., 2016)
```
Location: data/raw/E-MTAB-5061/
Files:
├── E-MTAB-5061.aggregated_filtered_normalised_counts.mtx (420 MB)
├── E-MTAB-5061.aggregated_filtered_normalised_counts.mtx_cols
├── E-MTAB-5061.aggregated_filtered_normalised_counts.mtx_rows
├── E-MTAB-5061.sdrf.txt (metadata)
└── E-MTAB-5061-normalised-files.zip
```
- **Cells:** ~3,500
- **Condition:** T2D (n=4) vs Normal (n=6)
- **Technology:** Smart-seq2

### GSE221156 (2023 Large Atlas)
```
Location: data/raw/GSE221156/extracted/
Files: 48 samples × 3 files each (barcodes, features, matrix)
├── GSM6846488_MS17001_*.gz
├── GSM6846489_MS17002_*.gz
├── ... (48 donors)
└── GSM6846535_MS19012_*.gz
```
- **Cells:** ~245,000
- **Donors:** 48 (17 T2D, 31 Normal)
- **Technology:** 10x Genomics

### GSE81547 (Enge et al., Aging Study)
```
Location: data/raw/GSE81547/
Files:
├── GSE81547_RAW.tar (176 MB)
├── GSE81547_family.soft.gz
└── extracted/
```
- **Cells:** ~2,500
- **Focus:** Aging pancreas

### GSE84133 (Baron et al.)
```
Location: data/raw/
Files:
└── baron_GSE84133.tar.gz (29 MB)
```
- **Cells:** ~8,500
- **Species:** Human + Mouse

---

## 3. Differential Expression Results

### Cell-Type Specific DEGs (T2D vs Normal)

| File | Cell Type | Top DE Genes |
|------|-----------|--------------|
| `DEG_Beta_T2D_vs_Normal.csv` | Beta | B2M↓, TFF3↑, SIX3↑, NEAT1↓ |
| `DEG_Alpha_T2D_vs_Normal.csv` | Alpha | GCG changes |
| `DEG_Delta_T2D_vs_Normal.csv` | Delta | SST pathway |
| `DEG_PP_T2D_vs_Normal.csv` | PP/Gamma | PPY related |
| `DEG_Epsilon_T2D_vs_Normal.csv` | Epsilon | GHRL pathway |
| `DEG_Acinar_T2D_vs_Normal.csv` | Acinar | Digestive enzymes |
| `DEG_Ductal_T2D_vs_Normal.csv` | Ductal | Transport genes |
| `DEG_Endothelial_T2D_vs_Normal.csv` | Endothelial | Vascular |
| `DEG_Stellate_T2D_vs_Normal.csv` | Stellate | ECM genes |
| `DEG_Immune_T2D_vs_Normal.csv` | Immune | Inflammatory |
| `DEG_Global_T2D_vs_Normal.csv` | All cells | Global changes |
| `DEG_PreDiabetes_vs_Normal.csv` | All | Prediabetes stage |

---

## 4. Multi-Omics Analysis Results

### Driver Gene Analysis
```
Location: results/multiomics_drivers/
├── final_therapeutic_targets_ranked.csv    # Top targets with scores
├── dedifferentiation_drivers_ranked.csv    # Dediff-specific drivers
├── deep_learning_gene_importance.csv       # DL feature importance
├── deep_learning_model_performance.csv     # Model metrics
└── multiomics_convergence.csv              # Cross-modal agreement
```

### Drug Discovery
```
Location: results/drug_validation/
├── final_drug_rankings.csv
├── integrated_candidate_ranking.csv
├── docking_results.csv
├── admet_predictions.csv
├── lincs_admet_results.csv
└── drug_target_network.csv
```

### Pathway Analysis
```
Location: results/final_integrated_analysis/
├── pathway_changes_by_celltype.csv
├── pathway_deep_analysis.csv
├── cell_cell_communication.csv
└── communication_changes_T2D.csv
```

---

## 5. Data for Deep Learning Workload Index

### Training Data (Available)

| Purpose | Files | Size |
|---------|-------|------|
| **Primary training** | segerstolpe_processed.h5ad | ~3,500 cells |
| **Large-scale training** | gse221156_atlas_processed.h5ad | ~245,000 cells |
| **Feature importance** | deep_learning_gene_importance.csv | Gene weights |
| **Labels (pseudo)** | DEG_Beta_T2D_vs_Normal.csv | T2D/Normal status |

### Validation Data (Available)

| Purpose | Files |
|---------|-------|
| **Cross-validation** | beta_trajectory.h5ad |
| **External validation** | GSE81547 (Enge aging) |
| **Spatial validation** | spatial_validation.h5ad |

### Literature-Based Prior Knowledge

| Source | File |
|--------|------|
| Gene signatures | `workload/docs/BETA_CELL_WORKLOAD_ARCHITECTURE.md` |
| Literature index | `workload/references/LITERATURE_INDEX.md` |

---

## 6. Estimated Total Data

| Metric | Value |
|--------|-------|
| **Total cells (processed)** | ~250,000+ |
| **Total cells (raw available)** | ~260,000+ |
| **Unique donors** | ~70 |
| **T2D donors** | ~21 |
| **Normal donors** | ~40+ |
| **Genes measured** | ~20,000 |
| **Cell types annotated** | 12+ |

---

## 7. Data Quality Notes

### Strengths
- Multiple independent datasets
- Both healthy and T2D conditions
- Cell-type annotations available
- Trajectory analysis completed
- Multi-omics integration done

### Gaps (May Need External Data)
- Limited prediabetes samples
- No ATAC-seq (chromatin accessibility)
- No proteomics at single-cell level
- Limited longitudinal data

---

## 8. Recommended Deep Learning Strategy

### Phase 1: Use Existing Data
1. Load `segerstolpe_processed.h5ad` + `gse221156_atlas_processed.h5ad`
2. Filter to beta cells (~18,000-20,000 cells)
3. Use T2D/Normal as labels
4. Use literature gene signatures as features

### Phase 2: Train Workload Index Model
1. Variational Autoencoder (VAE) for embedding
2. Supervised layer using T2D status
3. Extract workload scores from latent space

### Phase 3: Validation
1. Cross-validate on Enge aging data
2. Compare with literature-defined states
3. Correlate with known markers

---

*Inventory Version: 1.0*
*Last Updated: 2026-01-30*
