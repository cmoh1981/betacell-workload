# Deep Learning Architecture for Beta-Cell Workload Index

**Date:** 2026-01-30
**Objective:** Derive data-driven workload index using literature-guided deep learning

---

## 1. Overview

### Approach: Literature-Guided Deep Learning

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    WORKLOAD INDEX DEEP LEARNING PIPELINE                    │
└─────────────────────────────────────────────────────────────────────────────┘

┌──────────────────┐     ┌──────────────────┐     ┌──────────────────────────┐
│  LITERATURE      │     │  SINGLE-CELL     │     │  DEEP LEARNING           │
│  PRIOR KNOWLEDGE │     │  DATA            │     │  MODEL                   │
├──────────────────┤     ├──────────────────┤     ├──────────────────────────┤
│ • 60 curated     │     │ • 250,000+ cells │     │ • VAE with attention     │
│   genes          │ ──► │ • T2D/Normal     │ ──► │ • Literature-guided      │
│ • 3 pillars      │     │   labels         │     │   feature selection      │
│ • 5 states       │     │ • 12+ cell types │     │ • Multi-task learning    │
└──────────────────┘     └──────────────────┘     └──────────────────────────┘
                                                              │
                                                              ▼
                                                  ┌──────────────────────────┐
                                                  │  OUTPUTS                 │
                                                  ├──────────────────────────┤
                                                  │ • Continuous CWI score   │
                                                  │ • State classification   │
                                                  │ • Data-driven gene       │
                                                  │   importance             │
                                                  │ • Validated signatures   │
                                                  └──────────────────────────┘
```

---

## 2. Data Sources (From Inventory)

### Training Data
| Dataset | Cells | Purpose |
|---------|-------|---------|
| Segerstolpe (E-MTAB-5061) | ~3,500 | Primary training (labeled T2D/Normal) |
| GSE221156 Atlas | ~245,000 | Large-scale training |
| Beta trajectory | Variable | Pseudotime supervision |

### Validation Data
| Dataset | Cells | Purpose |
|---------|-------|---------|
| GSE81547 (Enge) | ~2,500 | External validation |
| Baron (GSE84133) | ~8,500 | Cross-dataset validation |

### Prior Knowledge
| Source | Content |
|--------|---------|
| Literature signatures | 60 genes across 3 pillars |
| DEG results | Cell-type specific changes |
| TF activities | Regulatory network |

---

## 3. Model Architecture

### 3.1 Literature-Guided Variational Autoencoder (LG-VAE)

```
INPUT LAYER (Gene Expression)
     │
     ▼
┌─────────────────────────────────────────────────────────────────┐
│  LITERATURE-GUIDED FEATURE SELECTION                            │
│  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐ ┌───────────┐  │
│  │ Pillar 1    │ │ Pillar 2    │ │ Pillar 3    │ │ Dediff    │  │
│  │ Biosynthetic│ │ Metabolic   │ │ Stress      │ │ Markers   │  │
│  │ (15 genes)  │ │ (20 genes)  │ │ (18 genes)  │ │ (7 genes) │  │
│  └──────┬──────┘ └──────┬──────┘ └──────┬──────┘ └─────┬─────┘  │
│         │               │               │              │        │
│         └───────────────┴───────┬───────┴──────────────┘        │
│                                 │                               │
│                          ┌──────┴──────┐                        │
│                          │ Attention   │                        │
│                          │ Mechanism   │                        │
│                          └──────┬──────┘                        │
└─────────────────────────────────┼───────────────────────────────┘
                                  │
                                  ▼
                    ┌─────────────────────────┐
                    │      ENCODER            │
                    │  FC(256) → FC(128)      │
                    │  → FC(64)               │
                    └───────────┬─────────────┘
                                │
                    ┌───────────┴───────────┐
                    │    LATENT SPACE       │
                    │    z ~ N(μ, σ²)       │
                    │    dim = 32           │
                    │                       │
                    │  ┌─────────────────┐  │
                    │  │ WORKLOAD DIM    │  │ ← Primary output (CWI)
                    │  │ (first 4 dims)  │  │
                    │  └─────────────────┘  │
                    └───────────┬───────────┘
                                │
          ┌─────────────────────┼─────────────────────┐
          ▼                     ▼                     ▼
┌─────────────────┐   ┌─────────────────┐   ┌─────────────────┐
│    DECODER      │   │  STATE HEAD     │   │  CONDITION HEAD │
│    (Recon)      │   │  (5-class)      │   │  (T2D/Normal)   │
└─────────────────┘   └─────────────────┘   └─────────────────┘
```

### 3.2 Multi-Task Learning Objectives

```python
LOSS = λ₁ × Reconstruction_Loss     # Preserve gene expression
     + λ₂ × KL_Divergence           # Regularize latent space
     + λ₃ × Condition_Loss          # Predict T2D vs Normal
     + λ₄ × State_Loss              # Classify 5 workload states
     + λ₅ × Literature_Alignment    # Align with prior knowledge

Where:
- λ₁ = 1.0   (reconstruction)
- λ₂ = 0.1   (KL regularization)
- λ₃ = 1.0   (condition prediction)
- λ₄ = 0.5   (state classification)
- λ₅ = 0.3   (literature alignment)
```

---

## 4. Feature Engineering

### 4.1 Literature-Based Feature Modules

```python
FEATURE_MODULES = {
    # Module 1: Biosynthetic Demand (15 genes)
    "biosynthetic": [
        "INS", "IAPP", "PCSK1", "PCSK2", "CPE", "CHGB", "SCG2",  # Production
        "HSPA5", "HSP90B1", "PDIA4", "PDIA6", "CALR", "CANX",    # ER folding
        "SLC30A8", "SNAP25"                                       # Secretory
    ],

    # Module 2: Metabolic Capacity (20 genes)
    "metabolic": [
        "GCK", "SLC2A2", "G6PC2", "PFKFB2",                      # Glucose
        "TFAM", "PPARGC1A", "MT-ND1", "MT-CO1", "MT-ATP6", "HADH", # Mito
        "PDX1", "MAFA", "NKX6-1", "UCN3", "NEUROD1", "NKX2-2", "PAX6", # Identity
        "PPARA", "PPARD", "HNF4A"                                 # Lipid
    ],

    # Module 3: Stress Response (18 genes)
    "stress": [
        "XBP1", "ATF6", "ERN1", "EIF2AK3",                       # UPR adaptive
        "DDIT3", "ATF4", "TRIB3", "BBC3", "GADD34",              # UPR terminal
        "NFE2L2", "SOD1", "SOD2", "GPX1", "CAT",                 # Oxidative
        "NFKB1", "TNF", "IL1B", "CCL2"                           # Inflammatory
    ],

    # Module 4: Dedifferentiation (7 genes)
    "dedifferentiation": [
        "ALDH1A3", "NEUROG3", "SOX9", "HES1", "GASTRIN",         # Up in dediff
        "LDHA", "HK1"                                             # Disallowed
    ]
}
```

### 4.2 Data-Driven Feature Expansion

In addition to literature genes, include:
- Top 100 DEGs from `DEG_Beta_T2D_vs_Normal.csv`
- Top 50 genes from `deep_learning_gene_importance.csv`
- All expressed TFs from `tf_activities.h5ad`

**Total features:** ~200-300 genes

---

## 5. Label Generation

### 5.1 Primary Labels (Supervised)

| Label Type | Source | Values |
|------------|--------|--------|
| Condition | Metadata | T2D=1, Normal=0 |
| Disease Stage | Annotation | PreDM=0.5 (if available) |

### 5.2 Pseudo-Labels (Semi-Supervised)

| Label Type | Method | Values |
|------------|--------|--------|
| Workload State | Literature scoring | S1-S5 (5 classes) |
| Pseudotime | Trajectory analysis | 0-1 continuous |
| Stress Level | UPR gene score | Low/Med/High |

### 5.3 Pseudo-Label Generation Algorithm

```python
def generate_workload_pseudolabels(adata, literature_signatures):
    """
    Generate pseudo-labels using literature-based scoring.
    """
    # Score each pillar
    biosyn_score = score_genes(adata, literature_signatures["biosynthetic"])
    metab_score = score_genes(adata, literature_signatures["metabolic"])
    stress_score = score_genes(adata, literature_signatures["stress"])
    dediff_score = score_genes(adata, literature_signatures["dedifferentiation"])

    # Compute CWI (Composite Workload Index)
    cwi = (biosyn_score / (metab_score + 0.1)) * (1 + stress_score) * (1 + dediff_score * 0.5)

    # Normalize to 0-5 range
    cwi_norm = (cwi - cwi.min()) / (cwi.max() - cwi.min()) * 5

    # Assign states
    states = pd.cut(cwi_norm,
                    bins=[-np.inf, 0.5, 1.0, 1.5, 2.5, np.inf],
                    labels=["S1_Resting", "S2_Active", "S3_Stressed",
                           "S4_Exhausted", "S5_Failing"])

    return cwi_norm, states
```

---

## 6. Training Strategy

### 6.1 Two-Phase Training

```
PHASE 1: Pre-training (Unsupervised)
├── Train VAE on all cells (~250k)
├── Learn general gene expression manifold
├── No condition labels used
└── Output: Pre-trained encoder

PHASE 2: Fine-tuning (Semi-supervised)
├── Use beta cells only (~18k)
├── Add condition prediction head
├── Add state classification head
├── Use literature pseudo-labels
└── Output: Final workload model
```

### 6.2 Training Parameters

```python
TRAINING_CONFIG = {
    # Phase 1: Pre-training
    "pretrain": {
        "epochs": 100,
        "batch_size": 256,
        "learning_rate": 1e-3,
        "latent_dim": 32,
        "dropout": 0.2
    },

    # Phase 2: Fine-tuning
    "finetune": {
        "epochs": 50,
        "batch_size": 128,
        "learning_rate": 1e-4,
        "freeze_encoder": False,
        "early_stopping": 10
    },

    # Architecture
    "encoder_layers": [256, 128, 64],
    "decoder_layers": [64, 128, 256],
    "latent_dim": 32,
    "workload_dims": 4,  # First 4 latent dims = workload

    # Loss weights
    "lambda_recon": 1.0,
    "lambda_kl": 0.1,
    "lambda_condition": 1.0,
    "lambda_state": 0.5,
    "lambda_literature": 0.3
}
```

---

## 7. Workload Index Extraction

### 7.1 From Latent Space to CWI

```python
def extract_workload_index(model, adata):
    """
    Extract continuous workload index from trained model.
    """
    # Get latent representation
    z_mean, z_var = model.encode(adata.X)

    # Workload dimensions (first 4)
    workload_latent = z_mean[:, :4]

    # Aggregate to single score
    # Weighted by learned importance
    weights = model.workload_weights  # Learned during training
    cwi = (workload_latent * weights).sum(axis=1)

    # Normalize to 0-5 range
    cwi_normalized = (cwi - cwi.min()) / (cwi.max() - cwi.min()) * 5

    return cwi_normalized
```

### 7.2 Data-Driven Gene Importance

```python
def get_gene_importance(model, gene_names):
    """
    Extract gene importance from attention weights.
    """
    # Get attention weights from first layer
    attention_weights = model.attention_layer.get_weights()

    # Map to gene names
    importance = pd.DataFrame({
        "gene": gene_names,
        "importance": attention_weights,
        "module": assign_to_module(gene_names)  # Pillar assignment
    })

    return importance.sort_values("importance", ascending=False)
```

---

## 8. Validation Strategy

### 8.1 Internal Validation

| Method | Metric | Target |
|--------|--------|--------|
| Reconstruction | MSE | < 0.1 |
| Condition prediction | AUROC | > 0.85 |
| State classification | Accuracy | > 0.75 |
| Cluster separation | Silhouette | > 0.3 |

### 8.2 External Validation

| Dataset | Validation Type |
|---------|-----------------|
| Enge (GSE81547) | Cross-dataset CWI correlation |
| Baron (GSE84133) | State distribution |
| Literature markers | Correlation with known genes |

### 8.3 Literature Alignment Validation

```python
def validate_literature_alignment(cwi, adata):
    """
    Validate CWI aligns with literature expectations.
    """
    validations = {}

    # 1. CWI should be higher in T2D
    if "condition" in adata.obs:
        t2d_cwi = cwi[adata.obs["condition"] == "T2D"].mean()
        normal_cwi = cwi[adata.obs["condition"] == "Normal"].mean()
        validations["t2d_vs_normal"] = t2d_cwi > normal_cwi

    # 2. CWI should correlate negatively with identity markers
    for gene in ["PDX1", "MAFA", "UCN3"]:
        if gene in adata.var_names:
            corr = np.corrcoef(cwi, adata[:, gene].X.flatten())[0, 1]
            validations[f"corr_{gene}"] = corr < 0  # Should be negative

    # 3. CWI should correlate positively with stress markers
    for gene in ["DDIT3", "ATF4", "ALDH1A3"]:
        if gene in adata.var_names:
            corr = np.corrcoef(cwi, adata[:, gene].X.flatten())[0, 1]
            validations[f"corr_{gene}"] = corr > 0  # Should be positive

    return validations
```

---

## 9. Implementation Plan

### Week 1: Data Preparation
- [ ] Load and merge h5ad files
- [ ] Filter to beta cells
- [ ] Generate pseudo-labels from literature
- [ ] Create train/val/test splits

### Week 2: Model Development
- [ ] Implement LG-VAE architecture
- [ ] Add attention mechanism
- [ ] Implement multi-task heads
- [ ] Set up training pipeline

### Week 3: Training
- [ ] Phase 1: Pre-train on all cells
- [ ] Phase 2: Fine-tune on beta cells
- [ ] Hyperparameter optimization
- [ ] Save best model

### Week 4: Validation
- [ ] Internal validation metrics
- [ ] External dataset validation
- [ ] Literature alignment check
- [ ] Gene importance analysis

### Week 5: Analysis
- [ ] Compare literature vs data-driven CWI
- [ ] Identify novel workload genes
- [ ] Refine gene signatures
- [ ] Generate figures

---

## 10. Expected Outputs

### Data Files
```
workload/results/deep_learning/
├── model_weights.pt              # Trained model
├── cwi_scores.csv                # Per-cell CWI
├── state_predictions.csv         # State assignments
├── gene_importance.csv           # Data-driven importance
├── validation_metrics.json       # Performance metrics
└── literature_comparison.csv     # Literature vs data-driven
```

### Key Deliverables
1. **Data-driven CWI** - Continuous workload score per cell
2. **Refined gene signatures** - Top genes from attention weights
3. **Validated states** - 5 workload states with markers
4. **External validation** - Cross-dataset performance

---

*Architecture Version: 1.0*
*Last Updated: 2026-01-30*
