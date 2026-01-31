# Beta-Cell Workload Architecture (Literature-Based)

## Defining Beta-Cell Workload During Hyperglycemia and Diabetes Progression

**Date:** 2026-01-30
**Version:** 2.0 (Literature-Integrated)
**Project:** Beta-Cell T2D Workload Analysis
**Literature Base:** `references/LITERATURE_INDEX.md` (17 primary sources)

---

## Literature Foundation

This architecture is derived from systematic literature review:

| Domain | Key Sources | Evidence Level |
|--------|-------------|----------------|
| Workload Definition | Prentki & Nolan JCI 2006 | Foundational |
| Disease Progression | Weir & Bonner-Weir Diabetes 2004; Diabetologia 2025 | Clinical + Mechanistic |
| ER Stress/UPR | Front Endocrinol 2021; IJMS 2022; Diabetes ADA 2021 | Molecular |
| Dedifferentiation | Nature EMM 2023; Cell Metab 2014; Nutrients 2021 | Cellular |
| Identity Markers | iScience 2022; Nat Metab 2023; PLOS Genet 2013 | Single-cell |
| Mitochondria | Nat Commun 2019; Cells 2025 | Metabolic |

---

## 1. Definition of Beta-Cell Workload

### 1.1 Literature-Based Definition

From **Prentki & Nolan, JCI 2006 (PMC5021190)**:

> *"The magnitude of the increased demand for insulin due to insulin resistance caused by excess caloric intake and physical inactivity exceeds the magnitude of β-cell mass expansion, resulting in an increase in β-cell workload. In individuals who are susceptible to T2DM, increased β-cell workload may lead to β-cell failure and the development of T2DM."*

From **Diabetologia 2025**:

> *"Early abnormalities in insulin secretion, rather than a reduction in beta cell mass, play a fundamental and primary role in early type 2 diabetes pathogenesis."*

### 1.2 Operational Definition

**Beta-cell workload** = The ratio of **insulin secretory demand** to **functional capacity**

```
                    ┌─────────────────────────────────────────┐
                    │     BETA-CELL WORKLOAD = DEMAND/CAPACITY │
                    └─────────────────────────────────────────┘
                                        │
          ┌─────────────────────────────┴─────────────────────────────┐
          ▼                                                           ▼
┌─────────────────────────┐                           ┌─────────────────────────┐
│   DEMAND (Numerator)    │                           │  CAPACITY (Denominator) │
├─────────────────────────┤                           ├─────────────────────────┤
│ • Insulin resistance    │                           │ • β-cell mass           │
│ • Chronic hyperglycemia │                           │ • ER folding capacity   │
│ • Lipotoxicity          │                           │ • Mitochondrial function│
│ • Decreased β-cell mass │                           │ • Antioxidant defense   │
│                         │                           │ • Identity/maturity     │
│ Evidence: JCI 2006      │                           │ Evidence: Nat Commun    │
└─────────────────────────┘                           └─────────────────────────┘
```

---

## 2. Five-Stage Progression Model

### Based on Weir & Bonner-Weir, Diabetes 2004 (PMID: 15561919)

> *"Five stages of evolving beta-cell dysfunction during progression to diabetes"*

| Stage | Name | Glucose | Insulin | β-Cell Status | Molecular State |
|-------|------|---------|---------|---------------|-----------------|
| **1** | Compensation | Normal | ↑↑ | Hyperfunction | Adaptive UPR |
| **2** | Stable Adaptation | Mild ↑ | ↑↑↑ | Stressed but stable | Chronic UPR |
| **3** | Early Decompensation | ↑↑ | Failing | Dysfunctional | Terminal UPR |
| **4** | Stable Decompensation | ↑↑↑ | ↓ | Dedifferentiated | Identity loss |
| **5** | Severe Decompensation | ↑↑↑↑ | ↓↓ | Failing/Apoptotic | Cell death |

### Clinical Timeline Evidence

From **Postgraduate Medicine 2020**:
> *"A T2D duration of 10 years appears to be the point at which beta-cell loss becomes irreversible."*

```
              TIME COURSE OF β-CELL FAILURE

Years:     0        5        10       15       20
           │        │        │        │        │
Stage 1 ───┼────────┤        │        │        │
(Comp.)    │        │        │        │        │
           │        │        │        │        │
Stage 2 ───│────────┼────────┤        │        │
(Stable)   │        │        │        │        │
           │        │        │        │        │
Stage 3 ───│────────│────────┼────────┤        │
(Early)    │        │        │        │        │
           │        │        │   POINT OF      │
Stage 4 ───│────────│────────│───NO RETURN─────┤
(Decomp.)  │        │        │        │        │
           │        │        │        │        │
Stage 5 ───│────────│────────│────────┼────────┤
(Severe)   │        │        │        │        │
           ▼        ▼        ▼        ▼        ▼
        NORMAL  PREDIAB   T2D      ADVANCED  FAILURE
```

---

## 3. Three-Pillar Molecular Framework

### Based on integrated literature findings

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                     BETA-CELL WORKLOAD PILLARS                               │
├─────────────────────┬─────────────────────┬──────────────────────────────────┤
│   PILLAR 1          │   PILLAR 2          │   PILLAR 3                       │
│   BIOSYNTHETIC      │   METABOLIC         │   STRESS RESPONSE                │
│   DEMAND            │   CAPACITY          │                                  │
├─────────────────────┼─────────────────────┼──────────────────────────────────┤
│ Literature:         │ Literature:         │ Literature:                      │
│ Front Endocrinol    │ Nat Commun 2019     │ IJMS 2022; Diabetes 2021         │
│ 2024; JCI 2019      │ PLOS Genet 2013     │ Front Endocrinol 2021            │
├─────────────────────┼─────────────────────┼──────────────────────────────────┤
│ • Proinsulin synth  │ • Glucose sensing   │ • UPR (adaptive/terminal)        │
│ • ER protein folding│ • OXPHOS/ATP        │ • Oxidative stress               │
│ • Insulin processing│ • β-cell identity   │ • Inflammation                   │
│ • Vesicle secretion │ • Lipid metabolism  │ • Senescence                     │
└─────────────────────┴─────────────────────┴──────────────────────────────────┘
```

---

## 4. Literature-Derived Gene Signatures

### 4.1 PILLAR 1: Biosynthetic Demand

**Source:** Front Endocrinol 2024; JCI 2019

| Module | Genes | Literature Evidence |
|--------|-------|---------------------|
| **Insulin Production** | `INS`, `IAPP`, `PCSK1`, `PCSK2`, `CPE`, `CHGB`, `SCG2` | "Massive insulin storage, with huge biosynthetic capability" |
| **ER Folding** | `HSPA5` (BiP), `HSP90B1`, `PDIA4`, `PDIA6`, `CALR`, `CANX`, `DNAJB11` | "BiP/GRP78 binds UPR sensors under normal conditions" |
| **Secretory** | `SLC30A8`, `SNAP25`, `VAMP2`, `STX1A`, `SYT7`, `RAB3A` | "Vesicle exocytosis machinery" |

### 4.2 PILLAR 2: Metabolic Capacity

**Source:** Nat Commun 2019; PLOS Genet 2013; iScience 2022

| Module | Genes | Literature Evidence |
|--------|-------|---------------------|
| **Glucose Sensing** | `GCK`, `SLC2A2`, `G6PC2`, `PFKFB2` | "Glucose enters β-cells and is metabolized via glycolysis" |
| **Mitochondrial** | `TFAM`, `PPARGC1A`, `MT-ND1`, `MT-CO1`, `MT-ATP6`, `HADH` | "Multiple genes in OXPHOS are downregulated in T2D" |
| **β-Cell Identity** | `PDX1`, `MAFA`, `NKX6-1`, `UCN3`, `NEUROD1`, `NKX2-2`, `PAX6` | "Nkx6.1 is both necessary and sufficient to specify β-cells" |
| **Lipid Metabolism** | `SREBF1`, `PPARA`, `PPARD`, `HNF4A` | "Activates β-oxidation genes in response to elevated fatty acids" |

**Identity Marker Hierarchy (iScience 2022):**
> *"The expression levels of Mafa, Nkx6-1, Pdx1, and Ucn3 were significantly downregulated in the diabetic beta cell subcluster."*

```
β-CELL IDENTITY HIERARCHY (Literature-based)

CORE (Essential):        PDX1 → MAFA → NKX6-1 → NKX2-2
                              ↓
MATURATION:              UCN3 → SLC2A2 → GCK → G6PC2
                              ↓
FUNCTION:                INS → IAPP → PCSK1/2 → CPE
```

### 4.3 PILLAR 3: Stress Response

**Source:** IJMS 2022; Diabetes ADA 2021; Front Endocrinol 2021

#### UPR Pathway Components (Literature-derived)

From **IJMS 2022 (PMC9104816)**:
> *"The UPR is initiated by three transmembrane proteins: PERK, IRE1, and ATF6. In the absence of stress, these exist as inactive monomers bound to BiP/GRP78."*

| Sensor | Gene | Downstream | Function | Phase |
|--------|------|------------|----------|-------|
| **IRE1α** | `ERN1` | `XBP1s` | ↑ Chaperones, ERAD | Adaptive |
| **PERK** | `EIF2AK3` | `ATF4` → `DDIT3` | ↓ Translation, ↑ stress genes | Both |
| **ATF6** | `ATF6` | Cleaved ATF6 | ↑ ER biogenesis | Adaptive |

#### UPR Gene Signatures

| Category | Genes | Role | Citation |
|----------|-------|------|----------|
| **Adaptive UPR** | `XBP1`, `ATF6`, `ERN1`, `HSPA5`, `PDIA4` | Protective | Diabetes 2021 |
| **Terminal UPR** | `DDIT3` (CHOP), `ATF4`, `TRIB3`, `BBC3`, `GADD34` | Pro-apoptotic | IJMS 2022 |
| **Oxidative** | `NFE2L2`, `SOD1`, `SOD2`, `GPX1`, `CAT`, `TXN` | Antioxidant | WJD 2023 |
| **Inflammatory** | `NFKB1`, `TNF`, `IL1B`, `CCL2`, `IL6` | Immune | PMC10075035 |
| **Senescence** | `CDKN2A`, `CDKN1A`, `RB1` | Cell cycle arrest | IJMS 2022 |

**β-Cell Vulnerability (WJD 2023):**
> *"Pancreatic β cells are known to exhibit intrinsically low intracellular antioxidative capacity, making them especially vulnerable to ROS-induced injury."*

### 4.4 Dedifferentiation Markers

**Source:** Nature EMM 2023; Nutrients 2021; Cell Metab 2014

From **Nutrients 2021 (PMC8151793)**:
> *"Recent studies suggest that β-cell failure occurs mainly due to increased β-cell dedifferentiation rather than limited β-cell proliferation or increased β-cell death."*

| Category | Genes | Direction | Evidence |
|----------|-------|-----------|----------|
| **DOWNREGULATED** | `PDX1`, `MAFA`, `NKX6-1`, `UCN3`, `SLC2A2`, `NEUROD1`, `FOXO1` | ↓ in T2D | iScience 2022 |
| **UPREGULATED** | `ALDH1A3`, `NEUROG3`, `SOX9`, `HES1`, `GASTRIN` | ↑ Progenitor | Nature EMM 2023 |
| **Disallowed** | `LDHA`, `HK1`, `MCT1` | ↑ in dysfunction | Literature review |

---

## 5. Workload States Classification

### Based on Literature Integration

| State | Name | CWI | Identity Markers | Stress Markers | Literature Basis |
|-------|------|-----|------------------|----------------|------------------|
| **S1** | Resting | <0.5 | UCN3++, MAFA++, PDX1++ | Low UPR | Normal β-cells (iScience 2022) |
| **S2** | Active | 0.5-1.0 | INS++, GCK+, PCSK1+ | Minimal | Weir Stage 1-2 |
| **S3** | Stressed | 1.0-1.5 | Identity intact, HSPA5+, XBP1+ | Adaptive UPR | Weir Stage 2-3 |
| **S4** | Exhausted | 1.5-2.5 | PDX1↓, MAFA↓, DDIT3+, ATF4++ | Terminal UPR | Weir Stage 3-4 |
| **S5** | Failing | >2.5 | ALDH1A3+, PDX1-, UCN3- | Chronic/Senescent | Weir Stage 4-5 |

### State-Specific Gene Signatures (Literature-Curated)

```python
# LITERATURE-BASED STATE MARKERS

STATE_MARKERS = {
    "S1_Resting": {
        "up": ["UCN3", "MAFA", "SLC2A2", "GCK", "PDX1"],
        "down": ["DDIT3", "ATF4", "ALDH1A3"],
        "source": "iScience 2022 - mature β-cell cluster"
    },
    "S2_Active": {
        "up": ["INS", "IAPP", "PCSK1", "PCSK2", "GCK"],
        "down": ["ALDH1A3", "NEUROG3"],
        "source": "Normal GSIS phenotype"
    },
    "S3_Stressed": {
        "up": ["HSPA5", "XBP1", "ATF6", "PDIA4", "ERN1"],
        "down": ["UCN3"],  # partial loss
        "source": "IJMS 2022 - adaptive UPR"
    },
    "S4_Exhausted": {
        "up": ["DDIT3", "ATF4", "TRIB3", "BBC3", "CDKN1A"],
        "down": ["MAFA", "PDX1", "INS", "UCN3"],
        "source": "Front Endocrinol 2021 - terminal UPR"
    },
    "S5_Failing": {
        "up": ["ALDH1A3", "GASTRIN", "SOX9", "NEUROG3", "HES1"],
        "down": ["PDX1", "MAFA", "NKX6-1", "UCN3", "NEUROD1"],
        "source": "Nature EMM 2023 - dedifferentiation"
    }
}
```

---

## 6. Composite Workload Index (CWI)

### Formula Based on Literature Concepts

```
CWI = (Biosynthetic_Demand / Metabolic_Capacity) × Stress_Factor × (1 + Dediff_Penalty)

Where:
- Biosynthetic_Demand = Mean expression of insulin production + ER load genes
- Metabolic_Capacity = Mean expression of identity + metabolic genes
- Stress_Factor = 1 + (Terminal_UPR - Adaptive_UPR) normalized
- Dediff_Penalty = Expression of ALDH1A3, NEUROG3, etc.
```

### Interpretation (Literature-Based)

| CWI Range | State | Prognosis | Evidence |
|-----------|-------|-----------|----------|
| < 1.0 | Compensated | Stable | JCI 2006: "adequate compensation" |
| 1.0 - 2.0 | Stressed | Reversible | Nature EMM 2023: "therapeutic window" |
| > 2.0 | Decompensated | Progressive | Postgrad Med: "10-year threshold" |

---

## 7. Final Gene Signature Summary

### MASTER GENE LIST (Literature-Curated, 60 genes)

```python
WORKLOAD_SIGNATURES_V2 = {
    # PILLAR 1: BIOSYNTHETIC DEMAND (15 genes)
    "biosynthetic_demand": {
        "insulin_production": ["INS", "IAPP", "PCSK1", "PCSK2", "CPE", "CHGB", "SCG2"],
        "er_folding": ["HSPA5", "HSP90B1", "PDIA4", "PDIA6", "CALR", "CANX"],
        "secretory": ["SLC30A8", "SNAP25"]
    },

    # PILLAR 2: METABOLIC CAPACITY (20 genes)
    "metabolic_capacity": {
        "glucose_sensing": ["GCK", "SLC2A2", "G6PC2", "PFKFB2"],
        "mitochondrial": ["TFAM", "PPARGC1A", "MT-ND1", "MT-CO1", "MT-ATP6", "HADH"],
        "beta_identity": ["PDX1", "MAFA", "NKX6-1", "UCN3", "NEUROD1", "NKX2-2", "PAX6"],
        "lipid": ["PPARA", "PPARD", "HNF4A"]
    },

    # PILLAR 3: STRESS RESPONSE (18 genes)
    "stress_response": {
        "upr_adaptive": ["XBP1", "ATF6", "ERN1", "EIF2AK3"],
        "upr_terminal": ["DDIT3", "ATF4", "TRIB3", "BBC3", "GADD34"],
        "oxidative": ["NFE2L2", "SOD1", "SOD2", "GPX1", "CAT"],
        "inflammatory": ["NFKB1", "TNF", "IL1B", "CCL2"]
    },

    # DEDIFFERENTIATION (7 genes)
    "dedifferentiation": {
        "progenitor_up": ["ALDH1A3", "NEUROG3", "SOX9", "HES1", "GASTRIN"],
        "disallowed": ["LDHA", "HK1"]
    }
}
```

---

## 8. Analysis Pipeline

### Phase 1: Data Integration
- Load Segerstolpe, GSE221156, Enge datasets
- Harmony batch correction
- Filter to β-cells (~18,000 cells)

### Phase 2: Workload Scoring
- Compute Pillar 1-3 scores per cell
- Calculate CWI
- Classify into 5 states

### Phase 3: Validation
- Compare Normal vs T2D state distributions
- Validate markers match literature
- Cross-dataset reproducibility

### Phase 4: Target Discovery
- Genes correlated with high CWI
- Pathway enrichment (ER stress, OXPHOS, identity)
- Druggable target identification

---

## 9. References

### Primary Sources (Cited in Architecture)

1. **Prentki M, Nolan CJ.** Islet beta cell failure in type 2 diabetes. *J Clin Invest.* 2006;116(7):1802-1812. [PMC5021190](https://pmc.ncbi.nlm.nih.gov/articles/PMC5021190/)

2. **Weir GC, Bonner-Weir S.** Five stages of evolving beta-cell dysfunction during progression to diabetes. *Diabetes.* 2004;53 Suppl 3:S16-21. PMID: 15561919

3. **Diabetologia 2025.** The role of the beta cell in type 2 diabetes: new findings from the last 5 years. [Link](https://link.springer.com/article/10.1007/s00125-025-06499-z)

4. **Front Endocrinol 2021.** Pathological β-Cell ER Stress in Type 2 Diabetes. [PMC8151793](https://pmc.ncbi.nlm.nih.gov/articles/PMC8151793/)

5. **IJMS 2022.** ER Stress and Pancreatic β-Cell Dysfunction and Senescence. [PMC9104816](https://pmc.ncbi.nlm.nih.gov/articles/PMC9104816/)

6. **Diabetes ADA 2021.** Living Dangerously: ER Stress Responses in β-Cells. [Link](https://diabetesjournals.org/diabetes/article/70/11/2431/123866/)

7. **Nature EMM 2023.** Reversing pancreatic β-cell dedifferentiation in T2D. [Link](https://www.nature.com/articles/s12276-023-01043-8)

8. **Nutrients 2021.** Brief Review of β-Cell Dedifferentiation Mechanisms. [PMC8151793](https://pmc.ncbi.nlm.nih.gov/articles/PMC8151793/)

9. **Cell Metab 2014.** Pancreatic β Cell Dedifferentiation and Redifferentiation. PMID: 24703696

10. **iScience 2022.** Single-cell RNA-seq of human and mouse islets. [PMC9626680](https://pmc.ncbi.nlm.nih.gov/articles/PMC9626680/)

11. **Nat Metab 2023.** Delineating mouse β-cell identity during lifetime and diabetes. [Link](https://www.nature.com/articles/s42255-023-00876-x)

12. **PLOS Genet 2013.** Nkx6.1 Controls Gene Network for β-Cell Identity. [Link](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003274)

13. **Nat Commun 2019.** Diabetes causes marked inhibition of mitochondrial metabolism. PMID: 31171772

14. **Cells 2025.** β-Cell Mitochondrial Dysfunction: Mechanisms and Strategies. [PMC12691418](https://pmc.ncbi.nlm.nih.gov/articles/PMC12691418/)

15. **WJD 2023.** β-cell dysfunction: inflammation and oxidative stress. [PMC10075035](https://pmc.ncbi.nlm.nih.gov/articles/PMC10075035/)

16. **JCI 2019.** β Cell dysfunction during progression to T2D. [Link](https://www.jci.org/articles/view/129188)

17. **Postgrad Med 2020.** Beta-cell failure: mechanisms, markers, clinical implications. [Link](https://www.tandfonline.com/doi/full/10.1080/00325481.2020.1771047)

---

*Architecture Version: 2.0 (Literature-Integrated)*
*Total References: 17 primary sources*
*Gene Signatures: 60 curated genes*
*Last Updated: 2026-01-30*
