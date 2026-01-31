# Drug Discovery Pipeline Results Summary

**Date:** 2026-01-31
**Project:** Beta-Cell Workload Modulator Discovery

---

## Executive Summary

This analysis identified drug candidates that could modulate beta-cell workload based on:
1. **LINCS L1000 Connectivity Mapping** - 19,811 compounds screened
2. **Target-Based Screening** - 6 MR-validated workload targets
3. **ADMET Prediction** - Drug-likeness assessment of top candidates

### Key Findings

| Metric | Value |
|--------|-------|
| Total LINCS compounds screened | 19,811 |
| Known T2D drugs in LINCS | 11/19 |
| MR-validated targets | 6 genes |
| ChEMBL target compounds | 10 |
| Candidates with good drug-likeness | 11/11 (100%) |

---

## Module A: LINCS Connectivity Mapping

### Workload Query Signature
- **Genes to DECREASE** (stress/dedifferentiation): 15 genes
  - DDIT3, ATF4, XBP1, HSPA5, ERN1, EIF2AK3, TRIB3, ATF6, CALR, CANX, HSP90B1
  - ALDH1A3, SOX9, HES1, NEUROG3

- **Genes to INCREASE** (capacity): 12 genes
  - INS, IAPP, PCSK1, PCSK2, CPE, SCG2, CHGB
  - GCK, SLC2A2, PDX1, MAFA, NKX6-1, ABCC8, KCNJ11

### Known T2D Drugs Found in LINCS
| Drug | Found | Signatures |
|------|-------|------------|
| Metformin | Yes | 1 |
| Glipizide | Yes | 1 |
| Glimepiride | Yes | 2 |
| Pioglitazone | Yes | 1 |
| Rosiglitazone | Yes | 1 |
| Sitagliptin | Yes | 1 |
| Tolbutamide | Yes | 1 |
| Repaglinide | Yes | 1 |
| Nateglinide | Yes | 2 |
| Acarbose | Yes | 2 |
| Miglitol | Yes | 1 |

### Top LINCS Hits (Novel Candidates)
1. **Kinase Inhibitors** - Multiple RHO-kinase, SRC-kinase, VEGF-receptor inhibitors
   - Potential mechanism: Reduce cellular stress pathways

---

## Module B: Target-Based Screening

### MR-Validated Workload Targets

| Gene | UniProt | MR OR | Effect | Therapeutic Action |
|------|---------|-------|--------|-------------------|
| **PDX1** | P52945 | 0.66 | Protective | Activate/Increase |
| **SLC2A2** | P11168 | 0.89 | Protective | Activate/Increase |
| **MAFA** | Q8NHW3 | 1.14 | Risk | Modulate |
| **GCK** | P35557 | - | Key metabolic | Activate |
| **HSPA5** | P11021 | - | Stress marker | Support/Chaperone |
| **DDIT3** | P35638 | - | Stress marker | Inhibit |

### Prioritized Drug Candidates

| Rank | Compound | Target | Action | Status | Score |
|------|----------|--------|--------|--------|-------|
| 1 | **Dorzagliatin** | GCK | Activator | Approved (China) | 6.0 |
| 2 | **4-PBA** | HSPA5 | Chemical chaperone | Approved | 5.5 |
| 3 | **TUDCA** | HSPA5 | Chemical chaperone | Approved | 5.5 |
| 4 | **Glibenclamide** | ABCC8 | Inhibitor | Approved | 5.0 |
| 5 | **Diazoxide** | ABCC8 | Opener | Approved | 5.0 |
| 6 | **BRD7552** | PDX1 | Expression inducer | Research | 3.5 |
| 7 | **Piragliatin** | GCK | Activator | Phase II | 3.0 |
| 8 | **MK-0941** | GCK | Activator | Phase II | 3.0 |
| 9 | **Salubrinal** | DDIT3 | Indirect inhibitor | Research | 2.5 |
| 10 | **Azoramide** | HSPA5 | ER proteostasis | Research | 2.0 |

---

## Module C: ADMET Analysis

### Drug-Likeness Assessment (Ranked by QED Score)

| Compound | QED | MW | LogP | Lipinski | Absorption | BBB |
|----------|-----|-----|------|----------|------------|-----|
| **Piragliatin** | 0.792 | 339 | 3.10 | PASS | Good | Likely |
| **4-PBA** | 0.769 | 178 | 2.73 | PASS | Good | Likely |
| **Diazoxide** | 0.737 | 231 | 1.68 | PASS | Good | Likely |
| **Glipizide** | 0.607 | 445 | 2.68 | PASS | Good | Unlikely |
| **Dapagliflozin** | 0.582 | 409 | 1.84 | PASS | Good | Unlikely |
| **Dorzagliatin** | 0.463 | 432 | 4.97 | PASS | Good | Unlikely |
| **Empagliflozin** | 0.454 | 457 | 3.24 | PASS | Good | Unlikely |
| **TUDCA** | 0.398 | 500 | 3.40 | PASS | Good | Unlikely |
| **Glibenclamide** | 0.382 | 516 | 3.55 | PASS (1 viol) | Moderate | Unlikely |
| **Salubrinal** | 0.369 | 362 | 3.59 | PASS | Good | Likely |
| **Metformin** | 0.249 | 129 | -1.03 | PASS | Good | Unlikely |

---

## Top Recommended Candidates

### Tier 1: Approved Drugs (Immediate Repurposing Potential)

1. **Dorzagliatin** (GCK activator)
   - Status: Approved in China for T2D
   - Mechanism: Directly activates glucokinase, restores glucose sensing
   - Workload relevance: Addresses metabolic capacity (GCK pathway)
   - QED: 0.463

2. **4-PBA / TUDCA** (ER stress modulators)
   - Status: Approved for other indications
   - Mechanism: Chemical chaperones reduce ER stress
   - Workload relevance: Addresses stress component of workload
   - QED: 0.769 / 0.398

3. **Diazoxide** (KATP opener)
   - Status: Approved
   - Mechanism: Opens KATP channels, reduces insulin secretion demand
   - Workload relevance: Directly reduces beta-cell workload
   - QED: 0.737

### Tier 2: Research Compounds (Novel Mechanisms)

4. **BRD7552** (PDX1 inducer)
   - Status: Research
   - Mechanism: Increases PDX1 expression
   - Workload relevance: MR-validated protective target (OR=0.66)
   - Note: Aligns with MR findings - higher PDX1 = reduced T2D risk

5. **Salubrinal** (DDIT3/CHOP inhibitor)
   - Status: Research tool
   - Mechanism: Blocks CHOP induction via eIF2α
   - Workload relevance: Reduces pro-apoptotic stress signaling

---

## Therapeutic Strategy

Based on workload framework, optimal treatment combines:

```
┌─────────────────────────────────────────────────────────────────┐
│                 WORKLOAD REDUCTION STRATEGY                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  1. REDUCE DEMAND          2. INCREASE CAPACITY                 │
│  ────────────────          ─────────────────────                │
│  • SGLT2 inhibitors        • GCK activators                     │
│    (Dapagliflozin)           (Dorzagliatin)                     │
│  • Diazoxide               • PDX1 inducers                      │
│    (reduce secretion)        (BRD7552)                          │
│                                                                 │
│  3. REDUCE STRESS          4. PREVENT DEDIFFERENTIATION         │
│  ────────────────          ────────────────────────────         │
│  • 4-PBA, TUDCA            • Maintain metabolic identity        │
│    (ER chaperones)         • Support NKX6-1, MAFA               │
│  • Salubrinal              • Avoid chronic stress               │
│    (CHOP inhibition)                                            │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## Files Generated

| File | Description |
|------|-------------|
| `lincs_ranked_compounds.csv` | Top 100 LINCS compounds by connectivity score |
| `workload_modulator_candidates.csv` | 24 identified workload modulator candidates |
| `known_t2d_drugs_lincs.csv` | T2D drugs presence in LINCS |
| `workload_query_signature.csv` | 27-gene workload query signature |
| `chembl_target_compounds.csv` | 10 compounds targeting workload genes |
| `prioritized_target_compounds.csv` | Ranked target-based candidates |
| `workload_drug_targets.csv` | 6 MR-validated drug targets |
| `admet_analysis.csv` | Full ADMET predictions |
| `admet_report.txt` | ADMET summary report |

---

---

## Module D: Network Pharmacology

### Network Statistics
- **Protein-protein interactions**: 5,402 edges from STRING
- **Network nodes**: 2,859 proteins
- **Pathways analyzed**: 21,889 (MSigDB C2, Hallmark, GO BP)
- **Enriched pathways (FDR<0.05)**: 1,012

### Drug Network Pharmacology Scores

| Rank | Drug | NP Score | Targets | Relevant Pathways |
|------|------|----------|---------|-------------------|
| 1 | **BRD7552** | 105.1 | PDX1 | 50 |
| 2 | **Dorzagliatin** | 91.4 | GCK | 43 |
| 3 | **Piragliatin** | 91.4 | GCK | 43 |
| 4 | **Salubrinal** | 90.9 | EIF2AK3, DDIT3 | 41 |
| 5 | **4-PBA** | 82.4 | HSPA5 | 38 |
| 6 | **TUDCA** | 82.4 | HSPA5 | 38 |
| 7 | **Glibenclamide** | 69.3 | ABCC8 | 32 |
| 8 | **Diazoxide** | 69.3 | ABCC8 | 32 |

### Top Synergistic Drug Combinations

| Rank | Combination | Synergy Score | Mechanism |
|------|-------------|---------------|-----------|
| 1 | **Salubrinal + Diazoxide** | 114.2 | Stress + Capacity |
| 2 | **Salubrinal + Glibenclamide** | 114.2 | Stress + Capacity |
| 3 | **Dorzagliatin + Salubrinal** | 103.8 | Capacity + Stress |
| 4 | **BRD7552 + Salubrinal** | 99.1 | Capacity + Stress |
| 5 | **Glibenclamide + Metformin** | 95.5 | Both Approved |

### Top Enriched Pathways (Workload Genes)

| Pathway | Genes | FDR |
|---------|-------|-----|
| HP_REDUCED_PANCREATIC_BETA_CELLS | 6 | 2.8e-15 |
| REACTOME_REGULATION_OF_GENE_EXPRESSION_IN_BETA_CELLS | 7 | 1.3e-14 |
| KEGG_MATURITY_ONSET_DIABETES_OF_THE_YOUNG | 7 | 3.5e-14 |
| HALLMARK_PANCREAS_BETA_CELLS | 7 | 6.7e-13 |
| KEGG_TYPE_II_DIABETES_MELLITUS | 7 | 1.7e-12 |

---

## Integrated Drug Discovery Results

### Final Prioritized Candidates

Based on all four modules (LINCS, Target-based, ADMET, Network Pharmacology):

| Drug | Target | ADMET QED | NP Score | Status | Recommendation |
|------|--------|-----------|----------|--------|----------------|
| **BRD7552** | PDX1 | - | 105.1 | Research | TOP (MR-validated) |
| **Dorzagliatin** | GCK | 0.463 | 91.4 | Approved | HIGH (clinical) |
| **4-PBA** | HSPA5 | 0.769 | 82.4 | Approved | HIGH (repurpose) |
| **TUDCA** | HSPA5 | 0.398 | 82.4 | Approved | HIGH (repurpose) |
| **Salubrinal** | DDIT3 | 0.369 | 90.9 | Research | MEDIUM (combo) |
| **Diazoxide** | ABCC8 | 0.737 | 69.3 | Approved | MEDIUM (workload) |

### Recommended Combination Strategies

```
┌──────────────────────────────────────────────────────────────────┐
│           OPTIMAL WORKLOAD MODULATION STRATEGY                   │
├──────────────────────────────────────────────────────────────────┤
│                                                                  │
│  COMBINATION 1: Maximum Synergy (Research Setting)               │
│  ─────────────────────────────────────────────────               │
│  Salubrinal (stress) + Diazoxide (capacity)                      │
│  Synergy Score: 114.2                                            │
│                                                                  │
│  COMBINATION 2: Clinically Feasible (Both Approved)              │
│  ─────────────────────────────────────────────────               │
│  Dorzagliatin (GCK) + 4-PBA (ER stress)                          │
│  Capacity enhancement + Stress reduction                         │
│                                                                  │
│  COMBINATION 3: MR-Guided (Research Priority)                    │
│  ─────────────────────────────────────────────                   │
│  BRD7552 (PDX1↑) + TUDCA (HSPA5)                                 │
│  Directly targets MR-validated gene (OR=0.66)                    │
│                                                                  │
└──────────────────────────────────────────────────────────────────┘
```

---

## Files Generated

| File | Description |
|------|-------------|
| `lincs_ranked_compounds.csv` | Top 100 LINCS compounds by connectivity score |
| `workload_modulator_candidates.csv` | 24 identified workload modulator candidates |
| `known_t2d_drugs_lincs.csv` | T2D drugs presence in LINCS |
| `workload_query_signature.csv` | 27-gene workload query signature |
| `chembl_target_compounds.csv` | 10 compounds targeting workload genes |
| `prioritized_target_compounds.csv` | Ranked target-based candidates |
| `workload_drug_targets.csv` | 6 MR-validated drug targets |
| `admet_analysis.csv` | Full ADMET predictions |
| `admet_report.txt` | ADMET summary report |
| `network_drug_effects.csv` | Network pharmacology drug scores |
| `synergistic_combinations.csv` | 55 drug combination analyses |
| `workload_pathway_enrichment.csv` | Pathway enrichment results |
| `network_pharmacology_report.txt` | Full NP analysis report |

---

## Next Steps

1. **Experimental Validation**
   - Test top candidates in beta-cell lines (MIN6, INS-1)
   - Measure workload signature genes before/after treatment
   - Priority: BRD7552 (PDX1), Dorzagliatin (GCK)

2. **Combination Studies**
   - Test synergistic pairs identified by network analysis
   - Example: Salubrinal + Diazoxide (Score: 114.2)
   - Clinically feasible: Dorzagliatin + 4-PBA

3. **Clinical Data Mining**
   - Analyze existing T2D clinical trials for workload biomarkers
   - Check if patients on 4-PBA/TUDCA show improved beta-cell function

4. **MR-Guided Drug Development**
   - BRD7552 directly increases PDX1 (MR OR=0.66, protective)
   - Strong genetic support for therapeutic mechanism

---

*Generated by Beta-Cell Workload Drug Discovery Pipeline*
*Modules: LINCS Connectivity | Target-Based Screening | ADMET | Network Pharmacology*

