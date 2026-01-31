# Table 2: Therapeutic Candidates from Integrative Analysis

## Table 2. Drug candidates identified through LINCS connectivity mapping and target validation

### Panel A: Novel Mechanism Candidates (Tier 1 - High Priority)

| Rank | Compound/Class | Mechanism | Connectivity Score | Target Genes | Evidence Level | Development Stage |
|------|----------------|-----------|-------------------|--------------|----------------|-------------------|
| 1 | RHO-kinase inhibitors | Cytoskeletal regulation | 1.00 | ROCK1, ROCK2 | Preclinical | Phase II (other indications) |
| 2 | SRC-kinase inhibitors | Focal adhesion signaling | 1.00 | SRC, YES1 | Preclinical | Approved (oncology) |
| 3 | VEGF-receptor inhibitors | Vascular signaling | 1.00 | VEGFR1-3 | Preclinical | Approved (oncology) |
| 4 | Phenformin | Mitochondrial complex I | 0.95 | AMPK pathway | Historical | Withdrawn (lactic acidosis) |
| 5 | Aurora kinase inhibitors | Cell cycle regulation | 0.92 | AURKA, AURKB | Preclinical | Phase II (oncology) |

### Panel B: Known Adjacent Mechanism Candidates (Tier 2 - Medium Priority)

| Rank | Compound/Class | Mechanism | Connectivity Score | Target Genes | Evidence Level | Development Stage |
|------|----------------|-----------|-------------------|--------------|----------------|-------------------|
| 6 | HDAC inhibitors | Epigenetic modulation | 0.87 | HDAC1-6 | Strong | Approved (oncology) |
| 7 | HSP90 inhibitors | Protein chaperone | 0.82 | HSP90AA1 | Moderate | Phase III (oncology) |
| 8 | PI3K inhibitors | Insulin signaling | 0.78 | PIK3CA/CB | Moderate | Approved (oncology) |
| 9 | JAK inhibitors | Inflammatory signaling | 0.75 | JAK1-3 | Moderate | Approved (autoimmune) |
| 10 | mTOR inhibitors | Nutrient sensing | 0.71 | MTOR | Strong | Approved (transplant) |

### Panel C: Benchmark Validation (Tier 3 - Existing T2D Drugs)

| Rank | Drug | Mechanism | Connectivity Score | Primary Target | Clinical Status |
|------|------|-----------|-------------------|----------------|-----------------|
| 11 | Metformin | AMPK activation | 0.68 | AMPK/Complex I | First-line T2D |
| 12 | Pioglitazone | PPARγ agonist | 0.65 | PPARG | Second-line T2D |
| 13 | Glipizide | Sulfonylurea | 0.62 | ABCC8/KCNJ11 | Second-line T2D |
| 14 | Sitagliptin | DPP4 inhibitor | 0.58 | DPP4 | Second-line T2D |
| 15 | Glimepiride | Sulfonylurea | 0.55 | ABCC8/KCNJ11 | Second-line T2D |
| 16 | Repaglinide | Meglitinide | 0.52 | ABCC8/KCNJ11 | Second-line T2D |
| 17 | Acarbose | Alpha-glucosidase inhibitor | 0.48 | Intestinal enzymes | Adjunct therapy |

### Panel D: Novel DisGeNET-Derived Targets

| Rank | Gene | DisGeNET Score | Function | Druggability | Existing Compounds |
|------|------|----------------|----------|--------------|-------------------|
| 1 | IRS1 | 0.907 | Insulin receptor substrate | High | Research compounds |
| 2 | HNF4A | 0.725 | Transcription factor | Moderate | No approved drugs |
| 3 | AKT2 | 0.721 | Kinase signaling | High | MK-2206 (Phase II) |
| 4 | HNF1B | 0.692 | Transcription factor | Low | No approved drugs |
| 5 | TCF7L2 | 0.648 | Wnt signaling TF | Low | Research compounds |
| 6 | MAPK8IP1 | 0.620 | JNK scaffold protein | Moderate | Research compounds |
| 7 | CAPN10 | 0.616 | Calpain protease | Moderate | Calpain inhibitors |
| 8 | IL6 | 0.593 | Inflammatory cytokine | High | Tocilizumab (approved) |
| 9 | INSR | 0.590 | Insulin receptor | High | Insulin (approved) |
| 10 | SLC2A4 | 0.551 | GLUT4 transporter | Low | No direct modulators |

---

## Table Legend

**Table 2. Drug candidates identified through LINCS connectivity mapping and target validation.**

**(A)** Tier 1 novel mechanism candidates showing highest connectivity scores (ability to reverse workload stress gene signatures) in LINCS L1000 database. These compounds operate through mechanisms not currently targeted in T2D therapy and represent high-priority candidates for further investigation. RHO-kinase inhibitors showed preclinical beta-cell protective effects. SRC-kinase inhibitors modulate glucose-stimulated insulin secretion pathways.

**(B)** Tier 2 candidates with known adjacent mechanisms that could be repurposed for T2D. HDAC inhibitors show potential for PDX1/MAFA reactivation. HSP90 and mTOR inhibitors may modulate protein homeostasis and nutrient sensing relevant to beta-cell stress.

**(C)** Tier 3 benchmark validation showing recovery of existing T2D drugs in connectivity analysis, validating our computational approach. Detection of metformin, sulfonylureas, and thiazolidinediones confirms biological relevance of identified signatures.

**(D)** Novel therapeutic targets identified from DisGeNET gene-disease associations among workload pathway components. DisGeNET scores represent integrated evidence from GWAS, animal models, and literature. Druggability assessed based on protein class, structural information, and existing compound availability.

---

## Column Definitions

- **Connectivity Score**: LINCS L1000 connectivity score measuring signature reversal (-1 to +1 scale; higher positive values indicate stronger reversal of stress signature)
- **Evidence Level**: Strength of preclinical/mechanistic evidence for T2D relevance
- **Development Stage**: Current drug development status for any indication
- **DisGeNET Score**: Gene-disease association score (0-1 scale; higher = stronger T2D association)
- **Druggability**: Assessment of target tractability based on protein class and existing modulators

## Methods Summary

LINCS L1000 connectivity mapping queried workload stress signatures (genes differentially expressed at CWI > 0.6 vs CWI < 0.3, n=847 genes) against 473,647 perturbation signatures across 27,927 compounds. Connectivity scores were computed using the Connectivity Map Query algorithm. Drug-target annotations derived from ChEMBL database. DisGeNET associations from processed dataset (82,833 total associations) filtered for T2D (DOID:9352).
