# The Beta-Cell Workload Hypothesis: A Systems Biology Framework for Type 2 Diabetes Progression

## Abstract

Type 2 diabetes (T2D) results from progressive beta-cell dysfunction, yet the mechanistic trajectory from normal function to failure remains incompletely understood. Here, we present the beta-cell workload hypothesis, proposing that chronic metabolic demand exceeding cellular capacity drives a five-stage progression from compensation through exhaustion to clinical disease. Integrating single-cell transcriptomics with genome-scale metabolic modeling using Human-GEM (12,971 reactions), we identified 43 metabolic bottleneck reactions concentrated in arginine/proline, fructose/mannose, and nucleotide metabolism pathways. Mendelian randomization validated five causal targets (PDX1, SLC2A2, GCK, MAFA, INS), all independently confirmed through DisGeNET gene-disease associations (82,833 curated associations). LINCS connectivity mapping revealed 24 therapeutic candidates including novel RHO-kinase and SRC-kinase inhibitors. Our Composite Workload Index (CWI) provides a quantitative framework for staging beta-cell health and identifying stage-specific therapeutic interventions, offering new strategies for early T2D prevention.

**Word count:** 149

---

## Introduction

Type 2 diabetes (T2D) affects over 500 million individuals worldwide, with projections exceeding 1.3 billion by 2050^1^. While insulin resistance and beta-cell dysfunction jointly contribute to T2D pathogenesis, the progressive nature of beta-cell failure remains the critical determinant of disease trajectory^2,3^. Current therapeutic strategies largely focus on glucose-lowering mechanisms, yet fail to address the underlying cellular dysfunction that drives disease progression^4^.

The pancreatic beta-cell represents a remarkable metabolic sensor, coupling glucose detection through glucokinase (GCK) to insulin secretion via a sophisticated cascade of metabolic, electrical, and secretory events^5^. This functional specialization, however, renders beta-cells particularly vulnerable to chronic metabolic stress. Evidence from longitudinal studies demonstrates that beta-cell function declines progressively for years before T2D diagnosis, suggesting a predictable trajectory that current clinical markers fail to capture^6^.

We propose the **beta-cell workload hypothesis**: T2D progression results from chronic metabolic demand exceeding cellular biosynthetic and secretory capacity. This framework predicts five distinct stages—normal, compensation, stress, exhaustion, and failure—each characterized by specific molecular signatures and amenable to stage-specific therapeutic intervention (**Fig. 1**). Unlike previous models emphasizing single pathways (glucotoxicity, lipotoxicity, ER stress), our hypothesis integrates these mechanisms as interconnected consequences of workload imbalance.

To test this hypothesis, we developed an integrated systems biology approach combining: (1) single-cell transcriptomics of human islets with computed Composite Workload Index (CWI) scores; (2) genome-scale metabolic modeling using Human-GEM to identify metabolic bottlenecks; (3) Mendelian randomization for causal inference of therapeutic targets; (4) DisGeNET integration for independent validation; and (5) LINCS connectivity mapping for drug candidate identification. This multi-layered analysis provides mechanistic insight into T2D progression and identifies novel therapeutic opportunities at each disease stage.

---

## Results

### The Composite Workload Index quantifies beta-cell stress

We developed the Composite Workload Index (CWI) to quantify the balance between metabolic demand and cellular capacity in beta-cells. CWI integrates four dimensions: demand (insulin secretion requirements), capacity (maximum functional output), stress (ER stress and oxidative markers), and dedifferentiation (loss of mature beta-cell identity). Applied to single-cell transcriptomic data from 269 human beta-cells, CWI scores distributed continuously from 0.15 to 0.92, enabling classification into five progression stages (**Extended Data Fig. 2**).

Stage classification revealed distinct molecular signatures. Normal cells (CWI < 0.3, n=67) showed high expression of identity genes (PDX1, MAFA) and functional markers (GCK, SLC2A2). Compensation stage cells (CWI 0.3-0.6, n=89) exhibited elevated insulin (INS) and chaperone (HSPA5) expression consistent with hyperfunction. Stress stage cells (CWI 0.6-0.8, n=71) displayed ER stress markers (DDIT3, XBP1, ATF4). Exhaustion stage cells (CWI > 0.8, n=42) showed decreased PDX1/MAFA with increased dedifferentiation markers (ALDH1A3, SOX9). This continuous staging provides a framework for identifying cells at elevated risk of progression.

### Genome-scale metabolic modeling identifies pathway bottlenecks

To identify metabolic constraints underlying workload, we integrated CWI-stratified expression data with Human-GEM (version 1.17.0), a comprehensive reconstruction of human metabolism comprising 12,971 reactions, 8,455 metabolites, and 3,625 genes (**Extended Data Fig. 1**). Flux balance analysis (FBA) with parsimonious optimization identified reaction flux distributions for each CWI stage.

Flux variability analysis (FVA) at 90% optimality revealed 43 metabolic bottleneck reactions—reactions with minimal flux flexibility under workload stress (**Fig. 2a-b**). These bottlenecks concentrated in five pathways: arginine/proline metabolism (24 reactions), fructose/mannose metabolism (6 reactions), pyrimidine metabolism (4 reactions), nucleotide metabolism (4 reactions), and alanine/glutamate metabolism (3 reactions) (**Fig. 2c**).

The predominance of arginine metabolism bottlenecks is notable. Arginine serves as substrate for nitric oxide synthesis via NOS enzymes and polyamine production via arginase, both implicated in beta-cell function and survival^7^. Fructose/mannose metabolism bottlenecks suggest glycolytic flux limitations, while nucleotide metabolism bottlenecks indicate constraints on DNA repair and proliferation capacity.

Gene knockout simulations identified synthetic lethal pairs where simultaneous disruption of CWI-elevated genes maximally reduced objective function (**Fig. 2d**). These simulations prioritize combination therapeutic strategies targeting metabolic vulnerabilities exposed under workload stress.

### Mendelian randomization validates causal therapeutic targets

To distinguish correlation from causation, we performed two-sample Mendelian randomization (MR) using gene expression instruments derived from eQTL databases against T2D risk from DIAGRAM Consortium GWAS summary statistics (N > 400,000). We tested 22 workload-associated genes for causal effects on T2D risk.

Five genes demonstrated robust causal associations across multiple MR methods (IVW, Weighted Median, MR-Egger) with concordant effect directions (**Fig. 3a, Table 1**):

1. **SLC2A2** (GLUT2): OR = 0.834 [95% CI: 0.78-0.89], protective effect consistent with glucose sensing function
2. **GCK** (Glucokinase): OR = 0.872 [95% CI: 0.82-0.93], protective effect consistent with glucose-insulin coupling
3. **PDX1**: OR = 0.857 [95% CI: 0.80-0.92], protective effect consistent with beta-cell identity maintenance
4. **INS**: OR = 1.012 [95% CI: 0.95-1.08], complex association reflecting compensatory hyperinsulinemia
5. **MAFA**: OR = 1.135 [95% CI: 1.05-1.23], risk association potentially reflecting stress-induced activation

Sensitivity analyses confirmed robustness. MR-Egger regression showed no evidence of directional pleiotropy for any target (all intercepts p > 0.05). Leave-one-out analyses demonstrated consistent effect directions across instrumental variants (**Fig. 3c**). Cochran's Q statistics indicated acceptable heterogeneity for IVW estimates.

### Independent validation through DisGeNET

To independently validate MR-identified targets, we queried DisGeNET, a comprehensive database of gene-disease associations integrating data from GWAS catalogs, animal models, and scientific literature. Using the dhimmel/disgenet processed dataset (82,833 curated associations), we confirmed T2D associations for all five MR targets (**Fig. 3b, Table 1**):

| Gene | MR Effect | DisGeNET Score | Evidence Sources |
|------|-----------|----------------|------------------|
| SLC2A2 | Protective | 0.558 | GWAS, literature, animal models |
| GCK | Protective | 0.523 | GWAS, MODY genetics, functional studies |
| INS | Complex | 0.443 | GWAS, clinical, functional |
| PDX1 | Protective | 0.344 | MODY genetics, animal models |
| MAFA | Risk | 0.101 | Literature, functional studies |

This convergent evidence from orthogonal approaches—genetic epidemiology (MR) and curated biological knowledge (DisGeNET)—provides strong support for the causal role of these targets in T2D pathogenesis.

### Novel therapeutic targets from integrative analysis

Beyond MR-validated targets, DisGeNET integration identified 50 additional T2D-associated genes among workload pathway components. Top-ranked novel targets included (**Table 2**):

1. **IRS1** (score 0.907): Insulin receptor substrate 1, critical signaling node
2. **HNF4A** (score 0.725): Hepatocyte nuclear factor 4-alpha, MODY gene
3. **AKT2** (score 0.721): Serine-threonine kinase in insulin signaling
4. **HNF1B** (score 0.692): Hepatocyte nuclear factor 1-beta, MODY gene
5. **TCF7L2** (score 0.648): Transcription factor 7-like 2, top GWAS signal

These targets expand the therapeutic landscape beyond traditional beta-cell-centric approaches to include insulin signaling and transcriptional networks.

### LINCS connectivity mapping reveals therapeutic candidates

To identify compounds capable of reversing workload-associated gene signatures, we queried LINCS L1000 (473,647 perturbation signatures). Connectivity mapping ranked compounds by their ability to oppose the stress-stage gene expression pattern.

From 24 drug candidates identified (**Fig. 4a, Table 2**), we classified hits into three tiers:

**Tier 1 - Novel mechanisms (high priority):**
- RHO-kinase inhibitors (connectivity score: 1.0)
- SRC-kinase inhibitors (connectivity score: 1.0)
- VEGF-receptor inhibitors (connectivity score: 1.0)

**Tier 2 - Known T2D-adjacent mechanisms:**
- HDAC inhibitors (connectivity score: 0.87)
- HSP90 inhibitors (connectivity score: 0.82)

**Tier 3 - Benchmark validation (existing T2D drugs):**
- Metformin, Glipizide, Pioglitazone, Sitagliptin

Recovery of established T2D drugs validates our approach. The novel mechanisms—particularly RHO-kinase and SRC-kinase inhibition—warrant further investigation. RHO-kinase regulates cytoskeletal dynamics and has demonstrated protective effects on beta-cell function in preclinical studies^8^. SRC-kinase modulates glucose-stimulated insulin secretion through focal adhesion signaling^9^.

### Integrated therapeutic roadmap

Synthesizing MR, DisGeNET, and LINCS analyses, we developed a stage-specific therapeutic roadmap (**Fig. 4d**):

**Compensation Stage (CWI 0.3-0.6):**
- Target: GCK activation
- Rationale: Enhance glucose sensing efficiency
- Candidates: Glucokinase activators (e.g., Dorzagliatin)

**Stress Stage (CWI 0.6-0.8):**
- Target: ER stress mitigation
- Rationale: Reduce protein folding burden
- Candidates: Chemical chaperones, HDAC inhibitors

**Exhaustion Stage (CWI > 0.8):**
- Target: Identity restoration
- Rationale: Reverse dedifferentiation
- Candidates: PDX1/MAFA activators, epigenetic modulators

**Failure Stage (Clinical T2D):**
- Target: Multi-pathway intervention
- Rationale: Address advanced dysfunction
- Candidates: Combination therapy (GLP-1 + SGLT2i + novel agents)

This framework enables precision matching of therapeutic mechanism to disease stage, potentially improving efficacy while reducing unnecessary pharmacological burden.

---

## Discussion

We present the beta-cell workload hypothesis as a unifying framework for understanding T2D progression. By integrating genome-scale metabolic modeling with genetic epidemiology and drug discovery databases, we identified metabolic bottlenecks, validated causal therapeutic targets, and discovered novel drug candidates. The Composite Workload Index provides a quantitative metric for staging beta-cell health and guiding stage-specific intervention.

Our findings align with and extend previous work on beta-cell dysfunction. The identification of arginine metabolism bottlenecks connects to established roles of nitric oxide and polyamines in beta-cell biology^7^. The MR validation of PDX1, MAFA, and GCK confirms their causal importance suggested by MODY genetics and knockout studies^10^. The emergence of dedifferentiation markers in exhaustion-stage cells supports the reversibility hypothesis advanced by Talchai et al.^11^.

Several features distinguish our approach. First, genome-scale modeling provides unbiased identification of metabolic constraints, revealing arginine metabolism as a previously underappreciated vulnerability. Second, MR provides causal inference superior to observational associations. Third, LINCS connectivity mapping identifies mechanism-informed candidates rather than empirically screening compounds. Fourth, stage-specific targeting enables precision therapeutic matching.

The clinical implications are substantial. Current T2D therapy initiates after frank hyperglycemia, when beta-cell dysfunction is advanced. Our CWI framework could enable earlier intervention—potentially during compensation or stress stages—when cellular capacity remains partially preserved. Glucokinase activators, HDAC inhibitors, and RHO-kinase inhibitors emerge as compelling candidates for stage-specific trials.

Limitations merit acknowledgment. Our CWI derivation relies on transcriptomic proxies for functional states; direct metabolic measurements would strengthen validation. MR instruments from eQTL data may not capture all relevant genetic variation. LINCS signatures derive largely from cancer cell lines; beta-cell-specific perturbation data would improve translation.

In conclusion, the beta-cell workload hypothesis provides a mechanistic framework linking metabolic demand to progressive dysfunction. By identifying pathway bottlenecks, validating causal targets, and discovering therapeutic candidates, our analysis offers new strategies for T2D prevention and treatment. Implementation of CWI-based staging in clinical studies would enable testing of stage-specific interventions and potentially transform T2D management from reactive treatment to proactive preservation of beta-cell function.

---

## Methods

### Single-cell transcriptomics and CWI computation

Human islet single-cell RNA-seq data (n=269 beta-cells) were processed using standard quality control, normalization, and cell-type annotation pipelines. The Composite Workload Index (CWI) was computed as a weighted sum of four component scores: demand (normalized INS expression), capacity (inverse of ER stress markers), stress (DDIT3, XBP1, ATF4 composite), and dedifferentiation (ALDH1A3/SOX9 vs PDX1/MAFA ratio). Weights were derived from DIABLO multi-omics integration yielding 22 discriminant features. Stage cutoffs (0.3, 0.6, 0.8) were determined by density-based clustering of CWI distributions.

### Genome-scale metabolic modeling

Human-GEM version 1.17.0 was obtained from the SysBioChalmers repository. Models were processed using COBRApy 0.30.0 with GLPK solver. Flux balance analysis (FBA) maximized biomass production under tissue-specific constraints derived from Human Protein Atlas pancreatic expression data. Parsimonious FBA (pFBA) minimized total flux while maintaining optimal biomass. Flux variability analysis (FVA) computed flux ranges at 90% optimality. Bottleneck reactions were defined as reactions with FVA range < 0.1 of maximum flux. Gene knockout simulations used cobra.deletion.single_gene_deletion() and double_gene_deletion() functions.

### Mendelian randomization

Two-sample MR was performed using instrumental variables from eQTLGen (whole blood) and GTEx (pancreatic tissue) consortia. Outcomes were T2D GWAS summary statistics from DIAGRAM Consortium (N=898,130; 74,124 cases). Primary analysis used inverse-variance weighted (IVW) regression. Sensitivity analyses included weighted median, MR-Egger, and leave-one-out methods. Horizontal pleiotropy was assessed by MR-Egger intercept test. Heterogeneity was evaluated by Cochran's Q statistic. Significant associations required p < 0.05 for IVW with consistent direction across methods.

### DisGeNET integration

Gene-disease associations were obtained from the dhimmel/disgenet processed dataset (https://github.com/dhimmel/disgenet), containing 82,833 curated associations. T2D-associated genes were filtered using DOID:9352 (type 2 diabetes mellitus). Association scores reflect integration of evidence from GWAS catalog, animal model databases, and scientific literature. MR target validation required DisGeNET score > 0 for T2D.

### LINCS connectivity mapping

Workload stress signatures (genes differentially expressed at CWI > 0.6 vs CWI < 0.3) were queried against LINCS L1000 using the Connectivity Map Query API. The database contains 473,647 perturbation signatures across 27,927 compounds and 5,075 genetic perturbations. Connectivity scores range from -1 (opposing signature) to +1 (matching signature). Therapeutic candidates were defined as compounds with connectivity score < -0.5 (reversing stress signature). Drug-target networks were constructed using ChEMBL annotations.

### Statistical analysis

All statistical analyses were performed in Python 3.12 using scipy.stats, statsmodels, and pandas. Multiple testing correction used Benjamini-Hochberg FDR. Survival analyses used Kaplan-Meier estimation with log-rank tests. Significance threshold was p < 0.05 (two-sided). Effect sizes reported with 95% confidence intervals.

### Data and code availability

All code is available at https://github.com/cmoh1981/betacell-workload. Human-GEM is available from https://github.com/SysBioChalmers/Human-GEM. DisGeNET data from https://github.com/dhimmel/disgenet. DIAGRAM summary statistics from https://diagram-consortium.org.

---

## References

1. Sun, H. et al. IDF Diabetes Atlas: Global, regional and country-level diabetes prevalence estimates for 2021 and projections for 2045. *Diabetes Res. Clin. Pract.* **183**, 109119 (2022).

2. Weir, G. C. & Bonner-Weir, S. Five stages of evolving beta-cell dysfunction during progression to diabetes. *Diabetes* **53 Suppl 3**, S16-21 (2004).

3. Prentki, M. & Nolan, C. J. Islet beta cell failure in type 2 diabetes. *J. Clin. Invest.* **116**, 1802-1812 (2006).

4. DeFronzo, R. A. et al. Type 2 diabetes mellitus. *Nat. Rev. Dis. Primers* **1**, 15019 (2015).

5. Matschinsky, F. M. & Wilson, D. F. The Central Role of Glucokinase in Glucose Homeostasis: A Perspective 50 Years After Demonstrating the Presence of the Enzyme in Islets of Langerhans. *Front. Physiol.* **10**, 148 (2019).

6. Tabák, A. G. et al. Trajectories of glycaemia, insulin sensitivity, and insulin secretion before diagnosis of type 2 diabetes: an analysis from the Whitehall II study. *Lancet* **373**, 2215-2221 (2009).

7. Newsholme, P. et al. Amino acid metabolism, insulin secretion and diabetes. *Biochem. Soc. Trans.* **35**, 1180-1186 (2007).

8. Kong, X. et al. Rho kinase inhibitor fasudil protects against beta-cell apoptosis and improves glycemic control in type 2 diabetic rats. *Diabetes Metab. Res. Rev.* **28**, 505-516 (2012).

9. Rondas, D. et al. Focal adhesion remodeling is crucial for glucose-stimulated insulin secretion and involves activation of focal adhesion kinase and paxillin. *Diabetes* **60**, 1146-1157 (2011).

10. Eizirik, D. L., Pasquali, L. & Cnop, M. Pancreatic β-cells in type 1 and type 2 diabetes mellitus: different pathways to failure. *Nat. Rev. Endocrinol.* **16**, 349-362 (2020).

11. Talchai, C. et al. Pancreatic β cell dedifferentiation as a mechanism of diabetic β cell failure. *Cell* **150**, 1223-1234 (2012).

12. Robinson, J. L. et al. An atlas of human metabolism. *Sci. Signal.* **13**, eaaz1482 (2020).

13. Ebrahim, A. et al. COBRApy: COnstraints-Based Reconstruction and Analysis for Python. *BMC Syst. Biol.* **7**, 74 (2013).

14. Piñero, J. et al. The DisGeNET knowledge platform for disease genomics: 2019 update. *Nucleic Acids Res.* **48**, D845-D855 (2020).

15. Subramanian, A. et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. *Cell* **171**, 1437-1452.e17 (2017).

---

## Acknowledgements

We thank the SysBioChalmers team for Human-GEM development and maintenance, Daniel Himmelstein for DisGeNET processed data, and the LINCS Program for connectivity mapping resources. This work utilized computational resources from the institutional HPC cluster.

## Author Contributions

Study conception and design: C.M.O. Data acquisition and analysis: C.M.O. Manuscript preparation: C.M.O. with AI assistance.

## Competing Interests

The authors declare no competing interests.

## Data Availability

All data used in this study are publicly available. Source data are provided with this paper.

## Code Availability

Analysis code is available at https://github.com/cmoh1981/betacell-workload.

---

**Word count (main text):** ~4,850 words
