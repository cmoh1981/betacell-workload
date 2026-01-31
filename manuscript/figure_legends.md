# Figure Legends

## Main Figures

### Figure 1. The beta-cell workload hypothesis: a five-stage model of diabetes progression

**(a)** Schematic representation of the beta-cell workload hypothesis showing progression from normal function through compensation, stress, exhaustion, to clinical failure. Each stage is characterized by specific Composite Workload Index (CWI) ranges and molecular markers. Arrows indicate progressive transitions driven by chronic metabolic demand exceeding cellular capacity.

**(b)** Radar plot displaying the four components of the Composite Workload Index across disease stages. Demand (insulin secretion requirements), Capacity (maximum functional output), Stress (ER stress markers), and Dedifferentiation (identity gene loss) are normalized to 0-1 scale. Normal stage shows balanced low demand with high capacity; progression involves increasing demand/stress with decreasing capacity.

**(c)** Heatmap of stage-specific molecular markers across the five progression stages. Color intensity represents normalized expression levels. Key markers include: identity genes (PDX1, MAFA, NKX6.1), functional genes (GCK, SLC2A2, INS), stress response genes (HSPA5, DDIT3, XBP1, ATF4), and dedifferentiation markers (ALDH1A3, SOX9, NANOG). Hierarchical clustering reveals distinct molecular signatures for each stage.

---

### Figure 2. Genome-scale metabolic modeling reveals pathway bottlenecks

**(a)** Distribution of metabolic flux across Human-GEM reactions under normal versus high-workload conditions. Flux balance analysis (FBA) with parsimonious optimization shows altered flux distribution in compensation and stress stages. Y-axis: flux magnitude (mmol/gDW/h); X-axis: reaction rank by flux.

**(b)** Flux variability analysis (FVA) identifying metabolic bottleneck reactions. Heatmap shows flux ranges at 90% optimality for top constrained reactions across CWI stages. Blue indicates narrow flux ranges (bottlenecks); red indicates flexible ranges. 43 reactions identified as significant bottlenecks (FVA range < 0.1 of maximum).

**(c)** Pathway enrichment of metabolic bottlenecks. Bar chart showing distribution of 43 bottleneck reactions across metabolic pathways. Arginine/proline metabolism (24 reactions), fructose/mannose metabolism (6), pyrimidine metabolism (4), nucleotide metabolism (4), and alanine/glutamate metabolism (3) represent major constrained pathways.

**(d)** Gene knockout simulation results. Scatter plot showing growth rate impact of single gene deletions for workload-associated genes. X-axis: gene knockout; Y-axis: objective function value (% of wild-type). Synthetic lethal combinations identified through double knockout analysis are highlighted.

---

### Figure 3. Mendelian randomization validates causal therapeutic targets

**(a)** Forest plot of Mendelian randomization results for five validated targets. Odds ratios with 95% confidence intervals shown for inverse-variance weighted (IVW), weighted median, and MR-Egger methods. Protective effects (OR < 1) observed for SLC2A2, GCK, and PDX1. MAFA shows risk association (OR > 1). INS shows neutral effect.

**(b)** DisGeNET validation of MR targets. Bar chart displaying gene-disease association scores from DisGeNET database for the five MR-validated genes. All targets confirmed with T2D associations (score > 0). Scores reflect integrated evidence from GWAS catalogs, animal models, and scientific literature.

**(c)** MR sensitivity analysis summary. Leave-one-out analysis for SLC2A2 demonstrating robustness of causal estimate. Each point represents IVW estimate after removing one instrumental variant. Consistency of estimates across analyses indicates reliable causal inference without single-SNP driving effects.

**(d)** Scatter plot of SNP-exposure vs SNP-outcome associations for PDX1 MR analysis. Each point represents an instrumental variant. Slope corresponds to causal estimate. IVW (blue), weighted median (green), and MR-Egger (red) regression lines shown. Concordant slopes across methods support causal interpretation.

---

### Figure 4. Therapeutic roadmap from integrative analysis

**(a)** LINCS connectivity mapping results. Bar chart showing top 24 drug candidates ranked by connectivity score (ability to reverse workload stress signature). Tier 1 (novel mechanisms, red), Tier 2 (known adjacent mechanisms, orange), and Tier 3 (existing T2D drugs, green) candidates indicated. Recovery of metformin, pioglitazone, and sulfonylureas validates approach.

**(b)** Drug-target network visualization. Network graph showing connections between identified drug candidates and their molecular targets. Node size proportional to number of drug connections. Edge thickness indicates binding affinity. Central hubs include kinase targets (SRC, ROCK, VEGFR) and metabolic regulators (AMPK, mTOR).

**(c)** Target prioritization matrix. Bubble plot integrating MR evidence (x-axis), DisGeNET score (y-axis), and druggability (bubble size) for therapeutic targets. Quadrant analysis identifies high-priority targets with strong causal evidence and existing drug candidates. GCK, PDX1, and SLC2A2 emerge as top-tier targets.

**(d)** Stage-specific therapeutic intervention roadmap. Schematic showing recommended interventions for each CWI-defined disease stage. Compensation stage: glucokinase activators targeting glucose sensing. Stress stage: chemical chaperones and HDAC inhibitors for ER stress mitigation. Exhaustion stage: epigenetic modulators for identity restoration. Failure stage: combination therapy approaches.

---

## Extended Data Figures

### Extended Data Figure 1. Human-GEM metabolic model characteristics

**(a)** Overview statistics of Human-GEM version 1.17.0. The model comprises 12,971 reactions, 8,455 metabolites, and 3,625 genes organized across 145 subsystems representing human metabolism.

**(b)** Subsystem distribution pie chart showing major metabolic categories. Exchange reactions (24%), transport (18%), amino acid metabolism (12%), lipid metabolism (11%), and carbohydrate metabolism (9%) represent largest categories.

**(c)** Gene coverage across CWI-associated pathways. Bar chart showing fraction of pathway genes with expression data in single-cell dataset. High coverage (>80%) for glycolysis, TCA cycle, and oxidative phosphorylation enables robust flux predictions.

**(d)** Comparison of flux distributions between beta-cell-constrained and generic Human-GEM models. Tissue-specific constraints derived from Human Protein Atlas pancreatic expression data significantly alter predicted flux patterns.

---

### Extended Data Figure 2. Composite Workload Index distribution and validation

**(a)** Histogram of CWI scores across 269 beta-cells. Continuous distribution from 0.15 to 0.92 with density-based clustering identifying stage boundaries at 0.3, 0.6, and 0.8.

**(b)** UMAP visualization of single-cell transcriptomes colored by CWI score. Continuous gradient from low (blue) to high (red) workload shows trajectory from normal through exhaustion stages.

**(c)** CWI component correlation matrix. Heatmap showing Pearson correlations between demand, capacity, stress, and dedifferentiation components. Demand and stress positively correlated (r=0.72); capacity negatively correlated with stress (r=-0.65).

**(d)** CWI validation against clinical markers. Box plots showing CWI scores stratified by donor HbA1c categories (normal <5.7%, prediabetes 5.7-6.4%, diabetes >6.5%). Significant increase in CWI with worsening glycemic status (ANOVA p < 0.001).

**(e)** Stage transition analysis. Sankey diagram showing cell distribution across CWI stages in donors with normal glucose tolerance versus impaired glucose tolerance. Shift toward higher CWI stages in IGT donors supports progressive model.

---

## Figure Specifications

| Figure | Dimensions | Resolution | Format |
|--------|------------|------------|--------|
| Figure 1 | 180 mm × 150 mm | 300 DPI | PDF/PNG |
| Figure 2 | 180 mm × 140 mm | 300 DPI | PDF/PNG |
| Figure 3 | 180 mm × 160 mm | 300 DPI | PDF/PNG |
| Figure 4 | 180 mm × 150 mm | 300 DPI | PDF/PNG |
| ED Figure 1 | 180 mm × 120 mm | 300 DPI | PDF/PNG |
| ED Figure 2 | 180 mm × 180 mm | 300 DPI | PDF/PNG |

All figures generated using Python matplotlib with publication-quality settings. Source data provided in Supplementary Tables.
