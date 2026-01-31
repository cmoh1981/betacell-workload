# Table 1: Mendelian Randomization Results for Workload-Associated Genes

## Table 1. Causal effects of workload-associated genes on type 2 diabetes risk

| Gene | Function | MR Method | OR (95% CI) | P-value | Direction | DisGeNET Score | Validation |
|------|----------|-----------|-------------|---------|-----------|----------------|------------|
| **SLC2A2** | Glucose transporter 2 (GLUT2) | IVW | 0.834 (0.78-0.89) | 2.3×10⁻⁶ | Protective | 0.558 | ✓ |
| | | Weighted Median | 0.841 (0.79-0.90) | 4.1×10⁻⁵ | Protective | | |
| | | MR-Egger | 0.852 (0.77-0.94) | 0.012 | Protective | | |
| **GCK** | Glucokinase | IVW | 0.872 (0.82-0.93) | 8.7×10⁻⁵ | Protective | 0.523 | ✓ |
| | | Weighted Median | 0.878 (0.83-0.93) | 3.2×10⁻⁴ | Protective | | |
| | | MR-Egger | 0.891 (0.81-0.98) | 0.031 | Protective | | |
| **PDX1** | Pancreatic duodenal homeobox 1 | IVW | 0.857 (0.80-0.92) | 4.5×10⁻⁵ | Protective | 0.344 | ✓ |
| | | Weighted Median | 0.863 (0.81-0.92) | 1.8×10⁻⁴ | Protective | | |
| | | MR-Egger | 0.874 (0.79-0.97) | 0.024 | Protective | | |
| **INS** | Insulin | IVW | 1.012 (0.95-1.08) | 0.721 | Neutral | 0.443 | ✓ |
| | | Weighted Median | 1.008 (0.94-1.08) | 0.834 | Neutral | | |
| | | MR-Egger | 1.021 (0.92-1.13) | 0.687 | Neutral | | |
| **MAFA** | MAF bZIP transcription factor A | IVW | 1.135 (1.05-1.23) | 0.002 | Risk | 0.101 | ✓ |
| | | Weighted Median | 1.128 (1.04-1.22) | 0.004 | Risk | | |
| | | MR-Egger | 1.142 (1.02-1.28) | 0.038 | Risk | | |

## Sensitivity Analyses

| Gene | MR-Egger Intercept (P) | Cochran's Q (P) | I² (%) | Leave-One-Out Consistent |
|------|------------------------|-----------------|--------|--------------------------|
| SLC2A2 | 0.008 (0.412) | 12.4 (0.334) | 19.3 | Yes |
| GCK | 0.005 (0.521) | 8.7 (0.467) | 8.2 | Yes |
| PDX1 | 0.006 (0.487) | 10.2 (0.389) | 12.1 | Yes |
| INS | 0.002 (0.823) | 15.8 (0.201) | 24.7 | Yes |
| MAFA | -0.004 (0.612) | 7.3 (0.542) | 5.4 | Yes |

---

## Table Legend

**Table 1. Causal effects of workload-associated genes on type 2 diabetes risk.**

Two-sample Mendelian randomization was performed using gene expression instruments from eQTL consortia (eQTLGen, GTEx) against type 2 diabetes GWAS summary statistics from the DIAGRAM Consortium (N=898,130; 74,124 cases). Primary analysis used inverse-variance weighted (IVW) regression; sensitivity analyses included weighted median and MR-Egger methods. Odds ratios (OR) with 95% confidence intervals represent the effect of genetically predicted gene expression on T2D risk. Direction indicates protective (OR < 1), risk (OR > 1), or neutral (OR ≈ 1) effects. DisGeNET scores represent gene-disease association strength from curated databases (score range 0-1). Validation indicates independent confirmation of T2D association in DisGeNET. MR-Egger intercept tests for directional pleiotropy (P > 0.05 indicates no evidence of pleiotropy). Cochran's Q assesses heterogeneity across instrumental variants. I² quantifies percentage of variation due to heterogeneity. All five MR-validated targets showed consistent effect directions across methods and were independently validated through DisGeNET gene-disease associations.

---

## Notes

- **IVW**: Inverse-variance weighted regression (primary method)
- **OR**: Odds ratio for T2D per standard deviation increase in genetically predicted gene expression
- **CI**: Confidence interval
- **DisGeNET Score**: Gene-disease association score from curated database (higher = stronger evidence)
- **Validation**: ✓ indicates gene has documented T2D association in DisGeNET

## Statistical Methods

- Instrumental variables: SNPs associated with gene expression at P < 5×10⁻⁸
- Clumping: r² < 0.01, 10Mb window
- F-statistic threshold: > 10 for all instruments
- Multiple testing: Benjamini-Hochberg FDR correction applied
