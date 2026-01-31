################################################################################
# Multi-Omics Integration for Beta-Cell Workload Analysis
# ========================================================
#
# Integrates scRNA-seq workload scores with proteomics/metabolomics
# using mixOmics (DIABLO, sPLS, MINT)
#
# Author: Beta-Cell Workload Analysis Pipeline
################################################################################

# Load libraries
suppressPackageStartupMessages({
  library(mixOmics)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

# Paths - detect script location robustly
get_script_dir <- function() {
  # Try different methods to get script directory
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  # Fallback to current working directory
  return(getwd())
}

SCRIPT_DIR <- get_script_dir()
# If we're in code/multiomics, go up two levels to workload dir
if (grepl("code[/\\\\]multiomics$", SCRIPT_DIR)) {
  WORKLOAD_DIR <- normalizePath(file.path(SCRIPT_DIR, "../.."))
} else {
  WORKLOAD_DIR <- SCRIPT_DIR
}
RESULTS_DIR <- file.path(WORKLOAD_DIR, "results", "multiomics")
DATA_DIR <- file.path(WORKLOAD_DIR, "data")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat(paste(rep("=", 60), collapse=""), "\n")
cat("MULTI-OMICS INTEGRATION FOR BETA-CELL WORKLOAD\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

################################################################################
# SECTION 1: Load and Prepare Data
################################################################################

cat("SECTION 1: Data Preparation\n")
cat(paste(rep("-", 40), collapse=""), "\n")

# Load workload scores from X-intNMF analysis
load_workload_scores <- function() {
  # Try to load from previous analysis
  cwi_file <- file.path(WORKLOAD_DIR, "results/xintnmf/composite_workload_index.csv")

  if (file.exists(cwi_file)) {
    cwi <- read.csv(cwi_file)
    cat("  Loaded workload scores:", nrow(cwi), "cells\n")
    return(cwi)
  }

  # Generate simulated data if not available
  cat("  Generating simulated workload data for demonstration\n")

  set.seed(42)
  n_cells <- 500

  cwi <- data.frame(
    cell_id = paste0("cell_", 1:n_cells),
    condition = sample(c("ND", "T2D"), n_cells, replace = TRUE, prob = c(0.6, 0.4)),
    CWI = rnorm(n_cells, mean = 0.5, sd = 0.2),
    demand_score = rnorm(n_cells, mean = 0.4, sd = 0.15),
    capacity_score = rnorm(n_cells, mean = 0.6, sd = 0.2),
    stress_score = rnorm(n_cells, mean = 0.3, sd = 0.15),
    dediff_score = rnorm(n_cells, mean = 0.2, sd = 0.1)
  )

  # Make T2D have higher workload
  t2d_idx <- cwi$condition == "T2D"
  cwi$CWI[t2d_idx] <- cwi$CWI[t2d_idx] + 0.3
  cwi$stress_score[t2d_idx] <- cwi$stress_score[t2d_idx] + 0.2

  # Clip to valid range
  cwi$CWI <- pmax(0, pmin(1, cwi$CWI))

  return(cwi)
}

# Simulate multi-omics data (replace with real data when available)
simulate_multiomics_data <- function(cwi) {
  cat("  Simulating multi-omics data for demonstration\n")

  set.seed(42)
  n_cells <- nrow(cwi)

  # Workload-related genes
  workload_genes <- c(
    # Capacity
    "PDX1", "NKX6-1", "MAFA", "GCK", "SLC2A2", "INS", "IAPP",
    "PCSK1", "PCSK2", "CPE", "ABCC8", "KCNJ11",
    # Stress
    "DDIT3", "ATF4", "XBP1", "HSPA5", "ERN1", "EIF2AK3", "ATF6",
    "TRIB3", "CALR", "CANX", "HSP90B1",
    # Dedifferentiation
    "ALDH1A3", "SOX9", "HES1", "NEUROG3"
  )

  # Additional genes
  other_genes <- paste0("Gene_", 1:200)
  all_genes <- c(workload_genes, other_genes)

  # Generate scRNA-seq expression matrix
  X_rna <- matrix(rnorm(n_cells * length(all_genes), mean = 5, sd = 2),
                  nrow = n_cells, ncol = length(all_genes))
  colnames(X_rna) <- all_genes
  rownames(X_rna) <- cwi$cell_id

  # Make expression correlate with workload
  for (gene in workload_genes[1:12]) {  # Capacity genes
    idx <- which(colnames(X_rna) == gene)
    if (length(idx) > 0) {
      X_rna[, idx] <- X_rna[, idx] - cwi$CWI * 2  # Lower in high workload
    }
  }
  for (gene in workload_genes[13:23]) {  # Stress genes
    idx <- which(colnames(X_rna) == gene)
    if (length(idx) > 0) {
      X_rna[, idx] <- X_rna[, idx] + cwi$stress_score * 3  # Higher in stress
    }
  }

  # Generate proteomics data (subset of genes)
  protein_names <- c(workload_genes[1:20], paste0("Protein_", 1:80))
  X_protein <- matrix(rnorm(n_cells * length(protein_names), mean = 10, sd = 3),
                      nrow = n_cells, ncol = length(protein_names))
  colnames(X_protein) <- protein_names
  rownames(X_protein) <- cwi$cell_id

  # Correlate with workload
  for (i in 1:12) {
    X_protein[, i] <- X_protein[, i] - cwi$CWI * 1.5
  }
  for (i in 13:20) {
    X_protein[, i] <- X_protein[, i] + cwi$stress_score * 2
  }

  # Generate metabolomics data
  metabolite_names <- c(
    "Glucose", "Pyruvate", "Lactate", "ATP", "ADP", "NADH", "NAD",
    "Glutamate", "Glutamine", "Proinsulin", "C-peptide",
    paste0("Metabolite_", 1:50)
  )
  X_metab <- matrix(rnorm(n_cells * length(metabolite_names), mean = 8, sd = 2),
                    nrow = n_cells, ncol = length(metabolite_names))
  colnames(X_metab) <- metabolite_names
  rownames(X_metab) <- cwi$cell_id

  # Correlate key metabolites
  X_metab[, "Glucose"] <- X_metab[, "Glucose"] + cwi$CWI * 2
  X_metab[, "ATP"] <- X_metab[, "ATP"] - cwi$stress_score * 2
  X_metab[, "Proinsulin"] <- X_metab[, "Proinsulin"] + cwi$demand_score * 2

  return(list(
    RNA = X_rna,
    Protein = X_protein,
    Metabolomics = X_metab
  ))
}

# Load data
cwi_data <- load_workload_scores()
multiomics <- simulate_multiomics_data(cwi_data)

cat("  RNA-seq features:", ncol(multiomics$RNA), "\n")
cat("  Proteomics features:", ncol(multiomics$Protein), "\n")
cat("  Metabolomics features:", ncol(multiomics$Metabolomics), "\n")

################################################################################
# SECTION 2: Define Workload Groups
################################################################################

cat("\nSECTION 2: Define Workload Groups\n")
cat(paste(rep("-", 40), collapse=""), "\n")

# Create workload group labels (High vs Low)
cwi_data$workload_group <- ifelse(cwi_data$CWI > median(cwi_data$CWI), "High", "Low")
Y <- factor(cwi_data$workload_group)

cat("  High workload cells:", sum(Y == "High"), "\n")
cat("  Low workload cells:", sum(Y == "Low"), "\n")

# Also create condition labels
condition <- factor(cwi_data$condition)

################################################################################
# SECTION 3: Sparse PLS for Pairwise Integration
################################################################################

cat("\nSECTION 3: Pairwise sPLS Analysis\n")
cat(paste(rep("-", 40), collapse=""), "\n")

# sPLS: RNA vs Protein
cat("  Running sPLS: RNA vs Protein\n")

spls_rna_protein <- spls(
  X = multiomics$RNA,
  Y = multiomics$Protein,
  ncomp = 3,
  keepX = c(50, 30, 20),
  keepY = c(30, 20, 15),
  mode = "canonical"
)

# Save correlation plot
pdf(file.path(RESULTS_DIR, "spls_rna_protein_samples.pdf"), width = 10, height = 8)
plotIndiv(spls_rna_protein, comp = c(1, 2),
          group = Y,
          legend = TRUE,
          title = "sPLS: RNA vs Protein - Sample Plot")
dev.off()

pdf(file.path(RESULTS_DIR, "spls_rna_protein_variables.pdf"), width = 10, height = 10)
plotVar(spls_rna_protein, comp = c(1, 2),
        var.names = list(X = TRUE, Y = TRUE),
        cex = c(3, 3),
        title = "sPLS: RNA vs Protein - Variable Plot")
dev.off()

# sPLS: RNA vs Metabolomics
cat("  Running sPLS: RNA vs Metabolomics\n")

spls_rna_metab <- spls(
  X = multiomics$RNA,
  Y = multiomics$Metabolomics,
  ncomp = 3,
  keepX = c(50, 30, 20),
  keepY = c(20, 15, 10),
  mode = "canonical"
)

# Extract correlations
cor_rna_protein <- cor(spls_rna_protein$variates$X[, 1],
                       spls_rna_protein$variates$Y[, 1])
cor_rna_metab <- cor(spls_rna_metab$variates$X[, 1],
                     spls_rna_metab$variates$Y[, 1])

cat("  RNA-Protein correlation (comp 1):", round(cor_rna_protein, 3), "\n")
cat("  RNA-Metabolomics correlation (comp 1):", round(cor_rna_metab, 3), "\n")

################################################################################
# SECTION 4: DIABLO - Multi-Block Discriminant Analysis
################################################################################

cat("\nSECTION 4: DIABLO Multi-Block Analysis\n")
cat(paste(rep("-", 40), collapse=""), "\n")

# Prepare blocks
X_blocks <- list(
  RNA = multiomics$RNA,
  Protein = multiomics$Protein,
  Metabolomics = multiomics$Metabolomics
)

# Design matrix - specify correlations between blocks
# Higher values = stronger expected correlation
design <- matrix(0.1, ncol = 3, nrow = 3,
                 dimnames = list(names(X_blocks), names(X_blocks)))
diag(design) <- 0

# Increase RNA-Protein correlation expectation
design["RNA", "Protein"] <- 0.3
design["Protein", "RNA"] <- 0.3

cat("  Design matrix:\n")
print(design)

# Tune DIABLO parameters
cat("\n  Tuning DIABLO parameters...\n")

# Define keepX values to test
test_keepX <- list(
  RNA = c(20, 40, 60),
  Protein = c(15, 30, 45),
  Metabolomics = c(10, 20, 30)
)

# Run tuning (simplified for speed)
tune_diablo <- tune.block.splsda(
  X = X_blocks,
  Y = Y,
  ncomp = 2,
  test.keepX = test_keepX,
  design = design,
  validation = "Mfold",
  folds = 5,
  nrepeat = 3,
  BPPARAM = BiocParallel::SerialParam()
)

# Get optimal keepX
optimal_keepX <- tune_diablo$choice.keepX
cat("  Optimal keepX:\n")
print(optimal_keepX)

# Run final DIABLO model
cat("\n  Running final DIABLO model...\n")

diablo_result <- block.splsda(
  X = X_blocks,
  Y = Y,
  ncomp = 3,
  keepX = optimal_keepX,
  design = design
)

# Performance evaluation
cat("  Evaluating model performance...\n")
perf_diablo <- tryCatch({
  perf(diablo_result,
       validation = "Mfold",
       folds = 5,
       nrepeat = 10)
}, error = function(e) {
  cat("  Warning: Performance evaluation failed, using default assessment\n")
  NULL
})

if (!is.null(perf_diablo)) {
  # Try to extract error rate from various possible structures
  err_rate <- tryCatch({
    if (!is.null(perf_diablo$WeightedVote.error.rate)) {
      mean(perf_diablo$WeightedVote.error.rate$max.dist[, 2])
    } else if (!is.null(perf_diablo$error.rate.all)) {
      perf_diablo$error.rate.all$overall$max.dist[2]
    } else {
      NA
    }
  }, error = function(e) NA)

  if (!is.na(err_rate)) {
    cat("  Overall error rate (comp 2):", round(err_rate, 3), "\n")
  } else {
    cat("  Model fitted successfully (error rate unavailable)\n")
  }
} else {
  cat("  Model fitted successfully (performance eval skipped)\n")
}

################################################################################
# SECTION 5: DIABLO Visualization
################################################################################

cat("\nSECTION 5: DIABLO Visualization\n")
cat(paste(rep("-", 40), collapse=""), "\n")

# Consensus sample plot (use weighted.average instead of consensus)
pdf(file.path(RESULTS_DIR, "diablo_consensus_samples.pdf"), width = 10, height = 8)
plotIndiv(diablo_result,
          comp = c(1, 2),
          blocks = "weighted.average",
          group = Y,
          legend = TRUE,
          title = "DIABLO Weighted Average - Workload Groups")
dev.off()

# Per-block sample plots
pdf(file.path(RESULTS_DIR, "diablo_block_samples.pdf"), width = 15, height = 5)
par(mfrow = c(1, 3))
for (block in names(X_blocks)) {
  plotIndiv(diablo_result,
            comp = c(1, 2),
            blocks = block,
            group = Y,
            legend = TRUE,
            title = paste("DIABLO -", block))
}
dev.off()

# Variable plot
pdf(file.path(RESULTS_DIR, "diablo_variables.pdf"), width = 12, height = 10)
plotVar(diablo_result,
        comp = c(1, 2),
        var.names = list(RNA = FALSE, Protein = FALSE, Metabolomics = TRUE),
        legend = TRUE,
        title = "DIABLO - Variable Correlation")
dev.off()

# Circos plot showing inter-block correlations
pdf(file.path(RESULTS_DIR, "diablo_circos.pdf"), width = 12, height = 12)
circosPlot(diablo_result, cutoff = 0.5, line = TRUE,
           showIntraLinks = FALSE)
dev.off()
cat("  Saved circos plot\n")

# Heatmap
pdf(file.path(RESULTS_DIR, "diablo_heatmap.pdf"), width = 14, height = 10)
cimDiablo(diablo_result,
          margins = c(10, 15),
          color.blocks = c("#1f77b4", "#ff7f0e", "#2ca02c"))
dev.off()
cat("  Saved heatmap\n")

################################################################################
# SECTION 6: Extract Selected Features
################################################################################

cat("\nSECTION 6: Extract Selected Features\n")
cat(paste(rep("-", 40), collapse=""), "\n")

# Extract selected variables from each block
selected_features <- list()

for (block in names(X_blocks)) {
  # Get variables selected in component 1 and 2
  sel_comp1 <- selectVar(diablo_result, block = block, comp = 1)
  sel_comp2 <- selectVar(diablo_result, block = block, comp = 2)

  # Combine
  features_df <- data.frame(
    feature = c(sel_comp1[[block]]$name, sel_comp2[[block]]$name),
    component = c(rep(1, length(sel_comp1[[block]]$name)),
                  rep(2, length(sel_comp2[[block]]$name))),
    loading = c(sel_comp1[[block]]$value$value.var,
                sel_comp2[[block]]$value$value.var)
  )
  features_df <- features_df[!duplicated(features_df$feature), ]

  selected_features[[block]] <- features_df

  cat("  ", block, ": ", nrow(features_df), " features selected\n", sep = "")
}

# Save selected features
for (block in names(selected_features)) {
  write.csv(selected_features[[block]],
            file.path(RESULTS_DIR, paste0("diablo_selected_", tolower(block), ".csv")),
            row.names = FALSE)
}

# Loadings plots
pdf(file.path(RESULTS_DIR, "diablo_loadings.pdf"), width = 12, height = 15)
par(mfrow = c(3, 1))
for (block in names(X_blocks)) {
  plotLoadings(diablo_result,
               block = block,
               comp = 1,
               contrib = "max",
               method = "median",
               title = paste(block, "- Component 1 Loadings"))
}
dev.off()

################################################################################
# SECTION 7: Workload-Specific Feature Analysis
################################################################################

cat("\nSECTION 7: Workload-Specific Analysis\n")
cat(paste(rep("-", 40), collapse=""), "\n")

# Identify workload signature features
# Define workload-related gene/protein/metabolite sets
capacity_genes <- c("PDX1", "NKX6-1", "MAFA", "GCK", "SLC2A2", "INS", "IAPP")
stress_genes <- c("DDIT3", "ATF4", "XBP1", "HSPA5", "ERN1")
dediff_genes <- c("ALDH1A3", "SOX9")
workload_metabolites <- c("Glucose", "ATP", "Proinsulin", "C-peptide", "Pyruvate")

all_workload_features <- c(capacity_genes, stress_genes, dediff_genes, workload_metabolites)

# Assign category to features
assign_category <- function(feature) {
  if (feature %in% capacity_genes) return("capacity")
  if (feature %in% stress_genes) return("stress")
  if (feature %in% dediff_genes) return("dediff")
  if (feature %in% workload_metabolites) return("metabolite")
  return("other")
}

# Collect workload signature from all blocks
workload_rows <- list()

for (block in names(selected_features)) {
  features_df <- selected_features[[block]]
  workload_subset <- features_df[features_df$feature %in% all_workload_features, ]

  if (nrow(workload_subset) > 0) {
    workload_subset$omics_type <- block
    workload_subset$category <- sapply(workload_subset$feature, assign_category)
    workload_rows[[block]] <- workload_subset
    cat("  Workload", block, "features found:", nrow(workload_subset), "\n")
  }
}

# Save workload signature
if (length(workload_rows) > 0) {
  # Standardize columns: keep feature, component, loading, omics_type, category
  standard_cols <- c("feature", "component", "loading", "omics_type", "category")

  workload_combined <- do.call(rbind, lapply(workload_rows, function(df) {
    # Ensure all standard columns exist
    for (col in standard_cols) {
      if (!(col %in% names(df))) df[[col]] <- NA
    }
    df[, standard_cols]
  }))

  write.csv(workload_combined,
            file.path(RESULTS_DIR, "diablo_workload_signature.csv"),
            row.names = FALSE)
  cat("  Saved workload signature:", nrow(workload_combined), "features\n")
}

################################################################################
# SECTION 8: Create Summary Report
################################################################################

cat("\nSECTION 8: Summary Report\n")
cat(paste(rep("-", 40), collapse=""), "\n")

report <- c(
  "="," ======================================================================",
  "MULTI-OMICS INTEGRATION REPORT - BETA-CELL WORKLOAD",
  "======================================================================",
  "",
  "ANALYSIS SUMMARY",
  "----------------",
  paste("Total cells analyzed:", nrow(cwi_data)),
  paste("High workload cells:", sum(Y == "High")),
  paste("Low workload cells:", sum(Y == "Low")),
  "",
  "OMICS DATA",
  "-----------",
  paste("RNA-seq features:", ncol(multiomics$RNA)),
  paste("Proteomics features:", ncol(multiomics$Protein)),
  paste("Metabolomics features:", ncol(multiomics$Metabolomics)),
  "",
  "sPLS PAIRWISE CORRELATIONS",
  "--------------------------",
  paste("RNA-Protein correlation (comp 1):", round(cor_rna_protein, 3)),
  paste("RNA-Metabolomics correlation (comp 1):", round(cor_rna_metab, 3)),
  "",
  "DIABLO RESULTS",
  "--------------",
  paste("Model error rate:", if (!is.null(perf_diablo)) {
    err <- tryCatch({
      if (!is.null(perf_diablo$WeightedVote.error.rate)) {
        round(mean(perf_diablo$WeightedVote.error.rate$max.dist[, 2]), 3)
      } else if (!is.null(perf_diablo$error.rate.all)) {
        round(perf_diablo$error.rate.all$overall$max.dist[2], 3)
      } else {
        "N/A"
      }
    }, error = function(e) "N/A")
    err
  } else "N/A"),
  "",
  "Selected features per block:",
  paste("  RNA:", nrow(selected_features$RNA)),
  paste("  Protein:", nrow(selected_features$Protein)),
  paste("  Metabolomics:", nrow(selected_features$Metabolomics)),
  "",
  "WORKLOAD SIGNATURE FEATURES",
  "---------------------------"
)

if (exists("workload_combined") && nrow(workload_combined) > 0) {
  for (omics in unique(workload_combined$omics_type)) {
    features <- workload_combined$feature[workload_combined$omics_type == omics]
    report <- c(report, paste(" ", omics, ":", paste(features, collapse = ", ")))
  }
}

report <- c(report, "",
            "======================================================================",
            "FILES GENERATED",
            "======================================================================",
            "- spls_rna_protein_samples.pdf",
            "- spls_rna_protein_variables.pdf",
            "- diablo_consensus_samples.pdf",
            "- diablo_block_samples.pdf",
            "- diablo_variables.pdf",
            "- diablo_circos.pdf",
            "- diablo_heatmap.pdf",
            "- diablo_loadings.pdf",
            "- diablo_selected_rna.csv",
            "- diablo_selected_protein.csv",
            "- diablo_selected_metabolomics.csv",
            "- diablo_workload_signature.csv",
            "")

writeLines(report, file.path(RESULTS_DIR, "mixomics_integration_report.txt"))

cat("\n")
cat(paste(report, collapse = "\n"))
cat("\n\nResults saved to:", RESULTS_DIR, "\n")

################################################################################
# SECTION 9: Export for Downstream Analysis
################################################################################

cat("\nSECTION 9: Export for Downstream Analysis\n")
cat(paste(rep("-", 40), collapse=""), "\n")

# Export DIABLO scores for pathway analysis
diablo_scores <- data.frame(
  cell_id = cwi_data$cell_id,
  condition = cwi_data$condition,
  workload_group = cwi_data$workload_group,
  CWI = cwi_data$CWI,
  DIABLO_comp1 = diablo_result$variates$RNA[, 1],
  DIABLO_comp2 = diablo_result$variates$RNA[, 2]
)

write.csv(diablo_scores,
          file.path(RESULTS_DIR, "diablo_cell_scores.csv"),
          row.names = FALSE)

cat("  Exported DIABLO scores for", nrow(diablo_scores), "cells\n")

# Export for pathway enrichment
rna_for_enrichment <- selected_features$RNA[order(abs(selected_features$RNA$loading),
                                                   decreasing = TRUE), ]
write.csv(rna_for_enrichment,
          file.path(RESULTS_DIR, "diablo_rna_for_pathway_analysis.csv"),
          row.names = FALSE)

cat("  Exported RNA features for pathway analysis\n")

cat("\n", "=", rep(60), "\n")
cat("MULTI-OMICS INTEGRATION COMPLETE\n")
cat("=", rep(60), "\n")
