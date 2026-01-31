# Install required packages for mixOmics analysis
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install mixOmics from Bioconductor
if (!requireNamespace("mixOmics", quietly = TRUE)) {
  BiocManager::install("mixOmics", ask = FALSE, update = FALSE)
}

# Install other dependencies
packages <- c("ggplot2", "reshape2", "pheatmap", "RColorBrewer", "gridExtra")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

cat("All packages installed successfully!\n")
