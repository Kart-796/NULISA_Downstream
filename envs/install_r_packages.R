# =============================================================================
# Install the R packages used by analysis/01_differential_abundance_nulisaseqr.R.
#   Rscript envs/install_r_packages.R
# NULISAseqR is installed from GitHub (Alamar Biosciences); see
# https://nulisaseqr.alamarbio.com for platform-specific prerequisites.
# Tip: afterwards run renv::init() / renv::snapshot() and commit renv.lock.
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))    install.packages("devtools")

install.packages(setdiff(c("yaml", "tidyverse", "readxl", "ggplot2", "ggrepel"),
                         rownames(installed.packages())))
BiocManager::install(setdiff("ComplexHeatmap", rownames(installed.packages())), update = FALSE, ask = FALSE)
if (!requireNamespace("NULISAseqR", quietly = TRUE)) devtools::install_github("Alamar-Biosciences/NULISAseqR")
