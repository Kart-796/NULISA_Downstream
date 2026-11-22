# =============================================================================
# Shared configuration loader for all R scripts.
#
# Usage (at the top of every script, run from the repository root):
#   source("utils/load_config.R")
#   setwd(cfg_path("multiome_workdir"))
#   sample_info <- read_samples("multiome")
#
# Settings live in config/config.yaml (git-ignored; copy config/config.example.yaml).
# Set the environment variable PROJECT_CONFIG to use a different config file.
# Relative paths in the config are resolved against the repository root.
# =============================================================================
if (!requireNamespace("yaml", quietly = TRUE)) stop("Please install the 'yaml' R package")

.repo_root <- normalizePath(Sys.getenv("PROJECT_ROOT", getwd()), mustWork = TRUE)
.cfg_file  <- Sys.getenv("PROJECT_CONFIG", file.path(.repo_root, "config", "config.yaml"))
if (!file.exists(.cfg_file)) {
  stop("Config file not found: ", .cfg_file,
       "\nCopy config/config.example.yaml to config/config.yaml and edit it,",
       "\nand run scripts from the repository root (or set PROJECT_ROOT).")
}
cfg <- yaml::read_yaml(.cfg_file)

.resolve <- function(p) {
  p <- path.expand(p)
  if (!grepl("^(/|[A-Za-z]:)", p)) p <- file.path(.repo_root, p)
  normalizePath(p, mustWork = FALSE)
}

# Absolute path for a key under `paths:`; creates the folder if create = TRUE.
cfg_path <- function(key, create = FALSE) {
  p <- cfg$paths[[key]]
  if (is.null(p) || !nzchar(p)) stop("Missing 'paths: ", key, "' in ", .cfg_file)
  p <- .resolve(p)
  if (create) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

# Sample sheet for a key under `samples:` (see config/samples_*.csv).
read_samples <- function(key) {
  f <- cfg$samples[[key]]
  if (is.null(f)) stop("Missing 'samples: ", key, "' in ", .cfg_file)
  read.csv(.resolve(f), stringsAsFactors = FALSE, check.names = FALSE)
}
