# ONT library planning helper (portfolio-ready)
#
# Goal:
#   Given per-sample Qubit concentrations (ng/µL), find candidate 8-sample pools ("runs")
#   that can meet a target total DNA mass (default: 3 µg) while respecting a max pipette
#   volume per sample (default: 50 µL). Then (optionally) search for 4 runs that together
#   cover all samples at least once.
#
# Caveat:
#   The exact pooling rules vary between labs. This script implements a reasonable,
#   parameterized default. Adjust target_mass_ug, per_sample_mass_mode, and max_vol_uL
#   to match your SOP.

rm(list = ls())
set.seed(2025)

suppressPackageStartupMessages({
  library(gtools)
})

# --- Inputs ---
qubit <- c(
  4.49, 59.70, 89.60, 16.80, 85.70, 40.80, 74.40, 34.50,
  46.80, 1.65, 60.00, 56.00, 17.50, 92.60, 74.00, 15.70,
  60.00, 55.00, 58.00, 47.10, 56.00, 53.00, 49.20
)
sample_ids <- paste0("S", seq_along(qubit))

# --- Parameters ---
target_mass_ug <- 3.0     # total mass per 8-sample run
max_vol_uL     <- 50.0    # max volume per sample
pool_size      <- 8
runs_to_cover  <- 4

target_mass_ng <- target_mass_ug * 1000
per_sample_mass_ng <- target_mass_ng / pool_size

# Required volume per sample to contribute equal mass:
required_vol_uL <- function(conc_ng_uL, mass_ng) mass_ng / conc_ng_uL

is_valid_run <- function(indices) {
  conc <- qubit[indices]
  # guard zero/NA
  if (any(is.na(conc)) || any(conc <= 0)) return(FALSE)
  vols <- required_vol_uL(conc, per_sample_mass_ng)
  all(vols <= max_vol_uL)
}

# --- Step 1: Generate all 8-sample combos and keep those that are valid ---
combo_indices <- combinations(n = length(qubit), r = pool_size, repeats.allowed = FALSE)

valid_mask <- apply(combo_indices, 1, is_valid_run)
valid_runs <- combo_indices[valid_mask, , drop = FALSE]
cat(sprintf("Valid runs: %d / %d\n", nrow(valid_runs), nrow(combo_indices)))

# --- Step 2: Search for sets of 4 runs that cover all samples at least once ---
# Brute force is combinatorially huge. We use a randomized greedy search.
all_samples <- seq_along(qubit)

score_cover <- function(run_set) {
  cov <- sort(unique(as.vector(t(run_set))))
  length(intersect(cov, all_samples))
}

find_cover_sets <- function(n_iter = 2000) {
  best <- NULL
  best_cov <- -1
  for (i in seq_len(n_iter)) {
    # random start
    idx <- sample(seq_len(nrow(valid_runs)), runs_to_cover, replace = FALSE)
    rs <- valid_runs[idx, , drop = FALSE]
    cov <- score_cover(rs)
    if (cov > best_cov) {
      best <- rs
      best_cov <- cov
      if (best_cov == length(all_samples)) break
    }
  }
  list(best = best, covered = best_cov)
}

res <- find_cover_sets(n_iter = 5000)
cat(sprintf("Best coverage with %d runs: %d / %d samples\n",
            runs_to_cover, res$covered, length(all_samples)))

if (!is.null(res$best)) {
  for (j in seq_len(nrow(res$best))) {
    ids <- sample_ids[res$best[j, ]]
    vols <- required_vol_uL(qubit[res$best[j, ]], per_sample_mass_ng)
    cat(sprintf("\nRun %d: %s\n", j, paste(ids, collapse = ", ")))
    cat(sprintf("Volumes (uL): %s\n", paste(sprintf("%.2f", vols), collapse = ", ")))
  }
}

# Export valid runs (optional)
# write.table(valid_runs, file="valid_runs.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
