# Bracken genus-level figures + metadata-aware summaries
# - Scans bracken_output/*.genus.bracken
# - Reads metadata from "metadata/masterlist.tsv"
# - Produces: stacked barplots (per-sample, per-group), alpha diversity (Shannon, Simpson, richness),
#             PCoA (Bray-Curtis), heatmap of top genera
# - Writes a missing-data report; does not simulate or impute data

options(stringsAsFactors = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

# Paths
root_dir <- getwd()
bracken_dir <- file.path(root_dir, "bracken_output")
metadata_file <- file.path(root_dir, "metadata/masterlist.tsv")
out_dir <- file.path(root_dir, "diversity_results", "bracken_genus_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

report_path <- file.path(out_dir, "missing_data_report.md")
summary_path <- file.path(out_dir, "figures_generated.txt")

# Helper: safe package check
has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

use_gg <- has_pkg("ggplot2")
if (use_gg) {
  library(ggplot2)
}

# Helper: read a Bracken genus file as a named vector of relative abundance
read_bracken_genus <- function(path) {
  df <- tryCatch(read.delim(path, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  req_cols <- c("name", "fraction_total_reads")
  if (!all(req_cols %in% names(df))) return(NULL)
  x <- df$fraction_total_reads
  names(x) <- df$name
  x[is.na(x)] <- 0
  x
}

# Helper: Shannon, Simpson, Richness
alpha_metrics <- function(p) {
  # Ensure numeric vector p >= 0, sum may be <= 1 if truncation; renormalize for Shannon/Simpson
  p[p < 0] <- 0
  s <- sum(p)
  if (s > 0) p <- p / s
  p_pos <- p[p > 0]
  shannon <- if (length(p_pos) > 0) -sum(p_pos * log(p_pos)) else NA_real_
  simpson <- if (length(p) > 0) 1 - sum(p^2) else NA_real_
  richness <- sum(p > 0)
  c(Shannon = shannon, Simpson = simpson, Richness = richness)
}

# Helper: Bray-Curtis distance between two numeric vectors
bray_curtis <- function(x, y) {
  # x and y non-negative
  x[x < 0] <- 0; y[y < 0] <- 0
  sx <- sum(x); sy <- sum(y)
  if (sx == 0 && sy == 0) return(0)
  if (sx == 0 || sy == 0) return(1)
  # Work on relative abundances
  x <- x / sx; y <- y / sy
  1 - (2 * sum(pmin(x, y))) / (sum(x) + sum(y))
}

# Load Bracken genus files
log_msg("Scanning Bracken genus files in ", bracken_dir)
if (!dir.exists(bracken_dir)) {
  stop(sprintf("Bracken directory not found: %s", bracken_dir))
}
bracken_files <- list.files(bracken_dir, pattern = "\\.genus\\.bracken$", full.names = TRUE)

if (length(bracken_files) == 0) {
  msg <- sprintf("No genus-level Bracken files found in %s. Exiting.", bracken_dir)
  writeLines(c("# Missing Data Report", "", msg), con = report_path)
  stop(msg)
}

# Build abundance list per sample
ab_list <- list()
sample_codes <- character()
for (f in bracken_files) {
  bn <- basename(f)
  sc <- sub("_.*$", "", bn) # e.g., SAMPLE-01_R1 -> SAMPLE-01
  x <- read_bracken_genus(f)
  if (is.null(x)) {
    log_msg("Skipping (unreadable or missing cols): ", bn)
    next
  }
  ab_list[[sc]] <- x
  sample_codes <- c(sample_codes, sc)
}

sample_codes <- unique(sample_codes)
if (length(ab_list) == 0) {
  msg <- "No readable Bracken genus files with required columns (name, fraction_total_reads)."
  writeLines(c("# Missing Data Report", "", msg), con = report_path)
  stop(msg)
}

# Make union of genera and build matrix
all_genera <- unique(unlist(lapply(ab_list, names)))
mat <- matrix(0, nrow = length(all_genera), ncol = length(ab_list),
              dimnames = list(all_genera, names(ab_list)))
for (sc in names(ab_list)) {
  v <- ab_list[[sc]]
  mat[names(v), sc] <- as.numeric(v)
}

# Save abundance matrix
abundance_csv <- file.path(out_dir, "bracken_genus_relative_abundance.csv")
write.csv(round(mat, 8), abundance_csv, row.names = TRUE)

# Alpha diversity per sample
alpha_df <- do.call(rbind, lapply(seq_len(ncol(mat)), function(j) {
  m <- mat[, j]
  a <- alpha_metrics(m)
  data.frame(SampleCode = colnames(mat)[j], t(a), check.names = FALSE)
}))
alpha_csv <- file.path(out_dir, "alpha_diversity.csv")
write.csv(alpha_df, alpha_csv, row.names = FALSE)

# Read metadata
log_msg("Reading metadata: ", metadata_file)
md <- tryCatch(read.delim(metadata_file, check.names = FALSE, fill = TRUE), error = function(e) NULL)
md_ok <- TRUE
if (is.null(md) || nrow(md) == 0) {
  md_ok <- FALSE
  log_msg("Metadata not found or unreadable.")
}

sample_col <- NA_character_
group_col <- NA_character_
if (md_ok) {
  # Remove completely empty columns
  keep <- vapply(md, function(col) any(!is.na(col) & trimws(as.character(col)) != ""), logical(1))
  md <- md[, keep, drop = FALSE]
  cn <- names(md)
  cn_norm <- toupper(trimws(cn))
  sample_col <- cn[which.max(cn_norm == "SAMPLE CODE")] # returns "" if none
  if (!any(cn_norm == "SAMPLE CODE")) {
    md_ok <- FALSE
    log_msg("Metadata is missing 'SAMPLE CODE' column; will proceed without metadata joins.")
  } else {
    names(md)[which(cn_norm == "SAMPLE CODE")] <- "SAMPLE CODE"
    sample_col <- "SAMPLE CODE"
    # Try to identify a grouping column
    if (any(cn_norm == "SAMPLE TYPE")) {
      names(md)[which(cn_norm == "SAMPLE TYPE")] <- "SAMPLE TYPE"
      group_col <- "SAMPLE TYPE"
    } else {
      # Fallback: none
      group_col <- NA_character_
    }
  }
}

md_use <- NULL
if (md_ok) {
  md_use <- md[, c(sample_col, setdiff(names(md), sample_col)), drop = FALSE]
  names(md_use)[1] <- "SampleCode"
  # Trim whitespace
  md_use$SampleCode <- trimws(as.character(md_use$SampleCode))
  # Keep only non-empty SampleCodes
  md_use <- md_use[!is.na(md_use$SampleCode) & md_use$SampleCode != "", , drop = FALSE]
}

# Missing data report
missing_lines <- c("# Missing Data Report", "")
present_ab <- colnames(mat)
if (md_ok) {
  present_md <- unique(md_use$SampleCode)
  only_in_ab <- setdiff(present_ab, present_md)
  only_in_md <- setdiff(present_md, present_ab)
  missing_lines <- c(missing_lines,
                     sprintf("Samples in Bracken but missing in metadata (%d):", length(only_in_ab)),
                     if (length(only_in_ab)) paste0("- ", sort(only_in_ab)) else "- None",
                     "",
                     sprintf("Samples in metadata but missing Bracken genus outputs (%d):", length(only_in_md)),
                     if (length(only_in_md)) paste0("- ", sort(only_in_md)) else "- None",
                     "")
  if (!is.na(group_col)) {
    # Summarize missing group values
    md_subset <- md_use[md_use$SampleCode %in% present_ab, , drop = FALSE]
    if (!(group_col %in% names(md_subset))) {
      missing_lines <- c(missing_lines, "Grouping column 'SAMPLE TYPE' not found in metadata.")
    } else {
      grp_na <- md_subset$`SAMPLE TYPE`
      na_samps <- md_subset$SampleCode[is.na(grp_na) | trimws(as.character(grp_na)) == ""]
      missing_lines <- c(missing_lines,
                         sprintf("Samples missing SAMPLE TYPE among those with Bracken data (%d):", length(na_samps)),
                         if (length(na_samps)) paste0("- ", sort(na_samps)) else "- None",
                         "")
    }
  }
} else {
  missing_lines <- c(missing_lines, "Metadata file missing or lacks 'SAMPLE CODE'. Figures will be generated without metadata grouping where possible.")
}
writeLines(missing_lines, con = report_path)

# Convenience: sample order by SampleCode
sample_order <- sort(colnames(mat))
mat <- mat[, sample_order, drop = FALSE]

# Figure creation helpers
save_plot <- function(file, width=9, height=6, expr) {
  ext <- tolower(tools::file_ext(file))
  if (ext == "png") {
    png(file, width = width, height = height, units = "in", res = 150)
  } else if (ext == "pdf") {
    pdf(file, width = width, height = height)
  } else if (ext == "svg" && has_pkg("svglite")) {
    svglite::svglite(file, width = width, height = height)
  } else {
    png(paste0(file, ".png"), width = width, height = height, units = "in", res = 150)
  }
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

generated <- character()

# 1) Stacked barplot of top genera per sample
overall_ab <- rowMeans(mat, na.rm = TRUE)
 topN <- 20
keep_genera <- names(sort(overall_ab, decreasing = TRUE))[seq_len(min(topN, length(overall_ab)))]
mat_top <- rbind(mat[keep_genera, , drop = FALSE], Other = pmax(0, 1 - colSums(mat[keep_genera, , drop = FALSE])))

if (use_gg) {
  library(reshape2)
  df_long <- reshape2::melt(mat_top, varnames = c("Genus", "Sample"), value.name = "Abundance")
  df_long$Sample <- factor(df_long$Sample, levels = sample_order)
  p <- ggplot(df_long, aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", width = 0.9) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = sprintf("Top %d genera per sample (Bracken)", topN), x = "Sample", y = "Relative abundance") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  f1 <- file.path(out_dir, "stacked_bar_top_genera.png")
  ggsave(f1, p, width = 12, height = 6, dpi = 200)
  generated <- c(generated, f1)
} else {
  # Base R stacked barplot
  f1 <- file.path(out_dir, "stacked_bar_top_genera.png")
  save_plot(f1, width = 12, height = 6, expr = {
    par(mar = c(8, 4, 3, 1))
    barplot(mat_top, col = rainbow(nrow(mat_top)), las = 2, cex.names = 0.7,
            ylab = "Relative abundance", xlab = "Sample", legend.text = rownames(mat_top), args.legend = list(x = "topright", cex = 0.6))
    title(main = sprintf("Top %d genera per sample (Bracken)", topN))
  })
  generated <- c(generated, f1)
}

# 2) Alpha diversity boxplots by group (if metadata available)
if (md_ok) {
  alpha_merged <- merge(alpha_df, md_use, by.x = "SampleCode", by.y = "SampleCode", all.x = TRUE)
  # Save merged table
  write.csv(alpha_merged, file.path(out_dir, "alpha_with_metadata.csv"), row.names = FALSE)
  if (!is.na(group_col) && ("SAMPLE TYPE" %in% names(alpha_merged))) {
    alpha_merged$`SAMPLE TYPE` <- as.factor(alpha_merged$`SAMPLE TYPE`)
    # Shannon
    f2 <- file.path(out_dir, "alpha_shannon_by_sample_type.png")
    if (use_gg) {
      p2 <- ggplot(alpha_merged, aes(x = `SAMPLE TYPE`, y = Shannon, fill = `SAMPLE TYPE`)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.2, height = 0, size = 1.8, alpha = 0.9) +
        theme_minimal(base_size = 11) +
        labs(title = "Alpha diversity (Shannon) by sample type", x = "Sample type", y = "Shannon index") +
        theme(legend.position = "none")
      ggsave(f2, p2, width = 9, height = 5, dpi = 200)
    } else {
      save_plot(f2, width = 9, height = 5, expr = {
        par(mar = c(8, 4, 3, 1))
        boxplot(Shannon ~ `SAMPLE TYPE`, data = alpha_merged, las = 2, ylab = "Shannon index", col = "lightblue")
        title(main = "Alpha diversity (Shannon) by sample type")
      })
    }
    generated <- c(generated, f2)
  }
}

# 3) PCoA (Bray-Curtis) colored by group if available
if (ncol(mat) >= 2) {
  # Compute distance matrix
  n <- ncol(mat)
  d <- matrix(0, n, n, dimnames = list(colnames(mat), colnames(mat)))
  for (i in seq_len(n)) for (j in seq_len(n)) if (j > i) {
    d[i, j] <- d[j, i] <- bray_curtis(mat[, i], mat[, j])
  }
  # cmdscale requires 'dist' object upper triangle
  dd <- as.dist(d)
  pc <- tryCatch(cmdscale(dd, k = 2, eig = TRUE), error = function(e) NULL)
  if (!is.null(pc)) {
    pts <- as.data.frame(pc$points)
    names(pts) <- c("PCoA1", "PCoA2")
    pts$SampleCode <- rownames(pc$points)
    if (md_ok) {
      pts <- merge(pts, md_use[, c("SampleCode", if (!is.na(group_col)) "SAMPLE TYPE" else NULL), drop = FALSE], by = "SampleCode", all.x = TRUE)
    }
    f3 <- file.path(out_dir, "pcoa_bray_curtis.png")
    if (use_gg) {
      p3 <- ggplot(pts, aes(x = PCoA1, y = PCoA2, label = SampleCode)) +
        {if (!is.null(pts$`SAMPLE TYPE`)) aes(color = `SAMPLE TYPE`) else aes()} +
        geom_point(size = 2) +
        geom_text(vjust = -1, size = 2.7) +
        theme_minimal(base_size = 11) +
        labs(title = "PCoA (Bray-Curtis) at genus level", x = "PCoA1", y = "PCoA2")
      ggsave(f3, p3, width = 7, height = 6, dpi = 200)
    } else {
      save_plot(f3, width = 7, height = 6, expr = {
        par(mar = c(4, 4, 2, 1))
        plot(pts$PCoA1, pts$PCoA2, pch = 19, xlab = "PCoA1", ylab = "PCoA2", main = "PCoA (Bray-Curtis)")
        text(pts$PCoA1, pts$PCoA2, labels = pts$SampleCode, pos = 3, cex = 0.6)
      })
    }
    generated <- c(generated, f3)
  }
}

# 4) Heatmap of top 50 genera
 topH <- 50
keep_h <- names(sort(overall_ab, decreasing = TRUE))[seq_len(min(topH, length(overall_ab)))]
mat_h <- mat[keep_h, , drop = FALSE]
f4 <- file.path(out_dir, "heatmap_top_genera.png")
save_plot(f4, width = 10, height = 10, expr = {
  # Scale rows to improve contrast
  m <- mat_h
  # Avoid division by zero
  rs <- rowSums(m)
  rs[rs == 0] <- 1
  m_scaled <- m / rs
  par(mar = c(6, 8, 2, 2))
  heatmap(m_scaled, Colv = NA, scale = "none", col = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100),
          xlab = "Samples", ylab = "Genera")
})
generated <- c(generated, f4)

# 5) Group-aggregated composition (if metadata available)
if (md_ok && !is.na(group_col) && ("SAMPLE TYPE" %in% names(md_use))) {
  md_subset <- md_use[, c("SampleCode", "SAMPLE TYPE"), drop = FALSE]
  md_subset$SampleCode <- as.character(md_subset$SampleCode)
  md_subset$`SAMPLE TYPE` <- as.character(md_subset$`SAMPLE TYPE`)
  md_subset <- md_subset[md_subset$SampleCode %in% colnames(mat) & !is.na(md_subset$`SAMPLE TYPE`) & trimws(md_subset$`SAMPLE TYPE`) != "", , drop = FALSE]
  if (nrow(md_subset) > 0) {
    groups <- unique(md_subset$`SAMPLE TYPE`)
    agg <- sapply(groups, function(g) {
      samp <- md_subset$SampleCode[md_subset$`SAMPLE TYPE` == g]
      if (length(samp) == 0) return(rowMeans(mat[, 0, drop = FALSE]))
      rowMeans(mat[, samp, drop = FALSE], na.rm = TRUE)
    })
    # Normalize columns to sum to 1
    cs <- colSums(agg)
    cs[cs == 0] <- 1
    agg <- sweep(agg, 2, cs, "/")
    # Keep top genera again
    agg_top <- rbind(agg[keep_genera, , drop = FALSE], Other = pmax(0, 1 - colSums(agg[keep_genera, , drop = FALSE])))
    if (use_gg) {
      library(reshape2)
      df_g <- reshape2::melt(agg_top, varnames = c("Genus", "Group"), value.name = "Abundance")
      p5 <- ggplot(df_g, aes(x = Group, y = Abundance, fill = Genus)) +
        geom_bar(stat = "identity", width = 0.8) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        theme_minimal(base_size = 11) +
        labs(title = "Group-aggregated composition (genus)", x = "Sample type", y = "Relative abundance")
      f5 <- file.path(out_dir, "stacked_bar_group_aggregated.png")
      ggsave(f5, p5, width = 9, height = 5.5, dpi = 200)
      generated <- c(generated, f5)
    } else {
      f5 <- file.path(out_dir, "stacked_bar_group_aggregated.png")
      save_plot(f5, width = 9, height = 5.5, expr = {
        par(mar = c(8, 4, 3, 1))
        barplot(agg_top, col = rainbow(nrow(agg_top)), las = 2, cex.names = 0.8,
                ylab = "Relative abundance", xlab = "Sample type", legend.text = rownames(agg_top), args.legend = list(x = "topright", cex = 0.6))
        title(main = "Group-aggregated composition (genus)")
      })
      generated <- c(generated, f5)
    }
  }
}

# Write summary of generated figures
if (length(generated) == 0) {
  writeLines("No figures were generated (insufficient data). See missing_data_report.md.", con = summary_path)
} else {
  writeLines(c("Figures generated:", paste0("- ", generated)), con = summary_path)
}

log_msg("Done. Outputs in:", out_dir)
