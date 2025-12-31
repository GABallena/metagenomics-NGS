library(ggplot2)
library(dplyr)

# Parse CLI arguments or env vars for I/O
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else Sys.getenv("CONTIGS_DIR", "spades_output")
output_dir <- if (length(args) >= 2) args[2] else Sys.getenv("CONTIGS_PLOT_DIR", ".")

data_dir <- normalizePath(data_dir, mustWork = FALSE)
output_dir <- normalizePath(output_dir, mustWork = FALSE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to extract contig data from a FASTA file
parse_contigs <- function(fasta_file) {
  # Print debug info
  cat("Processing file:", fasta_file, "\n")
  
  if (!file.exists(fasta_file)) {
    warning("File does not exist:", fasta_file)
    return(NULL)
  }
  
  contigs <- data.frame(Length = numeric(), Coverage = numeric(), Sample = character(), 
                       stringsAsFactors = FALSE)
  
  # Extract sample name from the path
  sample_name <- basename(dirname(fasta_file))
  cat("Sample name:", sample_name, "\n")
  
  # Read the entire file
  lines <- readLines(fasta_file, warn = FALSE)
  
  # Process headers
  headers <- lines[grep("^>", lines)]
  if (length(headers) == 0) {
    warning("No contig headers found in:", fasta_file)
    return(NULL)
  }
  
  for (header in headers) {
    # Extract length and coverage using regex
    length_match <- regexpr("length_\\d+", header)
    cov_match <- regexpr("cov_[0-9.]+", header)
    
    if (length_match > 0 && cov_match > 0) {
      length_str <- substr(header, length_match, 
                         length_match + attr(length_match, "match.length") - 1)
      cov_str <- substr(header, cov_match, 
                       cov_match + attr(cov_match, "match.length") - 1)
      
      length <- as.numeric(sub("length_", "", length_str))
      coverage <- as.numeric(sub("cov_", "", cov_str))
      
      contigs <- rbind(contigs, 
                      data.frame(Length = length, 
                                Coverage = coverage, 
                                Sample = sample_name,
                                stringsAsFactors = FALSE))
    }
  }
  
  if (nrow(contigs) == 0) {
    warning("No valid contigs found in:", fasta_file)
    return(NULL)
  }
  
  cat("Found", nrow(contigs), "contigs in", sample_name, "\n")
  return(contigs)
}

# Find all relevant FASTA files
fasta_files <- list.files(data_dir, 
                         pattern = "*_final_filtered_contigs.fa$", 
                         full.names = TRUE, 
                         recursive = TRUE)

cat("Found", length(fasta_files), "FASTA files\n")
if (length(fasta_files) == 0) {
  stop("No FASTA files found in: ", data_dir, "\n",
       "Pass the contig directory as the first CLI argument or set CONTIGS_DIR.")
}

# Parse all contigs from the files
all_contigs_list <- lapply(fasta_files, parse_contigs)
all_contigs_list <- Filter(Negate(is.null), all_contigs_list)  # Remove NULL results

if (length(all_contigs_list) == 0) {
  stop("No valid contig data found in any files")
}

# Combine all data
all_contigs <- bind_rows(all_contigs_list)

# Verify we have data
cat("Total contigs:", nrow(all_contigs), "\n")
if (nrow(all_contigs) == 0) {
  stop("No contigs data after combining")
}

# Filter data
all_contigs <- all_contigs %>%
  filter(Length >= 200)

# Calculate summary statistics for each sample
summary_stats <- all_contigs %>%
  group_by(Sample) %>%
  summarise(
    mean_coverage = mean(Coverage),
    sd_coverage = sd(Coverage),
    median_coverage = median(Coverage)
  ) %>%
  arrange(desc(median_coverage))  # Sort by median coverage

# Create Manhattan-style plot
p <- ggplot(all_contigs, aes(x = Sample, y = Coverage)) +
  geom_jitter(aes(color = Length), 
              alpha = 0.6, 
              size = 1.5,
              width = 0.3,  # Control horizontal jitter
              height = 0) + # No vertical jitter to preserve true coverage values
  scale_color_viridis_c(option = "magma", 
                        trans = "log10",
                        name = "Contig Length\n(log10 scale)") +
  scale_y_log10(labels = scales::comma) +
  labs(
    title = "Contig Coverage Distribution by Sample",
    subtitle = "Points colored by contig length",
    x = "Sample",
    y = "Coverage Depth (log10 scale)"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "ghostwhite", color = NA),
    panel.grid.major.x = element_blank(),  # Remove vertical gridlines
    panel.grid.minor.x = element_blank(),  # Remove minor vertical gridlines
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_line(color = "gray95"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10)
  )

# Save the plot
ggsave(file.path(output_dir, "contig_coverage_distribution.png"), 
       plot = p,
       width = 12,  # Made wider to accommodate all samples
       height = 7, 
       dpi = 300, 
       bg = "white")

# Function to create a violin plot with points
create_violin_plot <- function(data) {
  ggplot(data, aes(x = reorder(Sample, Coverage, FUN = median), y = Coverage)) +
    geom_violin(aes(fill = Sample), alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA) +
    geom_jitter(size = 0.5, alpha = 0.3, width = 0.2) +
    scale_y_log10(labels = scales::comma) +
    scale_fill_viridis_d(option = "mako") +
    labs(
      title = "Contig Coverage Distribution",
      subtitle = "Violin plot with embedded boxplot and individual points",
      x = "Sample",
      y = "Coverage Depth (log10 scale)"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "ghostwhite", color = NA),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

# Function to create a length vs coverage scatter plot
create_scatter_plot <- function(data) {
  ggplot(data, aes(x = Length, y = Coverage)) +
    geom_point(aes(color = Sample), alpha = 0.7, size = 1) +
    geom_smooth(method = "loess", color = "black", size = 0.5, alpha = 0.2) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    scale_color_viridis_d(option = "turbo") +
    facet_wrap(~Sample, scales = "free") +
    labs(
      title = "Contig Length vs Coverage",
      subtitle = "Relationship between contig length and coverage depth",
      x = "Contig Length (log10 scale)",
      y = "Coverage Depth (log10 scale)"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "ghostwhite", color = NA),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "none"
    )
}

# Function to create a density ridgeline plot
create_ridgeline_plot <- function(data) {
  library(ggridges)
  ggplot(data, aes(x = Length, y = Sample, fill = Sample)) +
    geom_density_ridges2(alpha = 0.7, scale = 2) +
    scale_x_log10(labels = scales::comma) +
    scale_fill_viridis_d(option = "cividis") +
    labs(
      title = "Contig Length Distribution",
      subtitle = "Density ridgeline plot showing length distribution by sample",
      x = "Contig Length (log10 scale)",
      y = "Sample"
    ) +
    theme_ridges() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
}

# Create and save all plots
violin_plot <- create_violin_plot(all_contigs)
scatter_plot <- create_scatter_plot(all_contigs)
ridgeline_plot <- create_ridgeline_plot(all_contigs)

# Save the new plots
ggsave(file.path(output_dir, "contig_violin_plot.png"), 
       plot = violin_plot,
       width = 12, 
       height = 7, 
       dpi = 300, 
       bg = "white")

ggsave(file.path(output_dir, "contig_scatter_plot.png"), 
       plot = scatter_plot,
       width = 15, 
       height = 10, 
       dpi = 300, 
       bg = "white")

ggsave(file.path(output_dir, "contig_ridgeline_plot.png"), 
       plot = ridgeline_plot,
       width = 12, 
       height = 8, 
       dpi = 300, 
       bg = "white")
