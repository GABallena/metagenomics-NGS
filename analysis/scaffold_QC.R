# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(reshape2)
library(cowplot)
library(scales)
library(ggExtra)
library(RColorBrewer)

# setwd("/path/to/project")  # user sets working directory externally


# Function to parse FASTA headers and extract metrics
parse_fasta_headers <- function(fasta_file) {
  data <- readLines(fasta_file)
  headers <- data[grepl('^>', data)]
  metrics <- lapply(headers, function(header) {
    match <- regexpr('NODE_\\d+_length_(\\d+)_cov_([0-9.]+)', header)
    if (match > 0) {
      values <- regmatches(header, regexec('NODE_\\d+_length_(\\d+)_cov_([0-9.]+)', header))[[1]]
      return(data.frame(Length = as.numeric(values[2]), Coverage = as.numeric(values[3]), stringsAsFactors = FALSE))
    } else {
      return(NULL)
    }
  })
  metrics <- Filter(Negate(is.null), metrics)
  if (length(metrics) == 0) {
    return(data.frame(Length = numeric(0), Coverage = numeric(0)))
  }
  metrics <- do.call(rbind, metrics)
  return(metrics)
}

# Publication theme settings
publication_theme <- theme_bw() +  # Changed from theme_minimal() to theme_bw()
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90"),
    panel.background = element_rect(fill = "white", color = "black"),  # Added explicit white background
    plot.background = element_rect(fill = "white", color = NA)        # Added white plot background
  )

# Single function to create and save combined plot
create_qc_plot <- function(scaffold_df, contig_df, output_dir, sample_name) {
  # Prepare combined data
  combined_df <- rbind(
    scaffold_df %>% mutate(Source = "Scaffolds"),
    contig_df %>% mutate(Source = "Contigs")
  )
  
  # Define consistent colors
  plot_colors <- brewer.pal(3, "Set1")[1:2]
  
  # Length Violin Plot
  length_plot <- ggplot(combined_df, aes(x = Source, y = Length, fill = Source)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_jitter(color = "grey30", size = 0.3, alpha = 0.4, width = 0.15) +
    scale_fill_manual(values = plot_colors) +
    scale_y_log10(labels = comma_format()) +
    labs(x = "Assembly Type", y = "Contig Length (bp)") +
    publication_theme +
    theme(legend.position = "none")

  # Coverage Violin Plot
  coverage_plot <- ggplot(combined_df, aes(x = Source, y = Coverage, fill = Source)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_jitter(color = "grey30", size = 0.3, alpha = 0.4, width = 0.15) +
    scale_fill_manual(values = plot_colors) +
    scale_y_log10(labels = comma_format()) +
    labs(x = "Assembly Type", y = "Coverage (X)") +
    publication_theme +
    theme(legend.position = "none")

  # Replace Scatter Plot with Dumbbell Plot
  # Prepare data for dumbbell plot
  dumbbell_data <- data.frame(
    Length_Contig = contig_df$Length,
    Length_Scaffold = scaffold_df$Length,
    Coverage_Contig = contig_df$Coverage,
    Coverage_Scaffold = scaffold_df$Coverage
  ) %>%
    # Add row identifier
    mutate(ID = row_number())

  # Create long format for plotting
  length_dumbbell <- ggplot(dumbbell_data, aes(y = ID)) +
    geom_segment(aes(x = Length_Contig, xend = Length_Scaffold, 
                    y = ID, yend = ID),
                color = "grey70") +
    geom_point(aes(x = Length_Contig), color = plot_colors[2], size = 2) +
    geom_point(aes(x = Length_Scaffold), color = plot_colors[1], size = 2) +
    scale_x_log10(labels = comma_format()) +
    labs(x = "Contig Length (bp)", y = "") +
    publication_theme +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  coverage_dumbbell <- ggplot(dumbbell_data, aes(y = ID)) +
    geom_segment(aes(x = Coverage_Contig, xend = Coverage_Scaffold, 
                    y = ID, yend = ID),
                color = "grey70") +
    geom_point(aes(x = Coverage_Contig), color = plot_colors[2], size = 2) +
    geom_point(aes(x = Coverage_Scaffold), color = plot_colors[1], size = 2) +
    scale_x_log10(labels = comma_format()) +
    labs(x = "Coverage (X)", y = "") +
    publication_theme +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  # Create a legend for dumbbell plots
  legend_data <- data.frame(
    x = 1:2,
    y = c(1,1),
    Type = c("Contigs", "Scaffolds")
  )
  
  legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Type)) +
    geom_point(size = 3) +
    scale_color_manual(values = rev(plot_colors)) +
    theme_void() +
    theme(legend.position = "top")

  # Combine dumbbell plots
  dumbbell_combined <- plot_grid(
    get_legend(legend_plot),
    plot_grid(length_dumbbell, coverage_dumbbell, ncol = 2, align = 'h'),
    ncol = 1,
    rel_heights = c(0.1, 1)
  )

  # Delta Analysis
  scaffold_means <- colMeans(scaffold_df[c("Length", "Coverage")])
  contig_means <- colMeans(contig_df[c("Length", "Coverage")])
  percent_changes <- ((scaffold_means - contig_means) / contig_means) * 100
  
  delta_plot <- ggplot(data.frame(
    Metric = c("Length", "Coverage"),
    Change = percent_changes
  )) +
    geom_bar(aes(x = Metric, y = Change, fill = Metric), stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = plot_colors) +
    labs(y = "Percent Change (%)") +
    publication_theme +
    theme(legend.position = "none")

  # Combine all plots
  combined_plot <- plot_grid(
    plot_grid(length_plot, coverage_plot, 
             labels = c("A", "B"), 
             ncol = 2),
    dumbbell_combined,
    delta_plot,
    labels = c("", "C", "D"),
    ncol = 1,
    rel_heights = c(1, 1.2, 1)
  )

  # Add title
  title <- ggdraw() + 
    draw_label(paste("Assembly QC Metrics:", sample_name),
              fontface = 'bold',
              size = 16)

  final_plot <- plot_grid(
    title, combined_plot,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )

  # Save combined plot only
  ggsave(
    filename = file.path(output_dir, paste0(sample_name, "_qc_summary.png")),
    plot = final_plot,
    width = 12,
    height = 20,
    dpi = 600
  )
}

# Main function to process files
process_files <- function(scaffolds_dir, contigs_dir, output_dir) {
  tryCatch({
    scaffolds_dir <- normalizePath(scaffolds_dir)
    contigs_dir <- normalizePath(contigs_dir)
    output_dir <- normalizePath(output_dir)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    scaffold_files <- list.files(scaffolds_dir, pattern = "contigs\\.fasta$", full.names = TRUE, recursive = TRUE)
    contig_files <- list.files(contigs_dir, pattern = "_final_filtered_contigs\\.fa$", full.names = TRUE, recursive = TRUE)

    if (length(scaffold_files) == 0) {
      stop("No scaffold files found in the specified directory")
    }
    if (length(contig_files) == 0) {
      stop("No contig files found in the specified directory")
    }

    for (scaffold_file in scaffold_files) {
      tryCatch({
        sample_name <- basename(dirname(scaffold_file))
        contig_file <- contig_files[grepl(sample_name, contig_files, fixed = TRUE)]

        if (length(contig_file) == 0) {
          warning(paste("No matching contig file found for sample:", sample_name))
          next
        }

        message(paste("Processing sample:", sample_name))

        scaffold_df <- parse_fasta_headers(scaffold_file)
        contig_df <- parse_fasta_headers(contig_file[1])

        if (nrow(scaffold_df) == 0 || nrow(contig_df) == 0) {
          warning(paste("No valid data found for sample:", sample_name))
          next
        }

        create_qc_plot(scaffold_df, contig_df, output_dir, sample_name)
      }, error = function(e) {
        warning(paste("Error processing sample", sample_name, ":", e$message))
      })
    }
  }, error = function(e) {
    stop(paste("Error in process_files:", e$message))
  })
}


# For interactive testing, uncomment and modify these paths:
scaffolds_dir <- "scaffolds"
contigs_dir <- "filtered_contigs"
output_dir <- "contig-extender_QC"

process_files(scaffolds_dir, contigs_dir, output_dir)

