# Load necessary libraries
library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(vegan)

# ALPHA DIVERSITY PLOTTING (Ridgeline and Violin Plots)

# Read in all alpha diversity results from bootstraps
alpha_dir <- Sys.getenv("ALPHA_DIR", "alpha_diversity_bootstraps")
files <- list.files(path = alpha_dir, pattern = "*.tsv", full.names = TRUE)

# Initialize an empty data frame to store all results
alpha_diversity_data <- data.frame()

# Loop over each file and add data to the dataframe
for (file in files) {
  df <- read.table(file, header = TRUE, sep = "\t")
  df$File <- basename(file)  # Add column to indicate which bootstrap replicate
  alpha_diversity_data <- rbind(alpha_diversity_data, df)
}

# Reshape the data to long format for ggplot
alpha_long <- alpha_diversity_data %>%
  pivot_longer(cols = c("shannon", "simpson", "pielou_e", "fisher_alpha", "chao1", "berger_parker"),
               names_to = "Metric", values_to = "Value")

# Plot ridgeline plots for alpha diversity
ggplot(alpha_long, aes(x = Value, y = SampleID, fill = Metric)) +
  geom_density_ridges(alpha = 0.8) +
  scale_fill_viridis_d() +
  facet_wrap(~ Metric, scales = "free") +
  theme_minimal() +
  labs(title = "Alpha Diversity Ridgeline Plot Across Bootstraps", x = "Diversity Index Value", y = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot violin plots for alpha diversity
ggplot(alpha_long, aes(x = SampleID, y = Value, fill = Metric)) +
  geom_violin(trim = FALSE) +
  geom_jitter(shape = 16, position = position_jitter(0.1), color = "black") +
  facet_wrap(~ Metric, scales = "free") +
  theme_minimal() +
  labs(title = "Alpha Diversity Violin Plots Across Bootstraps", y = "Diversity Index Value", x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# BETA DIVERSITY PLOTTING (Heatmap with Bootstrap Support)

# Read in all beta diversity distance matrices
beta_dir <- Sys.getenv("BETA_DIR", "beta_diversity_bootstraps")
beta_files <- list.files(path = beta_dir, pattern = "*.tsv", full.names = TRUE)

# Initialize a list to store all distance matrices
all_matrices <- list()

# Loop over each file and store the distance matrices
for (file in beta_files) {
  matrix <- as.matrix(read.table(file, header = TRUE, sep = "\t", row.names = 1))
  all_matrices[[file]] <- matrix
}

# Compute the average distance matrix across bootstraps
mean_matrix <- Reduce("+", all_matrices) / length(all_matrices)

# Plot heatmap with support (higher values mean more consistent distances)
pheatmap(mean_matrix,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "blue"))(100),
         main = "Beta Diversity Bootstrap Support Heatmap",
         display_numbers = TRUE,  # Display numeric values on the heatmap
         fontsize_number = 10)

# OPTIONAL: NMDS Plot for Beta Diversity
nmds_result <- metaMDS(mean_matrix, distance = "bray")

# Plot NMDS results
nmds_df <- as.data.frame(nmds_result$points)
nmds_df$Sample <- rownames(nmds_df)

ggplot(nmds_df, aes(x = MDS1, y = MDS2, label = Sample)) +
  geom_point(size = 4) +
  geom_text(nudge_x = 0.1, nudge_y = 0.1) +
  theme_minimal() +
  labs(title = "NMDS of Beta Diversity", x = "NMDS1", y = "NMDS2")
