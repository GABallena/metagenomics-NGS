# Load necessary libraries
library(pavian)
library(gridExtra)
library(ggplot2)

# Specify the input files (replace these with your actual file paths)
kraken2_summary_files <- list.files(path = "kraken2_output", pattern = "_kraken2_summary.txt", full.names = TRUE)

# Set the output PDF file
pdf("pavian_output.pdf")

# Generate the plots and tables
for (summary_file in kraken2_summary_files) {
  # Load the Kraken2 report
  kraken_report <- read_kraken(summary_file)
  
  # Generate the classification summary
  classification_summary <- pavian::report(kraken_report)
  
  # Plot the classification summary
  print(classification_summary)
  
  # Add more plots/tables as needed
  # Example: ggplot2 bar plot
  gg <- ggplot(kraken_report, aes(x = taxon, y = count)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        labs(title = "Taxonomic Classification", x = "Taxon", y = "Count")
  print(gg)
}

# Close the PDF device
dev.off()
