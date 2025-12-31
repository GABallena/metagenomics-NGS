# Load required packages
library(ape)
library(ggtree)
library(phytools)

# Read and process the phylogenetic tree
tree <- read.tree(snakemake@input[["tree"]])
clusters <- read.table(snakemake@input[["clusters"]], header=FALSE)

# Create tree plot with mobilome information
p <- ggtree(tree) + 
    geom_tiplab(size=2) +
    theme_tree2() +
    xlim(0, max(node.depth.edgelength(tree))*1.2)

# Save plot
ggsave(snakemake@output[["plot"]], p, width=12, height=8)