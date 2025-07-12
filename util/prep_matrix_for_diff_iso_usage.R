#!/bin/env Rscript
library(tidyverse)
library(argparse)

# Create argument parser
parser <- ArgumentParser(description = "Process isoform matrix data for differential isoform usage analysis")

# Add arguments
parser$add_argument("-i", "--input", 
                    required = TRUE,
                    help = "Input file path (tab-separated matrix file)")

parser$add_argument("-o", "--output",
                    required = TRUE, 
                    help = "Output file path for processed matrix")

# Parse arguments
args <- parser$parse_args()

# Read input file
df <- read.csv(args$input, header = TRUE, row.names = 1, sep = "\t")

# Process data
cluster_names <- colnames(df)
df$transcript_id <- rownames(df)
df$gene_id <- sub("\\^[^\\^]+$", "", df$transcript_id)
df <- df %>% select(gene_id, transcript_id, all_of(cluster_names))

# Write output file
write.table(df, file = args$output, sep = "\t", quote = FALSE, row.names = FALSE)

# Print completion message
cat("Processing complete!\n")
cat("Input file:", args$input, "\n")
cat("Output file:", args$output, "\n")