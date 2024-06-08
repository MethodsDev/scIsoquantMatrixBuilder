#!/usr/bin/env Rscript



suppressPackageStartupMessages(library("argparse"))

parser = ArgumentParser()
parser$add_argument("--dat", help="input data file", required=TRUE, nargs=1)
parser$add_argument("--output_prefix", help="output prefix", required=TRUE, nargs=1)

args = parser$parse_args()    
dat_filename = args$dat
output_prefix = args$output_prefix

library(tidyverse)
library(Matrix)


message("-reading in ", dat_filename)
data = read.csv(dat_filename, sep="\t", header=F)
colnames(data) = c('cell_barcode', 'feature_id', 'UMI_counts')

data = data %>% filter(feature_id != '.')

message("-building sparse matrix")
feature_ids = data %>% select(feature_id) %>% unique() %>% pull(feature_id)
cell_barcodes = data %>% select(cell_barcode) %>% unique() %>% pull(cell_barcode)

data$feature_id = factor(data$feature_id, levels=feature_ids)
data$cell_barcode = factor(data$cell_barcode, levels=cell_barcodes)

data$feature_id = as.numeric(data$feature_id)
data$cell_barcode = as.numeric(data$cell_barcode)

sparseM = sparseMatrix(j=data$cell_barcode, i=data$feature_id, x=data$UMI_counts)

outdirname = paste0(output_prefix, "-sc_matrix_from_isoquant")

if (! dir.exists(outdirname)) {
  dir.create(outdirname, recursive = TRUE)
}


Matrix::writeMM(obj =  sparseM, file =  paste0(outdirname, "/matrix.mtx"))
write(feature_ids, file = paste0(outdirname, "/features.tsv"))
write(cell_barcodes, file = paste0(outdirname, "/barcodes.tsv"))

system(paste0("gzip ", outdirname, "/*")) 

message("done")

quit(save = "no", status = 0, runLast = FALSE)



### Info for loading into Seurat:
#
#  library(Seurat)
#
#  data = Read10X(data.dir="sample_id.genes-sc_matrix_from_isoquant/",
#               gene.column = 1,
#               cell.column = 2,
#               unique.features = TRUE,
#               strip.suffix = FALSE)
#
#
#  seurat_obj <- CreateSeuratObject(counts = data, project = "sample_id", min.cells = 3, min.features = 200)
#
##########################

