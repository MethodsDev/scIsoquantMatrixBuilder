#!//usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("tidyverse"))

parser = ArgumentParser()
parser$add_argument("--sc_isoform_pseudobulk_matrix_for_diffusage", help="sc isoform pseudobulk count matrix with gene_id and transcript_id columns", required=TRUE, nargs=1)
parser$add_argument("--gene_isoform_exon_intron_patterns", help="file generated by gtf_to_isoform_exon_splice_patterns.py on gtf file", required=TRUE, nargs=1)
parser$add_argument("--output_matrix", help="name for output matrix file", required=TRUE, nargs=1)

args = parser$parse_args()
input_isoform_cluster_count_matrix = args$sc_isoform_pseudobulk_matrix_for_diffusage
gene_isoform_exon_intron_patterns_file = args$gene_isoform_exon_intron_patterns

input_counts_matrix = read.csv(input_isoform_cluster_count_matrix, header=T, sep="\t")

#input_counts_matrix$transcript_id = sapply(input_counts_matrix$transcript_id, function(x) { y = str_split(x, "\\^")[[1]]; y[length(y)] })

input_counts_matrix$transcript_id = sapply(input_counts_matrix$transcript_id, function(x) { str_replace_all(x, "-", "_") })

gene_isoform_exon_intron_patterns = read.csv(gene_isoform_exon_intron_patterns_file, header=T, sep="\t")
    


counts_matrix = left_join(input_counts_matrix,
                           gene_isoform_exon_intron_patterns %>% select(transcript_id, introns_hashcode),
                          by='transcript_id')

# sum across splice patterns.
counts_matrix = counts_matrix %>% mutate(introns_hashcode = ifelse(is.na(introns_hashcode) | introns_hashcode == "", transcript_id, introns_hashcode))


counts_matrix_sp_sums = counts_matrix %>% select(-transcript_id) %>% rename(transcript_id = introns_hashcode) %>%
    group_by(gene_id, transcript_id) %>% summarize(across(everything(), sum))
    
#counts_matrix_sp_sums = counts_matrix %>% select(gene_id, transcript_id, all_of(cluster_names))


write.table(counts_matrix_sp_sums, file=args$output_matrix, sep="\t", quote=F, row.names=FALSE)

quit(save = "no", status = 0, runLast = FALSE)

