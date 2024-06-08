#!/usr/bin/env python

import sys, os, re
from collections import defaultdict

# formatting of input file (tab-delimited):
# read_id, cell_barcode, UMI, gene_id, transcript_id, gene_symbol


def main():

    usage = "usage: {} read_assignments.tsv.CB_UMI_info.sorted.tsv > read_assignments.tsv.CB_UMI_info.sorted.w_reads_split.tsv".format(sys.argv[0])

    if len(sys.argv) < 2:
        exit(usage)


    prev_molecule = ""
    prev_reads_list = list()
    
    sorted_reads_info_filename = sys.argv[1]
    with open(sorted_reads_info_filename, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            read_name = vals[0]
            if read_name != prev_molecule:
                process_split_reads(prev_reads_list)
                prev_reads_list.clear()

            prev_reads_list.append(vals)
            prev_molecule = read_name

    # get last ones
    process_split_reads(prev_reads_list)

    sys.exit(0)


def process_split_reads(reads_list):

    if len(reads_list) < 1:
        return

    if len(reads_list) == 1:
        read = reads_list.pop()
        read.append("1.0")
        print("\t".join(read))
        return

    # divide read among genes
    gene_to_isoforms = defaultdict(list)
    for read_alignment in reads_list:
        gene_id = read_alignment[3]
        gene_to_isoforms[gene_id].append(read_alignment)

    num_genes = len(gene_to_isoforms)
    for isoform_list in gene_to_isoforms.values():
        num_isoforms = len(isoform_list)
        read_count_per_isoform = 1 / num_genes / num_isoforms
        for isoform in isoform_list:
            isoform.append(f"{read_count_per_isoform:.2f}")
            print("\t".join(isoform))

    return

        
    



if __name__=='__main__':
    main()

