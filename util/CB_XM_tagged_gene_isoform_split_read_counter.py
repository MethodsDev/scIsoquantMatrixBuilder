#!/usr/bin/env python

import sys, os, re
from collections import defaultdict

# formatting of input file (tab-delimited):
"""
molecule/35609701  AAACCTGAGCGACGTA  CCTACATACTTTTCTTATAT  ENSG00000003756.17  ENST00000469838.5        RBM5                1.0
molecule/11036232  AAACCTGAGCGACGTA  TACTCTCGCATTTCTTATAT  ENSG00000004059.11  ENST00000000233.10       ARF5                1.0
molecule/1874156   AAACCTGAGCGACGTA  TACTCTCGCATTTCTTATAT  ENSG00000004059.11  ENST00000000233.10       ARF5                1.0
molecule/24811702  AAACCTGAGCGACGTA  TACTCTCGCATTTCTTATAT  ENSG00000004059.11  ENST00000000233.10       ARF5                1.0
molecule/26292636  AAACCTGAGCGACGTA  GTCCATCAGCTTTCTTATAT  ENSG00000005893.16  ENST00000371335.4        LAMP2               1.0
molecule/12121149  AAACCTGAGCGACGTA  CAACATACGGTTTCTTATAT  ENSG00000006744.19  ENST00000338034.9        ELAC2               1.0
molecule/26085177  AAACCTGAGCGACGTA  CAACATACGGTTTCTTATAT  ENSG00000006744.19  ENST00000338034.9        ELAC2               1.0
molecule/47586568  AAACCTGAGCGACGTA  CAACATACGGTTTCTTATAT  ENSG00000006744.19  ENST00000338034.9        ELAC2               1.0
molecule/1700206   AAACCTGAGCGACGTA  TTTGAGACACTTTCTTATAT  ENSG00000007168.14  transcript225.chr17.nic  ENSG00000007168.14  1.0
molecule/6953292   AAACCTGAGCGACGTA  TTTGAGACACTTTCTTATAT  ENSG00000007168.14  transcript225.chr17.nic  ENSG00000007168.14  1.0
"""


def main():

    usage = "usage: {} read_assignments.tsv.CB_UMI_info.sorted.w_reads_split.cell_gene_iso_sorted.tsv output_prefix".format(
        sys.argv[0]
    )

    if len(sys.argv) != 3:
        exit(usage)

    prev_cell_gene_tok = ""
    prev_reads_list = list()

    sorted_info_filename = sys.argv[1]
    output_prefix = sys.argv[2]

    ofhs = dict()
    ofhs["ofh_cell_gene_counts"] = open(output_prefix + ".cell_gene_counts.tsv", "wt")
    ofhs["ofh_cell_isoform_counts"] = open(
        output_prefix + ".cell_isoform_counts.tsv", "wt"
    )
    ofhs["ofh_cell_unambig_read_isoform_counts"] = open(
        output_prefix + ".cell_unambig_isoform_counts.tsv", "wt"
    )

    with open(sorted_info_filename, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")

            cell_barcode = vals[1]
            gene_id = vals[3]

            cell_gene_tok = cell_barcode + "^" + gene_id

            if cell_gene_tok != prev_cell_gene_tok:
                process_cell_gene_reads(prev_reads_list, ofhs)
                prev_reads_list.clear()

            prev_reads_list.append(vals)
            prev_cell_gene_tok = cell_gene_tok

    # get last ones
    process_cell_gene_reads(prev_reads_list, ofhs)

    for ofh in ofhs.values():
        ofh.close()

    sys.exit(0)


def process_cell_gene_reads(reads_list, ofhs):

    if len(reads_list) < 1:
        return

    gene = reads_list[0][3]

    gene_sym = reads_list[0][5]
    # update gene sym according to known case (Isoquant doesnt indicate it for the nnic types)
    for read in reads_list:
        gene_sym_test = reads_list[0][5]
        if gene_sym_test != gene:
            gene_sym = gene_sym_test
            break

    cell = reads_list[0][1]

    gene_sum_reads = 0

    iso_total_read_count = defaultdict(float)
    iso_total_unambig_count = defaultdict(float)

    for read in reads_list:
        (
            read_id,
            cell_barcode,
            UMI,
            gene_id,
            transcript_id,
            gene_symbol,
            frac_read_assigned,
        ) = read
        frac_read_assigned = float(frac_read_assigned)

        assert gene_id == gene
        if gene_id != gene_symbol:
            assert (
                gene_symbol == gene_sym
            ), "Error conflicting gene symbol info for gene_id {} with {} vs {}\nfor reads:{}".format(
                gene, gene_sym, gene_symbol, reads_list
            )

        gene_sum_reads += frac_read_assigned
        iso_total_read_count[transcript_id] += frac_read_assigned
        if frac_read_assigned == 1.0:
            iso_total_unambig_count[transcript_id] += 1.0

    # add to outputs
    gene_tok = "^".join([gene_sym, gene])
    print(
        "\t".join([cell, gene_tok, f"{gene_sum_reads:.2f}"]),
        file=ofhs["ofh_cell_gene_counts"],
    )

    for isoform, count in iso_total_read_count.items():
        isoform_tok = "^".join([gene_sym, gene_id, isoform])
        print(
            "\t".join([cell, isoform_tok, f"{count:.2f}"]),
            file=ofhs["ofh_cell_isoform_counts"],
        )

    for isoform, count in iso_total_unambig_count.items():
        isoform_tok = "^".join([gene_sym, gene_id, isoform])
        print(
            "\t".join([cell, isoform_tok, f"{count:.2f}"]),
            file=ofhs["ofh_cell_unambig_read_isoform_counts"],
        )

    return


if __name__ == "__main__":
    main()
