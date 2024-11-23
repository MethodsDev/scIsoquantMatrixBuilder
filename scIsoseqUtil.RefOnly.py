#!/usr/bin/env python3

import sys, os, re

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "pylib"])
)
from Pipeliner import Pipeliner, Command
import logging
import argparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="generate sparse matrices from sc-IsoQuant outputs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--sample_id", type=str, required=True, help="sample_id")
    parser.add_argument(
        "--ref_annot_gtf", type=str, required=True, help="reference annotation gtf"
    )
    parser.add_argument(
        "--transcript_read_assignments",
        type=str,
        required=True,
        help="isoquant read_assignments.tsv.gz file",
    )

    args = parser.parse_args()

    sample_id = args.sample_id

    transcript_read_assignments_filename = args.transcript_read_assignments
    ref_annot_gtf = args.ref_annot_gtf

    utildir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "util")

    checkpoints_dir = "__checkpts"
    pipeliner = Pipeliner(checkpoints_dir)

    intermediates_dir = f"__{sample_id}.scIsoquant.tmpdir"
    if not os.path.exists(intermediates_dir):
        os.makedirs(intermediates_dir)

    # extract transcript, CB and UMI from bam
    read_CB_UMIs_filename = os.path.join(intermediates_dir, "read_trans_CB_umi.tsv")
    cmdstr = " ".join(
        [
            os.path.join(utildir, "read_trans_assignment_CB_UMI_extractor.py"),
            transcript_read_assignments_filename,
            f"| sort -u > {read_CB_UMIs_filename}",
        ]
    )
    cmdstr = f'/bin/bash -c "set -eof pipefail; {cmdstr}"'

    pipeliner.add_commands([Command(cmdstr, "read_trans_CB_UMIs.ok")])

    # get gene/trans mappings and gene name info:
    gene_trans_ids_filename = os.path.join(intermediates_dir, "gene_trans_info.tsv")
    cmdstr = " ".join(
        [
            os.path.join(utildir, "gtf_to_gene_trans_info.py"),
            ref_annot_gtf,
            f" | sort -u > {gene_trans_ids_filename}",
        ]
    )
    cmdstr = f'/bin/bash -c "set -eof pipefail; {cmdstr}"'

    pipeliner.add_commands([Command(cmdstr, "gene_trans_mappings.ok")])

    # merge read mappings with gene/trans and CB,UMI info
    read_CB_UMI_gene_trans_info_filename = os.path.join(
        intermediates_dir, "read_CB_UMI_gene_trans_info.tsv"
    )
    cmdstr = " ".join(
        [
            os.path.join(utildir, "merge_gene_trans_CB_UMI_counts.simpler.Rscript"),
            f"--read_trans_CB_UMIs {read_CB_UMIs_filename}",
            f"--gene_trans_map {gene_trans_ids_filename}",
            f"--output {read_CB_UMI_gene_trans_info_filename}",
        ]
    )

    pipeliner.add_commands([Command(cmdstr, "read_CB_UMI_gene_trans_info.ok")])

    # sort by read_name, gene
    read_info_sorted_by_readname_gene_filename = os.path.join(
        intermediates_dir, "read_info.sorted_rname_gene.tsv"
    )
    # read_id, cell_barcode, UMI, gene_id, transcript_id, gene_symbol
    cmdstr = f"cat {read_CB_UMI_gene_trans_info_filename} | sort -k1,1 -k4,4 >  {read_info_sorted_by_readname_gene_filename}"
    cmdstr = f'/bin/bash -c "set -eof pipefail; {cmdstr}"'

    pipeliner.add_commands([Command(cmdstr, "read_info_sorted_rname_gene.ok")])

    # split reads across genes and isoforms of genes
    reads_split_filename = os.path.join(
        intermediates_dir, "read_info.w_reads_split.tsv"
    )
    cmdstr = " ".join(
        [
            os.path.join(utildir, "CB_XM_tagged_gene_isoform_read_splitter.py"),
            read_info_sorted_by_readname_gene_filename,
            ">",
            reads_split_filename,
        ]
    )

    pipeliner.add_commands([Command(cmdstr, "reads_split.ok")])

    # sort by cell, gene, isoform:
    sorted_split_reads_filename = reads_split_filename + ".sorted"
    cmdstr = " ".join(
        [
            "sort -k2,2 -k4,4 -k5,5 ",
            reads_split_filename,
            ">",
            sorted_split_reads_filename,
        ]
    )
    pipeliner.add_commands([Command(cmdstr, "sorted_split_reads.ok")])

    # sum the gene and isoform UMIs
    output_prefix = f"{sample_id}.read_assignments"
    cmdstr = " ".join(
        [
            os.path.join(utildir, "CB_XM_tagged_gene_isoform_split_read_counter.py"),
            sorted_split_reads_filename,
            output_prefix,
        ]
    )

    pipeliner.add_commands([Command(cmdstr, "gene_isoform_reads_counted.ok")])

    # the above will generate three counts.tsv files: cell_gene_counts, cell_isofomr_counts, and cell_unambig_isoform_counts
    # where the cell_unambig_isoform_counts only counts those reads unambiguously assigned to specific isoforms.

    # finally, build the sparse matrices for use with Seurat - do this for each of the above counts files:
    # with the 2nd parameter below being the output prefix.

    cmdstr = " ".join(
        [
            os.path.join(
                utildir, "CB_XM_tagged_IsoQuant_read_assignments_to_sparse_matrix.R"
            ),
            f"--dat {output_prefix}.cell_gene_counts.tsv",
            f"--output_prefix {sample_id}.genes",
        ]
    )
    pipeliner.add_commands([Command(cmdstr, "genes.sparse_matrix.ok")])

    cmdstr = " ".join(
        [
            os.path.join(
                utildir, "CB_XM_tagged_IsoQuant_read_assignments_to_sparse_matrix.R"
            ),
            f"--dat {output_prefix}.cell_isoform_counts.tsv",
            f"--output_prefix {sample_id}.isoforms",
        ]
    )

    pipeliner.add_commands([Command(cmdstr, "isoforms.sparse_matrix.ok")])

    cmdstr = " ".join(
        [
            os.path.join(
                utildir, "CB_XM_tagged_IsoQuant_read_assignments_to_sparse_matrix.R"
            ),
            f"--dat {output_prefix}.cell_unambig_isoform_counts.tsv",
            f"--output_prefix {sample_id}.unambig_isoforms",
        ]
    )

    pipeliner.add_commands([Command(cmdstr, "unambiguous_isoforms.sparse_matrix.ok")])

    pipeliner.run()

    sys.exit(0)


if __name__ == "__main__":
    main()
