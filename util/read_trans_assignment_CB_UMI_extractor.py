#!/usr/bin/env python

import sys, os, re
import gzip


def main():

    usage = "\n\n\tusage: {} read_assignments.tsv.gz\n\n".format(sys.argv[0])

    if len(sys.argv) != 2:
        exit(usage)

    read_assignments_filename = sys.argv[1]

    with gzip.open(read_assignments_filename, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue

            line = line.rstrip()

            vals = line.split("\t")

            read_id = vals[0]
            gene_id = vals[4]
            trans_id = vals[3]
            gene_info = vals[8]

            if gene_id == ".":
                continue

            m = re.search(" ?CB=([^;]+); XM=([^;]+);", gene_info)
            if m is not None:
                cell_barcode = m.group(1)
                umi = m.group(2)

                print("\t".join([read_id, gene_id, trans_id, cell_barcode, umi]))
            else:
                raise RuntimeError(
                    "Error, couldn't extract cell barcode and umi from {}".format(
                        gene_info
                    )
                )

    exit(0)


if __name__ == "__main__":
    main()
