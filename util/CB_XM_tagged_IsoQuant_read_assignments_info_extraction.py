#!/usr/bin/env python3

import sys, os, re
import csv
import gzip

def main():

    usage = "usage: {} read_assignments.tsv.gz > gene_trans_CB_UMI.tsv\n\n".format(sys.argv[0])
    if len(sys.argv) < 2:
        exit(usage)


    read_assignments_file = sys.argv[1]

    with gzip.open(read_assignments_file, "rt") as fh:
        # skip first two lines
        next(fh)
        next(fh)
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:

            read_name = row['#read_id']
            gene_id = row['gene_id']
            transcript_id = row['isoform_id']

            m = re.search("CB=([^;]+); XM=([^;]+);", row['additional_info'])
            if m is not None:
                CB = m.group(1)
                UMI = m.group(2)
            else:
                raise RuntimeError("cannot extract CB and UMI from additional_info of row: {}".format(row))

            print("\t".join([read_name, gene_id, transcript_id, CB, UMI]))
            
    
    sys.exit(0)
    

if __name__=='__main__':
    main()
