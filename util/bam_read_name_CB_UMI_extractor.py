#!/usr/bin/env python

import sys, os, re
import pysam


def main():

    usage = "\n\n\tusage: {} reads.bam\n\n".format(sys.argv[0])

    if len(sys.argv) != 2:
        exit(usage)

    bam_path = sys.argv[1]
        
    samfile = pysam.AlignmentFile(bam_path, "rb", check_sq = False)
    for read in samfile:
        print("\t".join([read.query_name, read.get_tag("CB"), read.get_tag("XM")]))


    exit(0)


if __name__=='__main__':
    main()
    
        
 
