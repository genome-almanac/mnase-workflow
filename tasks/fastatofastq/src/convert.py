#!/usr/bin/env python

from __future__ import print_function

import sys
import os

from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

def main(argc, argv):

    if argc < 4:
        print("usage: convert.py input.fasta input.csqual output.fastq", file = sys.stderr)
        return 1

    fastain = argv[1]
    if argv[1].endswith(".gz"):
        os.system("gunzip %s" % argv[1])
        fastain = argv[1][:-3]
    qualin = argv[2]
    if argv[2].endswith(".gz"):
        os.system("gunzip %s" % argv[2])
        qualin = argv[2][:-3]
    
    with open(fastain, 'r') as fastafile:
        with open(qualin, 'r') as qualfile:
            rec_iter = PairedFastaQualIterator(fastafile, qualfile)
            with open(argv[3], 'w') as o:
                try:
                    SeqIO.write(rec_iter, o, "fastq")
                except:
                    print("warning: skipped a record which failed to write", file = sys.stderr)
                    pass

    return os.system("gzip %s" % argv[3])

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
