#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import gzip

COLORMAP = {
    'A': [ 'A', 'C', 'G', 'T' ],
    'C': [ 'C', 'A', 'T', 'G' ],
    'G': [ 'G', 'T', 'A', 'C' ],
    'T': [ 'T', 'G', 'C', 'A' ]
}

def decode(seq):
    if len(seq) == 0: return ""
    last = seq[0]; retval = ""
    if last not in COLORMAP: return ""
    for i in range(1, len(seq)):
        try:
            last = COLORMAP[last][int(seq[i])] if last != 'N' else 'N'
        except:
            last = 'N'
        retval += last
    return retval

def qual(seq):
    retval = ""
    for x in seq.strip().split():
        retval += chr(int(x) + 32)
    return retval

def writeseq(name, seq, qual, o):
    o.write("@%s\n%s\n+\n%s\n" % (name.strip(), seq, qual))

def main(argc, argv):

    if argc < 4:
        print("usage: convert.py input.fasta input.csqual output.fastq", file = sys.stderr)
        return 1
    fastain = argv[1]
    qualin = argv[2]
    seqs = {}; quals = {}
    
    with (gzip.open if fastain.endswith(".gz") else open)(fastain, 'r') as fastafile:
        with (gzip.open if qualin.endswith(".gz") else open)(qualin, 'r') as qualfile:
            with gzip.open(argv[3] + (".gz" if not argv[3].endswith(".gz") else ""), 'w') as o:
                fname = fastafile.readline(); qname = qualfile.readline(); i = 0
                while fname and fname.startswith('#'):
                    fname = fastafile.readline()
                while qname and qname.startswith('#'):
                    qname = qualfile.readline()
                while fname and qname:
                    i += 1
                    if i % 100000 == 0: print(i)
                    if qname == fname:
                        writeseq(fname, decode(fastafile.readline()), qual(qualfile.readline()), o)
                    else:
                        if qname in seqs:
                            writeseq(qname, seqs[qname], qual(qualfile.readline()), o)
                            del seqs[qname]
                        else:
                            quals[qname] = qual(qualfile.readline())
                        if fname in quals:
                            writeseq(fname, decode(fastafile.readline()), quals[fname], o)
                            del quals[fname]
                        else:
                            seqs[fname] = decode(fastafile.readline())
                    qname = qualfile.readline()
                    fname = fastafile.readline()

    return 0

if __name__ == "__main__":
    sys.exit(main(len(sys.argv), sys.argv))
