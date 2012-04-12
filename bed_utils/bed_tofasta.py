#!/usr/bin/env python

import sys
import os
from support.eta import ETA
import pysam


def bed_tofasta(fname, ref_fasta, min_size=50, stranded=True):
    if not os.path.exists('%s.fai' % ref_fasta):
        pysam.faidx(ref_fasta)

    fasta = pysam.Fastafile(ref_fasta)

    refs = set()
    with open('%s.fai' % ref_fasta) as f:
        for line in f:
            refs.add(line.split('\t')[0].strip())

    with open(fname) as f:
        fsize = os.stat(fname).st_size
        eta = ETA(fsize, fileobj=f, modulo=1000)

        for line in f:
            if line[0] == '#':
                continue
            eta.print_status()

            cols = line.strip().split('\t')

            ref = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            strand = cols[5]

            if end - start >= min_size and ref in refs:
                seq = fasta.fetch(ref, start, end)
                if stranded:
                    if strand == '-':
                        seq = revcomp(seq)
                    print '>%s:%d-%d[%s]\n%s' % (ref, start, end, strand, seq)
                else:
                    print '>%s:%d-%d\n%s' % (ref, start, end, seq)

        eta.done()
    fasta.close()

_compliments = {
'a': 't',
'A': 'T',
'c': 'g',
'C': 'G',
'g': 'c',
'G': 'C',
't': 'a',
'T': 'A',
'n': 'n',
'N': 'N'
}


def revcomp(seq):
    ret = []

    for s in seq:
        ret.append(_compliments[s])

    ret.reverse()
    return ''.join(ret)


def usage():
    print """\
Usage: bedutils tofasta {-min size} {-ns} bedfile ref.fasta

Outputs the sequences of each BED region to FASTA format.

Option:
-min  The minumum size of a region
-ns   Ignore the strand of a region (always return seq from the + strand)
"""

if __name__ == "__main__":

    min_size = 50
    bed = None
    ref = None
    stranded = True

    last = None
    for arg in sys.argv[1:]:
        if last == '-min':
            min_size = int(arg)
            last = None
        elif arg in ['-min']:
            last = arg
        elif arg == '-ns':
            stranded = False
        elif not bed and os.path.exists(arg):
            bed = arg
        elif not ref and os.path.exists(arg):
            ref = arg

    if not bed or not ref:
        usage()
        sys.exit(1)

    bed_tofasta(bed, ref, min_size=min_size, stranded=stranded)
