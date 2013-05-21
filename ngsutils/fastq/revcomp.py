#!/usr/bin/env python
## category General
## desc Reverse compliment a FASTQ file
'''
Reverse compliments each read in a FASTQ file.
(base-space only)
'''

import os
import sys

from ngsutils.fastq import FASTQ
import ngsutils.support


def fastq_revcomp(fastq, out=sys.stdout, quiet=False):
    if fastq.is_colorspace:
        sys.stderr.write("Reverse-complimenting only works on base-space files")
        sys.exit(1)

    for read in fastq.fetch(quiet=quiet):
        seq = ngsutils.support.revcomp(read.seq)
        qual = read.qual[::-1]

        read.clone(seq=seq, qual=qual).write(out)


def usage():
    print __doc__
    print "Usage: fastqutils revcomp filename.fastq{.gz}"
    sys.exit(1)

if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()

        if not fname and os.path.exists(arg):
            fname = arg

    if not fname:
        usage()

    fq = FASTQ(fname)
    fastq_revcomp(fq)
    fq.close()
