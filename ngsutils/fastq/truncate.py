#!/usr/bin/env python
## category General
## desc Truncates reads to a maximum length
'''
Truncates each read seq/qual values to be a maximum length
'''

import os
import sys

from ngsutils.fastq import FASTQ


def fastq_truncate(fastq, max_len, out=sys.stdout, quiet=False):
    for read in fastq.fetch(quiet=quiet):
        if fastq.is_colorspace and read.seq[0] in "atcgATGC":
            seq = read.seq[:max_len + 1]
            qual = read.qual[:max_len]
        else:
            seq = read.seq[:max_len]
            qual = read.qual[:max_len]

        read.clone(seq=seq, qual=qual).write(out)


def usage():
    print __doc__
    print "Usage: fastqutils truncate filename.fastq{.gz} max_len"
    sys.exit(1)

if __name__ == '__main__':
    fname = None
    maxlen = 0

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if not fname and os.path.exists(arg):
            fname = arg
        else:
            maxlen = int(arg)

    if not fname or not maxlen:
        usage()

    fq = FASTQ(fname)
    fastq_truncate(fq, maxlen)
    fq.close()
