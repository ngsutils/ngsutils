#!/usr/bin/env python
## category Conversion
## desc Converts color-space FASTQ file to encoded FASTQ
'''
Convert color-space FASTQ file to encoded color-space.

This is required for some tools that don't natively support color-space
files, such as BWA. For these tools, the color-space 01234 must be converted
into an encoded ACGTN.

0 -> A
1 -> C
2 -> G
3 -> T
4,5,6,. -> N
'''

import os
import sys

from ngsutils.fastq import FASTQ

_encoding = {
'0': 'A',
'1': 'C',
'2': 'G',
'3': 'T',
'4': 'N',
'5': 'N',
'6': 'N',
'.': 'N',
}


def encoded_seq(seq):
    ret = []
    for s in seq.upper():
        ret.append(_encoding[s])
    return ''.join(ret)


def fastq_csencode(fastq, out=sys.stdout, quiet=False):
    for read in fastq.fetch(quiet=quiet):
        if read.seq[0] in 'ATCG':  # linker prefix
            # skip 2 bases, unable to determine the proper convertion for
            # the first colorspace base

            seq = read.seq[2:]
            qual = read.qual[1:]
        else:
            seq = seq[1:]
            qual = read.qual[1:]

        out.write('@%s\n%s\n+\n%s\n' % (read.name, encoded_seq(seq), qual))


def usage():
    print __doc__
    print "Usage: fastqutils csencode filename.fastq{.gz}"
    sys.exit(1)


if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if os.path.exists(arg):
            fname = arg

    if not fname:
        usage()

    fq = FASTQ(fname)
    fastq_csencode(fq)
    fq.close()
