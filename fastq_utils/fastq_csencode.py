#!/usr/bin/env python
## category Conversion
## desc Converts color-space FASTQ file to encoded FASTQ
'''
Convert color-space FASTQ file to encoded color-space.

This is required for some tools that don't natively support color-space
files, such as BWA. For these tools, the color-space 01234 must be converted
into an encoded ACGTN.
'''

import os
import sys

from fastq_utils import read_fastq

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
    for s in seq:
        ret.append(_encoding[s])
    return ''.join(ret)


def fastq_csencode(fname):
    for name, seq, qual in read_fastq(fname, quiet=False):
        if seq:
            if seq[0] in 'ATCG':
                seq = seq[2:]
            else:
                seq = seq[1:]

            sys.stdout.write('%s\n%s\n+\n%s\n' % (name, encoded_seq(seq), qual[1:]))


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

    fastq_csencode(fname)
