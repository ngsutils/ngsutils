#!/usr/bin/env python
## category Conversion
## desc Converts qual values from Illumina to Sanger scale
'''
Converts Illumina Qual values to Sanger scale.
'''

import os
import sys

from fastq_utils import read_fastq


def usage():
    print __doc__
    print "Usage: fastqutils convertqual filename.fastq{.gz}"
    sys.exit(1)


def fastq_convertqual(fname):
    for name, seq, qual in read_fastq(fname):
        sys.stdout.write('@%s\n%s\n+\n%s\n' % (name, seq, ''.join([chr(ord(q) - 31) for q in qual])))

if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if os.path.exists(arg):
            fname = arg

    if not fname:
        usage()

    fastq_convertqual(fname)
