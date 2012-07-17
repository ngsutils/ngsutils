#!/usr/bin/env python
## category Conversion
## desc Converts to FASTA format (seq or qual)
'''
Convert FASTQ to FASTA.  Optionally outputs just the quality values.

Same format as SOLiD csfasta / qual files.
'''

import os
import sys

from fastq_utils import read_fastq


def export_seq(fname):
    for name, seq, quals in read_fastq(fname, quiet=False):
        sys.stdout.write('>%s\n%s\n' % (name[1:], seq))


def export_qual(fname):
    for name, seq, quals in read_fastq(fname, quiet=False):
        sys.stdout.write('>%s\n%s\n' % (name[1:], ' '.join([str(ord(x) - 33) for x in quals])))


def usage():
    print """Usage: fastqutils tofasta {-qual} filename.fastq{.gz}
Options:
  -qual    Export the quality values (space separated numbers)
"""
    sys.exit(1)

if __name__ == '__main__':
    qual = False
    fname = None

    for arg in sys.argv[1:]:
        if arg == '-qual':
            qual = True
        elif os.path.exists(arg):
            fname = arg

    if not fname:
        usage()

    if qual:
        export_qual(fname)
    else:
        export_seq(fname)
