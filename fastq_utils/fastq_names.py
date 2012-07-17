#!/usr/bin/env python
## category General
## desc Write out the read names
'''
Writes out the read names present in the FASTQ file.
'''

import os
import sys

from fastq_utils import read_fastq


def usage():
    print __doc__
    print "Usage: fastqutils names filename.fastq{.gz}"
    sys.exit(1)

if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if os.path.exists(arg):
            fname = arg

    if not fname:
        usage()

    for name, seq, qual in read_fastq(fname):
        sys.stdout.write('%s\n' % name.split()[0][1:])
