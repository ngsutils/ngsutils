#!/usr/bin/env python
## category General
## desc Write out the read names
'''
Writes out the read names present in the FASTQ file.
'''

import os
import sys

from ngsutils.fastq import FASTQ


def export_names(fastq, out=sys.stdout, quiet=False):
    for read in fastq.fetch(quiet=quiet):
        out.write('%s\n' % read.name.split()[0])


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

    fq = FASTQ(fname)
    export_names(fq)
    fq.close()
