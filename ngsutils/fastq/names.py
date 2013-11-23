#!/usr/bin/env python
## category General
## desc Write out the read names
'''
Writes out the read names present in the FASTQ file.
'''

import os
import sys

from ngsutils.fastq import FASTQ


def export_names(fastq, include_comment=False, out=sys.stdout, quiet=False):
    for read in fastq.fetch(quiet=quiet):
        if include_comment:
            out.write('%s%s%s\n' % (read.name, ' ' if read.comment else '', read.comment))
        else:
            out.write('%s\n' % read.name)


def usage():
    print __doc__
    print "Usage: fastqutils names {-comment} filename.fastq{.gz}"
    sys.exit(1)

if __name__ == '__main__':
    fname = None
    include_comment = False

    for arg in sys.argv[1:]:
        if arg == '-comment':
            include_comment = True
        elif os.path.exists(arg):
            fname = arg

    if not fname:
        usage()

    fq = FASTQ(fname)
    export_names(fq, include_comment)
    fq.close()
