#!/usr/bin/env python
## category Conversion
## desc Convert BAM reads to BED regions
'''
Convert BAM reads to BED regions
'''
import sys
import os
from ngsutils.bam import bam_iter
import pysam


def bam_tobed(fname):
    bamfile = pysam.Samfile(fname, "rb")

    for read in bam_iter(bamfile):
        if not read.is_unmapped:
            print '%s\t%s\t%s\t%s\t0\t%s' % (bamfile.getrname(read.rname), read.pos, read.aend, read.qname, '-' if read.is_reverse else '+')

    bamfile.close()


def usage():
    print __doc__
    print """\
Usage: bamutils tobed bamfile

Ouputs the read positions of all mapped reads in BED6 format.
"""

if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[-1]):
        usage()
        sys.exit(1)

    bam_tobed(sys.argv[1])
