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


def bam_tobed(fname, out=sys.stdout):
    bamfile = pysam.Samfile(fname, "rb")

    for read in bam_iter(bamfile):
        write_read(read, bamfile.getrname(read.rname), out)
    bamfile.close()


def write_read(read, chrom, out):
    if not read.is_unmapped:
        out.write('%s\t%s\t%s\t%s\t0\t%s\n' % (chrom, read.pos, read.aend, read.qname, '-' if read.is_reverse else '+'))


def usage():  # pragma: no cover
    print __doc__
    print """\
Usage: bamutils tobed bamfile

Ouputs the read positions of all mapped reads in BED6 format.
"""

if __name__ == "__main__":  # pragma: no cover
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[-1]):
        usage()
        sys.exit(1)

    bam_tobed(sys.argv[1])
