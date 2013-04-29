#!/usr/bin/env python
## category General
## desc Splits a BAM file into smaller pieces
"""
Splits a BAM file into smaller pieces

Given a BAM file, this script will split it into smaller BAM files with a
limit on the number of reads included.

Or it will also split a BAM file into a separate BAM file for each reference
that is included.
"""

import os
import sys
from ngsutils.bam import bam_iter
import pysam


def usage():
    print __doc__
    print """
Usage: bamutils split {-n num | -ref} in.bam out_template_name

out_template_name will be the template for the smaller BAM files.  They will
be named "out_template_name.N.bam" where out_template_name is the given
argument and N is the file number.

Options:
    -n      The number of reads to include in sub-files
            (default: 1000000)

    -ref    Split by references
"""
    sys.exit(1)


def bam_split(infile, out_template, read_count=1000000, reference=False, quiet=False):
    bamfile = pysam.Samfile(infile, "rb")
    outfile = None

    file_count = 0

    count = 0
    fname = ""
    lastref = -1
    for read in bam_iter(bamfile):
        if not outfile or (not reference and count >= read_count) or (reference and lastref != read.tid):
            if outfile:
                outfile.close()
            file_count += 1
            count = 0
            if reference:
                if read.tid >= 0:
                    fname = '%s.%s.bam' % (out_template, bamfile.getrname(read.tid))
                else:
                    fname = None
            else:
                fname = '%s.%s.bam' % (out_template, file_count)

            if fname:
                outfile = pysam.Samfile(fname, "wb", template=bamfile)
            else:
                outfile = None

        if outfile:
            outfile.write(read)
            count += 1

        lastref = read.tid

    bamfile.close()
    if outfile:
        outfile.close()
    if not quiet:
        sys.stderr.write("Split into %s files" % (file_count))


if __name__ == '__main__':
    infile = None
    outfile = None
    num = 1000000
    reference = False
    last = None

    for arg in sys.argv[1:]:
        if last == '-n':
            num = int(arg)
            last = None
        elif arg == '-ref':
            reference = True
        elif arg == '-h':
                usage()
        elif arg in ['-n']:
            last = arg
        elif not infile:
            infile = arg
        elif not outfile:
            outfile = arg

    if not infile or not outfile:
        usage()
    else:
        bam_split(infile, outfile, num, reference=reference)
