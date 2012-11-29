#!/usr/bin/env python
## category Conversion
## desc Converts qual values from Illumina to Sanger scale
'''
Converts Illumina Qual values to Sanger scale.
'''

import os
import sys

from ngsutils.fastq import FASTQ, convert_illumina_qual


def usage():
    print __doc__
    print "Usage: fastqutils convertqual filename.fastq{.gz}"
    sys.exit(1)


def fastq_convertqual(fastq, out=sys.stdout, quiet=False):
    if fastq.check_qualtype() != "Illumina" and not quiet:
        sys.stderr.write("\nWarning: Unable to verify that FASTQ file contains Illumia scaled quality values!\n\n")

    for read in fastq.fetch(quiet=quiet):
        out.write('@%s\n%s\n+\n%s\n' % (read.name, read.seq, convert_illumina_qual(read.qual)))

if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if os.path.exists(arg):
            fname = arg

    if not fname:
        usage()

    fq = FASTQ(fname)
    fastq_convertqual(fq)
    fq.close()
