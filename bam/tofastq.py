#!/usr/bin/env python
## category Conversion
## desc Convert BAM reads back to FASTQ sequences
'''
Convert BAM reads back to FASTQ sequences (mapped or unmapped)
'''

import sys
import os
from support.eta import ETA
import pysam


def bam_tofastq(sam_fname, colorspace=False, only_mapped=False, only_unmapped=False):
    if only_mapped == False and only_unmapped == False:
        return

    if sam_fname[-4:].lower() == '.bam':
        sam = pysam.Samfile(sam_fname, 'rb')
    else:
        sam = pysam.Samfile(sam_fname, 'r')

    eta = ETA(0, bamfile=sam)
    for read in sam:
        eta.print_status(extra=read.qname, bam_pos=(read.rname, read.pos))
        if only_mapped and read.is_unmapped:
            continue
        if only_unmapped and not read.is_unmapped:
            continue

        if colorspace:
            seq = read.opt('CS')
            qual = read.opt('CQ')
        else:
            seq = read.seq
            qual = read.qual

        sys.stdout.write('@%s\n%s\n+\n%s\n' % (read.qname, seq, qual))

    eta.done()


def usage():
    print __doc__
    print "Usage: bamutils tofastq [-cs] {-mapped} {-unmapped} file.[sam|bam]"
    sys.exit(1)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        usage()

    cs = False
    samf = None
    mapped = False
    unmapped = False

    for arg in sys.argv[1:]:
        if arg == '-cs':
            cs = True
        elif arg == '-unmapped':
            unmapped = True
        elif arg == '-mapped':
            mapped = True
        elif os.path.exists(arg):
            samf = arg
    if not samf:
        usage()

    if not unmapped and not mapped:
        mapped = True
        unmapped = True

    bam_tofastq(samf, cs, mapped, unmapped)
