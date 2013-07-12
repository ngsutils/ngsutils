#!/usr/bin/env python
## category Conversion
## desc Convert BAM reads back to FASTQ sequences
'''
Convert BAM reads back to FASTA/FASTQ sequences (mapped or unmapped)
'''

import sys
import os
from ngsutils.bam import bam_iter, bam_open


def bam_tofastx(fname, colorspace=False, show_mapped=True, show_unmapped=True, fastq=True):
    if show_mapped is False and show_unmapped is False:
        return

    sam = bam_open(fname)

    last_qname = None

    for read in bam_iter(sam):
        if last_qname == read.qname:
            continue

        show = False
        if show_mapped and not read.is_unmapped:
            show = True
        if show_unmapped and read.is_unmapped:
            show = True

        if not show:
            continue

        if fastq:
            write_fastq(read, colorspace=colorspace)
            last_qname = read.qname
        else:
            write_fasta(read, colorspace=colorspace)
            last_qname = read.qname



def write_fasta(read, out=sys.stdout, colorspace=False):
    if colorspace:
        seq = read.opt('CS')
    else:
        seq = read.seq

    out.write('>%s\n%s\n' % (read.qname, seq))


def write_fastq(read, out=sys.stdout, colorspace=False):
    if colorspace:
        seq = read.opt('CS')
        qual = read.opt('CQ')
    else:
        seq = read.seq
        qual = read.qual

    out.write('@%s\n%s\n+\n%s\n' % (read.qname, seq, qual))


def usage(fastq=True):
    print __doc__

    if fastq:
        print "Usage: bamutils tofastq {opts} file.bam"
    else:
        print "Usage: bamutils tofasta {opts} file.bam"

    print """
Options:
    -cs        Output color-space sequences
    -mapped    Only output mapped sequences
    -unmapped  Only output unmapped sequences
"""
    sys.exit(1)


def main(fastq=True):
    if len(sys.argv) == 1:
        usage()

    cs = False
    samf = None
    mapped = True
    unmapped = True

    for arg in sys.argv[1:]:
        if arg == '-cs':
            cs = True
        elif arg == '-unmapped':
            if not unmapped:
                usage()
            mapped = False
        elif arg == '-mapped':
            if not mapped:
                usage()
            unmapped = False
        elif os.path.exists(arg):
            samf = arg
    if not samf:
        usage()

    bam_tofastx(samf, cs, mapped, unmapped, fastq)

if __name__ == '__main__':
    main()
