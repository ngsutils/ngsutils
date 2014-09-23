#!/usr/bin/env python
## category Conversion
## desc Converts to BAM format (unmapped)
'''
Convert FASTQ to BAM. This doesn't perform any mapping, it simply stores the
read sequences in BAM format as unmapped reads. If given two files, the reads
will be correctly flagged as pairs.
'''

import os
import sys
import itertools
import pysam

from ngsutils.fastq import FASTQ


def export_bam(outbam, read1, read2, quiet=False):
    if read2:
        def gen():
            for r1, r2 in itertools.izip(read1.fetch(quiet=quiet), read2.fetch(quiet=True)):
                yield (r1, r2)
    else:
        def gen():
            for r1 in read1.fetch(quiet=quiet):
                yield (r1, None)


    for r1, r2 in gen():
        record1 = pysam.AlignedRead()
        record1.qname = r1.name
        record1.seq = r1.seq
        record1.qual = r1.qual

        if r2:
            record1.is_paired = True
            record1.is_read1 = True

            record2 = pysam.AlignedRead()
            record2.qname = r1.name
            record2.seq = r1.seq
            record2.qual = r1.qual
            record2.is_paired = True
            record2.is_read2 = True

        outbam.write(record1)
        if r2:
            outbam.write(record2)


def usage():
    print """Usage: fastqutils tobam {opts} outfile.bam read1.fastq{.gz} {read2.fastq}

Note: If two FASTQ files are given, they are assumed to be paired end reads.

Options:
  -f    Force overwriting output file

"""
    sys.exit(1)

if __name__ == '__main__':
    outname = None
    read1_fname = None
    read2_fname = None

    force = False

    for arg in sys.argv[1:]:
        if arg == '-f':
            force = True
        elif not outname:
            if not force and os.path.exists(arg):
                usage('Output file exists! (Use -f to force overwriting): %s' % arg)
            outname = arg
        elif not read1_fname and os.path.exists(arg):
            read1_fname = arg
        elif not read2_fname and os.path.exists(arg):
            read2_fname = arg

    if not outname or not read1_fname:
        usage()

    read1 = FASTQ(read1_fname)
    read2 = FASTQ(read2_fname) if read2_fname else None

    bam = pysam.Samfile(outname, 'wb')
    export_bam(bam, read1, read2)
    bam.close()

    read1.close()
    read2.close()
