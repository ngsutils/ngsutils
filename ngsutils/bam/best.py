#!/usr/bin/env python
## category General
## desc Filter out multiple mappings for a read, selecting only the best
"""
Given a BAM file sorted by read name with potentially multiple mappings, this
will remove all but the best mapping.

The value of the attribute/tag given will be used to determine which reads
should be kept and which should be discarded. The tag should be a numeric
(int/float) type. Multiple tags may be used. This defaults to 'AS+, NM-'.
"""

import os
import sys
import pysam
import ngsutils.bam


def usage(msg=None):
    if msg:
        print msg
    print __doc__
    print """
Usage: bamutils best {opts} input.bam output.bam

Options
  -tag VAL    Tag to use to determine from which file reads will be taken.
              (must be type :i or :f) You may have more than one of these,
              in which case they will be sorted in order. You can add a +/-
              at the end of the name to signify sort order (asc/desc). 
              [default: AS+, NM-]

  -fail filename.bam    Write all failed mappings to this file.

"""
    sys.exit(1)


def bam_best(infile, outfile, failfile=None, tags=['AS+', 'NM-'], quiet=False):
    inbam = pysam.Samfile(infile, "rb")
    outbam = pysam.Samfile('%s.tmp' % outfile, "wb", template=inbam)

    if failfile:
        failbam = pysam.Samfile('%s.tmp' % failfile, "wb", template=inbam)

    for reads in ngsutils.bam.bam_batch_reads(inbam, quiet=quiet):
        best_val = None
        best_reads = []
        failed = []

        for read in reads:
            if not read.is_unmapped:
                tag_val = []
                for tag in tags:
                    val = float(read.opt(tag[:2]))
                    if tag[-1] == '-':
                        val = -val
                    tag_val.append(val)

                if not best_val or tag_val > best_val:
                    if best_reads and failfile:
                        failed.extend(best_reads)

                    best_val = tag_val
                    best_reads = [read]

                elif tag_val == best_val:
                    best_reads.append(read)

                elif failfile:
                    failed.append(read)

        for read in best_reads:
            outbam.write(read)

        if failfile:
            for read in failed:
                failbam.write(read)

    outbam.close()
    os.rename('%s.tmp' % outfile, outfile)

    if failfile:
        failbam.close()
        os.rename('%s.tmp' % failfile, failfile)


if __name__ == '__main__':
    infile = None
    outfile = None
    failfile = None
    last = None
    tags = []

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-tag':
            tags.append(arg)
            last = None
        elif last == '-fail':
            failfile = arg
            last = None
        elif arg in ['-tag', '-fail']:
            last = arg
        elif not infile and os.path.exists(arg):
            infile = arg
        elif not outfile:
            outfile = arg
        else:
            usage('Unknown option: %s' % arg)

    if not tags:
        tags = ['AS+', 'NM-']

    if not infile or not outfile:
        usage()
    else:
        bam_best(infile, outfile, failfile, tags)
