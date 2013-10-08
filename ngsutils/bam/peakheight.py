#!/usr/bin/env python
## category General
## desc Find the size (max height, width) of given peaks (BED) in a BAM file
'''
Find the size of a set of defined peaks (BED) in a BED file.

If you have a set of defined peak regions (in the form of a BED file), this
will return the size characteristics of the region for a given BAM file.

The following fields will be added to the BED file:
    maximum peak size
    total unique reads
    total coverage (sum of coverage over each base)
    coverage density (total coverage / length)
'''
import sys
import os
from ngsutils.bam import bam_open


def bam_peakheight(bam, bed_fobj, mask=1796):
    for line in bed_fobj:
        cols = line.strip('\n').split('\t')
        coverage_acc = 0
        max_coverage = 0
        reads = set()
        ref = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        for pileup in bam.pileup(reference=ref, start=start, end=end, mask=mask):
            coverage_acc += pileup.n
            if pileup.n > max_coverage:
                max_coverage = pileup.n

            for read in pileup.pileups:
                reads.add(read.alignment.qname)

        cols.append(max_coverage)
        cols.append(len(reads))
        cols.append(coverage_acc)
        cols.append(float(coverage_acc) / (end-start))
        print '\t'.join([str(x) for x in cols])


def usage():
    print __doc__
    print """\
Usage: bamutils peakheight {options} bamfile peaks.bed
"""
    sys.exit(1)

if __name__ == "__main__":
    bam_fname = None
    bed_fname = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif not bam_fname and os.path.exists(arg):
            bam_fname = arg
        elif not bed_fname and os.path.exists(arg):
            bed_fname = arg
        else:
            print 'Unknown argument: %s' % arg
            usage()
    if not bam_fname or not bed_fname:
        usage()

    bam = bam_open(bam_fname)
    with open(bed_fname) as f:
        bam_peakheight(bam, f)

    bam.close()
