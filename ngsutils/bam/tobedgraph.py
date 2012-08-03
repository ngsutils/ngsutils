#!/usr/bin/env python
## category Conversion
## desc Convert BAM coverage to bedGraph (for visualization)
'''
Convert BAM coverage to bedGraph

Takes a BAM file and produces a bedGraph file.  This can optionally
normalize the counts by a given factor.

 See: http://genome.ucsc.edu/goldenPath/help/bedgraph.html
      http://genome.ucsc.edu/goldenPath/help/bigWig.html
'''

import sys
import os
from ngsutils.support.eta import ETA
import pysam


def write_bedgraph(chrom, start, end, count, normalize=None):
    if start and end and count:
        if normalize:
            sys.stdout.write('%s\t%s\t%s\t%s\n' % (chrom, start, end, int(normalize * count)))
        else:
            sys.stdout.write('%s\t%s\t%s\t%s\n' % (chrom, start, end, count))


def bam_tobedgraph(bam_name, strand=None, normalize=None):
    bamfile = pysam.Samfile(bam_name, "rb")
    eta = ETA(0, bamfile=bamfile)

    last_chrom = None
    last_count = 0
    last_start = None
    last_end = None

    for pileup in bamfile.pileup():
        # sys.stdin.readline()
        chrom = bamfile.getrname(pileup.tid)

        if chrom != last_chrom and last_count > 0 and last_end:
            write_bedgraph(last_chrom, last_start, last_end + 1, last_count, normalize)
            last_count = 0
            last_start = None
            last_end = None

        eta.print_status(extra='%s:%s' % (chrom, pileup.pos), bam_pos=(pileup.tid, pileup.pos))

        count = 0
        if strand is None:
            count = pileup.n
        else:
            for read in pileup.pileups:
                if not read.is_reverse and strand == '+':
                    count += 1
                elif read.is_reverse and strand == '-':
                    count += 1

        # print pileup.pos,count,last_start,last_end
        if count != last_count or (pileup.pos - last_end) > 1:
            if last_count > 0:
                write_bedgraph(last_chrom, last_start, last_end + 1, last_count, normalize)

            if count == 0:
                last_start = None
            # elif not last_end or (pileup.pos-last_end) > 1:
            #     last_start = pileup.pos
            else:
                last_start = pileup.pos

        last_end = pileup.pos
        last_count = count
        last_chrom = chrom

    write_bedgraph(last_chrom, last_start, last_end + 1, last_count, normalize)

    eta.done()
    bamfile.close()


def usage():
    print __doc__
    print """\
Usage: bamutils tobedgraph [-plus | -minus] {-norm N} bamfile

Options:
    -plus             only count reads on the plus strand
                      (default: count all reads)
    -minus            only count reads on the minus strand

    -norm VAL         the count at every position is calculated as:
                      floor(count * VAL).
"""
    sys.exit(1)

if __name__ == "__main__":
    bam = None
    strand = None
    norm = None

    last = None
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if last == '-norm':
            norm = float(arg)
            last = None
        elif arg in ['-norm']:
            last = arg
        elif arg == '-plus':
            strand = '+'
        elif arg == '-minus':
            strand = '-'
        elif not bam and os.path.exists(arg) and os.path.exists('%s.bai' % arg):
            bam = arg
        else:
            print "Unknown option or missing index: %s" % arg
            usage()

    if not bam:
        usage()

    bam_tobedgraph(bam, strand, norm)
