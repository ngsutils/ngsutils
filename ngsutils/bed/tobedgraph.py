#!/usr/bin/env python
## category Conversion
## desc BED to BedGraph
'''
Takes a BED file with overlapping regions and produces a BedGraph file.  This
can optionally normalize the counts by a given factor.

 See: http://genome.ucsc.edu/goldenPath/help/bedgraph.html
      http://genome.ucsc.edu/goldenPath/help/bigWig.html
'''

import sys
import os
from ngsutils.support.eta import ETA


def write_regions(regions, normalize=None):
    while regions:
        regions.sort()
        fragment = None

        for chrom, start, end in regions:
            if start != regions[0][1]:
                fragment = (regions[0][1], start)
                break

        if not fragment:
            for chrom, start, end in regions:
                fragment = (regions[0][1], end)
                break

        new_regions = []
        count = 0
        for chrom, start, end in regions:
            if fragment[0] <= start < fragment[1]:
                count += 1
                if end > fragment[1]:
                    new_regions.append((chrom, fragment[1], end))
            else:
                new_regions.append((chrom, start, end))

        if normalize:
            sys.stdout.write('%s\t%s\t%s\t%s\n' % (chrom, fragment[0], fragment[1], int(normalize * count)))
        else:
            sys.stdout.write('%s\t%s\t%s\t%s\n' % (chrom, fragment[0], fragment[1], count))

        regions = new_regions

    pass


def bed_tobedgraph(fname, only_strand=None, normalize=None):
    with open(fname) as f:
        eta = ETA(os.stat(fname).st_size, fileobj=f)

        regions = []
        last_end = 0
        last_chrom = None

        for line in f:
            chrom, start, end, name, score, strand = line.strip().split('\t')
            eta.print_status(extra='%s:%s' % (chrom, start))

            if only_strand and only_strand != strand:
                continue

            if last_chrom != chrom or start > last_end:
                write_regions(regions, normalize)
                last_chrom = chrom
                regions = []

            regions.append((chrom, start, end))
            if end > last_end:
                last_end = end

        if regions:
            write_regions(regions, normalize)
    eta.done()


def usage():
    print __doc__
    print """\
Usage: bedutils tobedgraph [-plus | -minus] {-norm N} bamfile

Options:
    -plus             only count reads on the plus strand
                      (default: count all reads)
    -minus            only count reads on the minus strand

    -norm VAL         the count at every position is calculated as:
                      floor(count * VAL).
"""
    sys.exit(1)

if __name__ == "__main__":
    bed = None
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
        elif not bed and os.path.exists(arg):
            bed = arg
        else:
            print "Unknown option or missing index: %s" % arg
            usage()

    if not bed:
        usage()

    bed_tobedgraph(bed, strand, norm)
