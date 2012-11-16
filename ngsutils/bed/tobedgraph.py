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
from ngsutils.bed import BedFile


def write_regions(regions, normalize=None, out=sys.stdout):
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
            out.write('%s\t%s\t%s\t%s\n' % (chrom, fragment[0], fragment[1], int(normalize * count)))
        else:
            out.write('%s\t%s\t%s\t%s\n' % (chrom, fragment[0], fragment[1], count))

        regions = new_regions


def bed_tobedgraph(bed, only_strand=None, normalize=None, out=sys.stdout):
    regions = []
    last_end = 0
    last_chrom = None

    for region in bed:
        if only_strand and only_strand != region.strand:
            continue

        if last_chrom != region.chrom or region.start > last_end:
            write_regions(regions, normalize, out)
            last_chrom = region.chrom
            regions = []

        regions.append((region.chrom, region.start, region.end))
        if region.end > last_end:
            last_end = region.end

    if regions:
        write_regions(regions, normalize, out)


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

    bed_tobedgraph(BedFile(bed), strand, norm)
