#!/usr/bin/env python
## category General
## desc Merges overlapping BED regions
'''
Reduces BED regions to overlap them. The BED file *must* be sorted in order
to merge them.
'''

import os
import sys
from ngsutils.bed import BedStreamer, BedRegion


def usage():
    print __doc__
    print """\
Usage: bedutils reduce {opts} bedfile

-extend num{,num}   Extend the BED region {num} bases 5' and 3'
                    Can either a single number (extend the same in both
                    direction) or a comma delimited pair of numbers where the
                    first number extends the region in the 5' direction and
                    the second number extends the region in the 3' direction.

-clip               Only extend the reads to find overlapping regions, don't
                    extend the edges of the regions.

-c                  Output the number of regions merged as the score (count).
                    Otherwise, the scores for all of the regions are added
                    together.

-nostrand           Ignore strand information when merging regions
"""
    sys.exit(1)


class MergedRegion(object):
    '''
    Manages regions to be merged together.
    '''
    def __init__(self, extend=(0, 0), clip=False, count=False, out=sys.stdout):
        '''
        extend is a tuple or list. The first is the 5' extension,
        the last is the 3' extension. These are strand specific.
        '''
        self.extend = extend
        self.clip = clip
        self.count = count
        self.out = out

        self._reset()

    def _reset(self):
        self.chrom = None
        self.extended_start = 0
        self.extended_end = 0
        self.score = 0
        self.strand = None
        self.region_start = 0
        self.region_end = 0
        self.members = []

    def add(self, region):
        if not region.strand or region.strand == '+':
            newstart = region.start - self.extend[0]
            newend = region.end + self.extend[1]
        else:
            newstart = region.start - self.extend[1]
            newend = region.end + self.extend[0]

        if newstart < 0:
            newstart = 0

        if self.chrom != region.chrom:
            self.write()
            self._reset()
        elif newstart >= self.extended_end:
            self.write()
            self._reset()

        if not self.extended_start:
            self.extended_start = newstart

        if self.clip:
            if not self.region_start:
                self.region_start = region.start

            self.region_end = region.end
        else:
            if not self.region_start:
                self.region_start = newstart

            self.region_end = newend

        self.chrom = region.chrom

        if not self.strand:
            self.strand = region.strand
        elif self.strand != region.strand:
            self.strand = '+'

        self.extended_start = min(self.extended_start, newstart)
        self.extended_end = max(self.extended_end, newend)
        if region.name:
            self.members.append(region.name)

        if self.count:
            self.score += 1
        else:
            self.score += region.score

    def write(self):
        if self.chrom and self.region_start and self.region_end:
            region = BedRegion(self.chrom, self.region_start, self.region_end, ','.join(sorted(set(self.members))), self.score, self.strand)
            region.write(self.out)
        self.region_start = 0
        self.region_end = 0


def bed_reduce(bed, extend=(0, 0), stranded=True, count=False, clip=False, out=sys.stdout):
    plus_region = MergedRegion(extend, clip, count, out)
    minus_region = MergedRegion(extend, clip, count, out)

    # these are just for checking that the file is sorted
    lchrom = None
    lstart = None
    lend = None
    lregion = None

    for region in bed:
        if lchrom == region.chrom:
            if region.start < lstart or (region.start == lstart and region.end < lend):
                print 'last    ', lregion
                print 'current ', region
                
                sys.stderr.write('BED file is not sorted!\n')
                sys.stderr.write('chrom: %s\t%s (= %s)\n' % (lchrom, region.chrom, (region.chrom == lchrom)))
                sys.stderr.write('start: %s\t%s (< %s)\n' % (lstart, region.start, (lstart < region.start)))
                sys.stderr.write('end: %s\t%s\n' % (lend, region.end))
                sys.exit(1)

        lchrom = region.chrom
        lstart = region.start
        lend = region.end

        lregion = region

        if not stranded or region.strand == '+':
            plus_region.add(region)
        else:
            minus_region.add(region)

    plus_region.write()
    minus_region.write()

if __name__ == '__main__':
    fname = None
    extend = (0, 0)
    stranded = True
    count = False
    last = None
    clip = False

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if last == '-extend':
            if ',' in arg:
                extend = [int(x) for x in arg.split(',')]
            else:
                extend = [int(arg), ] * 2
            last = None
        elif arg in ['-extend']:
            last = arg
        elif arg == '-clip':
            clip = True
        elif arg == '-nostrand':
            stranded = False
        elif arg == '-c':
            count = True
        elif not fname and (arg == '-' or os.path.exists(arg)):
            fname = arg
        else:
            print "Unknown option: %s" % arg
            usage()

    if not fname:
        usage()

    bed_reduce(BedStreamer(fname), extend, stranded, count, clip)
