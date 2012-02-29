#!/usr/bin/env python
'''
Reduces BED regions to overlap them. The BED file *must* be sorted in order
to merge them.
'''

import os
import sys


def usage():
    print __doc__
    print """\
Usage: bedutils reduce {opts} bedfile

-extend num{,num}   Extend the BED region {num} bases 5' and 3'
                    Can either a single number (extend the same in both
                    direction) or a comma delimited pair of numbers where the
                    first number extends the region in the 5' direction and
                    the second number extends the region in the 3' direction.

-c                  Output the number of regions merged as the score (count).
                    Otherwise, the scores for all of the regions are added
                    together.

-nostrand           Ignore strand information when merging regions
"""


class MergedRegion(object):
    '''
    Manages regions to be merged together.
    '''
    def __init__(self, extend=(0, 0), strand='+'):
        '''
        extend is a tuple or list. The first is the 5' extension,
        the last is the 3' extension. These are strand specific.
        '''
        self.extend = extend
        self.strand = strand
        self._reset()

    def _reset(self):
        self.chrom = None
        self.start = 0
        self.end = 0
        self.score = 0
        self.members = []

    def add(self, chrom, start, end, name, score, strand):
        if strand == '+':
            start = start - self.extend[0]
            end = end + self.extend[1]
        else:
            start = start - self.extend[1]
            end = end + self.extend[0]

        if start < 0:
            start = 0

        if self.chrom != chrom:
            self.write()
            self._reset()
        elif start > self.end:
            self.write()
            self._reset()

        if self.start == 0:
            self.start = start

        self.chrom = chrom
        self.start = min(self.start, start)
        self.end = max(self.end, end)
        self.members.append(name)
        self.score += score

    def write(self):
        if self.chrom and self.start and self.end:
            sys.stdout.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (self.chrom, self.start, self.end, ','.join(self.members), self.score, self.strand))


def bed_reduce(fname, extend=(0, 0), stranded=True, count=False):
    if fname == '-':
        f = sys.stdin
    else:
        f = open(fname)

    plus_region = MergedRegion(extend, '+')
    minus_region = MergedRegion(extend, '-')

    # these are just for checking that the file is sorted
    lchrom = None
    lstart = None
    lend = None

    for line in f:
        chrom, start, end, name, score, strand = line.strip().split('\t')
        start = int(start)
        end = int(end)
        score = int(score)

        if lchrom == chrom:
            if start < lstart or (start == lstart and end < lend):
                sys.stderr.write('BED file is not sorted!\n')
                sys.stderr.write('chrom: %s\t%s (= %s)\n' % (lchrom, chrom, (chrom == lchrom)))
                sys.stderr.write('start: %s\t%s (< %s)\n' % (lstart, start, (lstart < start)))
                sys.stderr.write('end: %s\t%s\n' % (lend, end))
                sys.exit(1)

        lchrom = chrom
        lstart = start
        lend = end

        if not stranded or strand == '+':
            plus_region.add(chrom, start, end, name, score, strand)
        else:
            minus_region.add(chrom, start, end, name, score, strand)

    plus_region.write()
    minus_region.write()

    if f != sys.stdin:
        f.close()

if __name__ == '__main__':
    fname = None
    extend = (0, 0)
    stranded = True
    count = False
    last = None

    for arg in sys.argv[1:]:
        if last == '-extend':
            if ',' in arg:
                extend = [int(x) for x in arg.split(',')]
            else:
                extend = [int(arg), ] * 2
            last = None
        elif arg in ['-extend']:
            last = arg
        elif arg == '-nostrand':
            stranded = False
        elif arg == '-c':
            count = True
        elif not fname and (arg == '-' or os.path.exists(arg)):
            fname = arg
        else:
            print "Unknown option: %s" % arg
            usage()
            sys.exit(1)

    if not fname:
        usage()
        sys.exit(1)

    bed_reduce(fname, extend, stranded, count)
