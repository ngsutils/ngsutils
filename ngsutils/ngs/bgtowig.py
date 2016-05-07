#!/usr/bin/env python
## category Conversion
## desc Converts bedGraph to Wiggle format
'''
Converts bedGraph format files to Wiggle (.wig) files. Writes the .wig file
to stdout.
'''

import os
import sys


def usage():
    print __doc__
    print "Usage: ngsutils bgtowig filename.bg"
    sys.exit(1)


def convert_bg_to_wig(fname):
    with open(fname) as f:
        cur_chrom = None
        for line in f:
            cols = line.strip().split('\t')
            if cols[0] != cur_chrom:
                cur_chrom = cols[0]
                print 'variableStep chrom=%s' % cur_chrom

            start = int(cols[1])
            end = int(cols[2])

            for i in xrange(start, end + 1):
                print '%s %s' % (i, cols[3])


if __name__ == '__main__':
    fname = None
    for arg in sys.argv[1:]:
        if not fname and os.path.exists(arg):
            fname = arg

    if not fname:
        usage()

    convert_bg_to_wig(fname)
