#!/usr/bin/env python
## category General
## desc Cleans a BED file (score should be integers)
'''
Cleanup a BED file, rounds scores to integers
'''

import os
import sys
from ngsutils.bed import BedFile, BedRegion


def usage():
    print __doc__
    print """\
Usage: bedutils clean bedfile

Converts the "score" field to be an integer

"""
    sys.exit(1)


def bed_clean(bed, out=sys.stdout):
    for region in bed:
        region.score = int(region.score)
        if region.score > 1000:
            region.score = 1000

        out.write('%s\n' % region)

if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg
        else:
            print "Unknown option: %s" % arg
            usage()

    if not fname:
        usage()

    bed_clean(BedFile(fname))
