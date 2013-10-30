#!/usr/bin/env python
## category General
## desc Extract the sizes of BED regions
'''
Extract the sizes of BED regions
'''

import os
import sys
from ngsutils.bed import BedFile


def usage():
    print __doc__
    print """\
Usage: bedutils sizes bedfile

"""
    sys.exit(1)


def bed_size(bed, out=sys.stdout):
    for region in bed:
        out.write('%s\n' % (region.end-region.start))

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

    bed_size(BedFile(fname))
