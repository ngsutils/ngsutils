#!/usr/bin/env python
## category General
## desc Extends BED regions (3')
'''
Extends BED regions (3' only)
'''

import os
import sys
from ngsutils.bed import BedFile, BedRegion


def usage():
    print __doc__
    print """\
Usage: bedutils extend {+}SIZE bedfile

SIZE is what the total size of the region should be.  The size of the region
will be extended or reduced to make the total length exactly SIZE. The region
is always adjusted at the 3' end, respective to the given strand.

If the first character of SIZE is '+', then the region is extended
SIZE bases, regardless of how long it is to start with.

"""
    sys.exit(1)


def bed_extend(bed, size, relative=False, out=sys.stdout):
    for region in bed:
        if region.strand == '+':
            start = region.start
            if relative:
                end = region.end + size
            else:
                end = region.start + size
        else:
            end = region.end
            if relative:
                start = region.start - size
            else:
                start = region.end - size

        if start < 0:
            start = 0

        out.write('%s\n' % BedRegion(region.chrom, start, end, region.name, region.score, region.strand))

if __name__ == '__main__':
    fname = None
    size = None
    relative = False
    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if not size:
            if arg[0] == '+':
                relative = True
                size = int(arg[1:])
            else:
                size = int(arg)
        elif not fname and os.path.exists(arg):
            fname = arg
        else:
            print "Unknown option: %s" % arg
            usage()

    if not fname or not size:
        usage()

    bed_extend(BedFile(fname), size, relative)
