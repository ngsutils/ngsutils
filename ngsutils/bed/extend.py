#!/usr/bin/env python
## category General
## desc Extends BED regions (3')
'''
Extends BED regions (5' and 3')
'''

import os
import sys
from ngsutils.bed import BedStreamer


def usage(msg=None):
    if msg:
        print 'ERROR: %s' % msg
    print __doc__
    print """\
Usage: bedutils extend {opts} bedfile

Options:
    -5 SIZE     Extend the region SIZE bases in the 5' direction
    -3 SIZE     Extend the region SIZE bases in the 5' direction

SIZE is how much the region should be extended in the corresponding
direction. If the first character of SIZE is '=', then the region is
extended in the corresponding direction and the length of the region
is set to be exactly SIZE. Only one direction is allowed when SIZE
starts with "=".

Note: using an exact size may make a region smaller!

"""
    sys.exit(1)


def bed_extend(stream, ext_5=None, ext_3=None, exact=False, out=sys.stdout):
    for region in stream:
        if region.strand == '+' or not region.strand:
            if ext_5:
                if exact:
                    start = region.end - ext_5
                else:
                    start = region.start - ext_5

            elif ext_3:
                if exact:
                    end = region.start + ext_3
                else:
                    end = region.end + ext_3

        else:
            if ext_5:
                if exact:
                    end = region.start + ext_5
                else:
                    end = region.end + ext_5

            elif ext_3:
                if exact:
                    start = region.end - ext_3
                else:
                    start = region.start - ext_3

        if start < 0:
            start = 0

        out.write('%s\n' % region.clone(start=start, end=end))

if __name__ == '__main__':
    fname = None
    ext_5 = None
    ext_3 = None
    exact = False

    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-5':
            if arg[0] == '=':
                exact = True
                ext_5 = int(arg[1:])
            else:
                ext_5 = int(arg)
            last = None
        elif last == '-3':
            if arg[0] == '=':
                exact = True
                ext_3 = int(arg[1:])
            else:
                ext_3 = int(arg)
            last = None
        elif arg in ['-5', '-3']:
            last = arg
        elif not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg
        else:
            print "Unknown option: %s" % arg
            usage()

    if not fname:
        usage('Missing input filename!')

    if not ext_5 and not ext_3:
        usage('Missing required option! -5 or -3 is required!')

    if ext_5 and ext_3 and exact:
        usage("You can't specify an exact length and 5' and 3' extensions!")

    bed_extend(BedStreamer(fname), ext_5, ext_3, exact)
