#!/usr/bin/env python
#!/usr/bin/env python
## category General
## desc Find overlapping BED regions from a query and target file
'''
Finds BED regions from the query file that overlap with regions from the target file.
'''

import os
import sys
from ngsutils.bed import BedStreamer, BedFile


def usage():
    print __doc__
    print """\
Usage: bedutils overlap reference query

-nostrand           Ignore strand information when merging regions
"""
    sys.exit(1)

def bed_reduce(target_bed, query_bed, stranded=True, exact=False, out=sys.stdout):
    for qregion in query_bed:
        for tregion in target_bed.fetch(qregion.chrom, qregion.start, qregion.end, qregion.strand if stranded else None):
            if not exact or (qregion.start == tregion.start and qregion.end == tregion.end):
                qregion.write(out)
                break


if __name__ == '__main__':
    fnames = []
    stranded = True
    exact = False

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if arg == '-nostrand':
            stranded = False
        elif arg == '-exact':
            exact = True
        elif os.path.exists(arg):
            fnames.append(arg)
        else:
            print "Unknown option: %s" % arg
            usage()

    if not fnames or len(fnames) < 2:
        usage()

    bed_reduce(BedFile(fnames[0]), BedStreamer(fnames[1]), stranded, exact)
