#!/usr/bin/env python
## category General
## desc Subtracts one set of BED regions from another
'''
Removes overlaping BED regions from another BED file.
'''

import os
import sys
from ngsutils.bed import BedFile
from ngsutils.support.ngs_utils import gzip_opener

def usage(msg=None):
    if msg:
        print msg
        print ""

    print __doc__
    print """\
Usage: bedutils subtract {opts} bedfile1 bedfile2

Returns: 
    bedfile1 - bedfile2

Options:
    -nostrand       Ignore strand information when merging regions
"""
    sys.exit(1)



def bed_subtract(fname, removebed, stranded=True, out=sys.stdout):
    with gzip_opener(fname) as fobj:
        for line in fobj:
            cols = line.strip().split('\t')
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            strand = cols[5] if stranded else None

            newstart = start
            newend = end
            skip = False

            for region in removebed.fetch(chrom, start, end, strand):
                if region.start < newstart < newend < region.end:
                    skip = True
                    break

                if region.start < newstart < region.end:
                    newstart = region.end
                elif region.start < newend < region.end:
                    newend = region.start


            if not skip:
                cols[1] = str(newstart)
                cols[2] = str(newend)
                out.write('%s\n' % '\t'.join(cols))


if __name__ == '__main__':
    fname1 = None
    fname2 = None
    stranded = True

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif arg == '-nostrand':
            stranded = False
        elif not fname1 and (arg == '-' or os.path.exists(arg)):
            fname1 = arg
        elif not fname2 and (arg == '-' or os.path.exists(arg)):
            fname2 = arg
        else:
            print "Unknown option: %s" % arg
            usage()

    if not fname1 or fname2:
        usage()

    if fname1 == '-' and fname2 == '-':
        usage("Both input files can't be from stdin!")

    bed2 = BedFile(fname2)
    bed_subtract(fname1, bed2, stranded)
    bed2.close()
