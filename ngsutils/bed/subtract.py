#!/usr/bin/env python
## category General
## desc Subtracts one set of BED regions from another
'''
Removes overlaping BED regions from another BED file.
'''

import os
import sys
from ngsutils.bed import BedFile, BedStreamer

def usage(msg=None):
    if msg:
        print 'ERROR: %s' % msg

    print __doc__
    print """\
Usage: bedutils subtract {opts} bedfile1 bedfile2

Returns: 
    bedfile1 - bedfile2

Options:
    -nostrand       Ignore strand information when merging regions
"""
    sys.exit(1)



def bed_subtract(mainbed, removebed, stranded=True, out=sys.stdout):
    for mainregion in mainbed:
            regionlist = [mainregion,]
            changed = False

            while regionlist:
                cur = regionlist[0]
                regionlist = regionlist[1:]
                skip = False

                newstart = cur.start
                newend = cur.end

                for region in removebed.fetch(cur.chrom, cur.start, cur.end, cur.strand if stranded else None):
                    if region.start < newstart < newend < region.end:
                        # subtract out the entire region
                        skip = True
                        break
                    
                    if newstart < region.start < region.end < newend:
                        # subtract a subregion - add the end to the list.
                        regionlist.append(cur.clone(start=region.end, end=newend, name='%s*' % cur.name))
                        newend = region.start
                        changed = True
                    elif region.start < newstart < region.end:
                        newstart = region.end
                        changed = True
                    elif region.start < newend < region.end:
                        newend = region.start
                        changed = True

                if not skip:
                    cur.start = newstart
                    cur.end = newend
                    if changed:
                        cur.name = "%s*" % cur.name

                    out.write('%s\n' % cur)


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

    if not fname1 or not fname2:
        usage()

    if fname1 == '-' and fname2 == '-':
        usage("Both input files can't be from stdin!")

    bed1 = BedStreamer(fname1)
    bed2 = BedFile(fname2)

    bed_subtract(bed1, bed2, stranded)

    bed1.close()
    bed2.close()
