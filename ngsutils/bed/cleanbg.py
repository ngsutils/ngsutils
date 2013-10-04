#!/usr/bin/env python
## category Misc
## desc Cleans up a bedgraph file
'''
Cleanup a bedgraph file, ensuring that the coordinates are w/in the same
range as the associated chrom.sizes file.
'''

import os
import sys
from ngsutils.bed import BedFile, BedRegion


def usage():
    print __doc__
    print """\
Usage: bedutils cleanbg bedfile chrom.sizes

Note: a samtools faidx (.fai) file also works as a chrom.sizes file

"""
    sys.exit(1)


def bedgraph_clean(bedgraph, chrom_sizes, out=sys.stdout):
    refs = {}
    with open(chrom_sizes) as f:
        for line in f:
            cols = line.strip().split('\t')
            refs[cols[0]] = int(cols[1])

    with open(bedgraph) as f:
        for line in f:
            cols = line.strip().split('\t')
            ref = cols[0]
            start = int(cols[1])
            end = int(cols[2])

            if start >= refs[ref]:
                # skip this... it makes no sense
                continue
            if end > refs[ref]:
                # truncate this record...
                end = refs[ref]

            out.write('%s\n' % '\t'.join([str (x) for x in cols]))

    
if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if not fname and os.path.exists(arg):
            fname = arg
        else:
            print "Unknown option: %s" % arg
            usage()

    if not fname:
        usage()

    bed_clean(BedFile(fname))
