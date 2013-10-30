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
from ngsutils.support import gzip_reader


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

    first = True
    extra  = ''

    for line in gzip_reader(bedgraph, callback=lambda: extra):
        if first:
            out.write(line)  # header
            first = False
            continue
        cols = line.strip().split('\t')
        ref = cols[0]
        start = int(cols[1])
        end = int(cols[2])

        extra = '%s:%s-%s' % (ref, start, end)

        if not ref in refs:
            continue

        if start >= refs[ref]:
            # skip this... it makes no sense
            continue
        if end > refs[ref]:
            # truncate this record...
            cols[2] = refs[ref]

        out.write('%s\n' % '\t'.join([str (x) for x in cols]))

    
if __name__ == '__main__':
    fname = None
    chrom_sizes = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg
        elif not chrom_sizes and os.path.exists(arg):
            chrom_sizes = arg
        else:
            print "Unknown option: %s" % arg
            usage()

    if not fname or not chrom_sizes:
        usage()

    bedgraph_clean(fname, chrom_sizes)
