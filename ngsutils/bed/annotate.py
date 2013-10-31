#!/usr/bin/env python
## category Conversion
## desc Annotate BED files by adding / altering columns
'''
Convert a BED3 file to a BED6 file with constant name, score, and strand.
This can also be used to replace values in a BED6 file, or add RGB colors based 
on the "name" field of a region.

If you add an RGB column, then the thickStart and thickEnd columns will also be
added and set to the region start and end respectfully (if they are missing).
'''

import os
import sys
from ngsutils.bed import BedFile


def usage():
    print __doc__
    print """\
Usage: bedutils annotate {-name val} {-score val} {-strand val} {-rgb name value} bedfile

"""
    sys.exit(1)


def bed_annotate(bed, name=None, score=None, strand=None, rgb=None, out=sys.stdout):
    for region in bed:
        if name:
            region.name = name
        if not region.name:
            region.name = '%s:%s-%s' % (region.chrom, region.start, region.end)
        if score:
            region.score = score
        if strand:
            region.strand = strand
        if rgb:
            region.rgb = rgb[region.name]

        region.write(out)

if __name__ == '__main__':
    fname = None
    name = None
    score = None
    strand = None
    rgb = {}
    rgb_name = None

    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if last == '-name':
            name = arg
            last = None
        elif last == '-score':
            score = arg
            last = None
        elif last == '-strand':
            strand = arg
            last = None
        elif last == '-rgb':
            if not rgb_name:
                rgb_name = arg
            else:
                rgb[rgb_name] = arg
                rgb_name = None
                last = None
        elif arg in ['-name', '-score', '-strand', '-rgb']:
            last = arg
        elif not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg
        else:
            print "Unknown option: %s" % arg
            usage()

    if not fname:
        usage()

    bed_annotate(BedFile(fname), name=name, score=score, strand=strand, rgb=rgb)
