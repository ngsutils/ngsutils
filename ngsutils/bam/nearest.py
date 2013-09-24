#!/usr/bin/env python
## category General
## desc Find the nearest annotated BED region for each read
## experimental
'''
For each read, find the nearest annotated region from a BED file.

For example, if you have a BED file with the transcription starting site (TSS)
for all genes in the genome, this will find the nearest TSS for each read, and
report the distance to the site.

'''

import sys
import os
import ngsutils.bam
from ngsutils.bed import BedFile


def find_nearest(bam, bed, maxdist=100000, out=sys.stdout):
    for read in ngsutils.bam.bam_iter(bam):
        if read.is_unmapped:
            continue

        dists = []  # will be an list or tuples (abs_val, signed_val)
                    # so that we can just sort the list by abs val, then
                    # get the signed value

        chrom = bam.getrname(read.tid)
        strand = '-' if read.is_reverse else '+'
        start = min(0, read.pos - maxdist)
        end = max(read.apos + maxdist, bam.lengths[read.tid])

        for region in bed.fetch(chrom, start, end, strand):
            if region.start <= read.pos <= region.end:
                # start is w/in region
                dists.append((0, 0))
                continue
            elif region.start <= read.apos <= region.end:
                # end is w/in region
                dists.append((0, 0))
                continue
            elif read.pos <= region.start <= region.end <= read.apos:
                # read spans the entire region
                dists.append((0, 0))
                continue
            elif region.end < read.pos:
                dists.append((read.pos - region.end), (read.pos - region.end))
            else:
                dists.append((region.start - read.apos, read.apos - region.start)) # should be negative

        if dists:
            dists.sort()
            out.write('%s\t%s\n' % (read.qname, dists[0][0]))
        else:
            out.write('%s\t*\n' % (read.qname,))


def usage(msg=None):
    if msg:
        sys.stdout.write('%s\n\n' % msg)
    sys.stdout.write(__doc__)
    sys.stdout.write('''\
Usage: bamutils nearest {-max val} filename.bam regions.bed

Options:
  -max    The maximal distance to look for a nearest region
          (default: 100K)

''')

    sys.exit(1)

if __name__ == '__main__':
    bam_fname = None
    bed_fname = None
    maxdist = 100000

    last = None

    for arg in sys.argv[1:]:
        if last == '-max':
            maxdist = int(arg)
            last = None
        elif arg in ['-max']:
            last = arg
        elif not bam_fname and os.path.exists(arg):
            bam_fname = arg
        elif not bed_fname and os.path.exists(arg):
            bed_fname = arg

    if not bam_fname or not bed_fname:
        usage()

    bed = BedFile(bed_fname)
    bam = ngsutils.bam.bam_open(bam_fname)

    find_nearest(bam, bed, maxdist)
