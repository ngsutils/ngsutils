#!/usr/bin/env python
## category General
## desc Find the nearest annotation for a BED file
## experimental
'''
For each BED region, find the nearest annotated region from a reference file.
'''

import sys
import os
import ngsutils.bam
from ngsutils.bed import BedFile


def find_nearest(inbed, refbed, maxdist=100000, restrict_name=False, nostrand=False, out=sys.stdout):
    for qregion in inbed:

        dists = []  # will be an list tuples: (abs_val of the distance, 'up/down') (respective to + strand)

        start = max(0, qregion.start - maxdist)
        end = qregion.end + maxdist

        for region in refbed.fetch(qregion.chrom, start, end, qregion.strand if not nostrand else None):
            if restrict_name:
                if restrict_name and qregion.name not in region.name:
                    continue

            if region.start <= qregion.start <= region.end:
                # start is w/in region
                dists.append((0, '', region))
                continue
            elif region.start <= qregion.end <= region.end:
                # end is w/in region
                dists.append((0, '', region))
                continue
            elif qregion.start <= region.start <= region.end <= qregion.end:
                # read spans the entire region
                dists.append((0,'', region))
                continue
            elif region.end < qregion.start:
                dists.append((qregion.start - region.end, 'up', region))
            else:
                dists.append((region.start - qregion.end, 'down', region))
            

        if dists:
            dists.sort()
            best_dist = None
            best_name = None
            for dist in dists:
                distance = dist[0]
                orient = dist[1]
                region = dist[2]

                if distance > 0:
                    if region.strand == '+' and orient == 'down':
                        distance = -distance
                    elif region.strand == '-' and orient == 'up':
                        distance = -distance
                    
                if (best_dist is None or distance == best_dist) and region.name != best_name:
                    out.write('%s\t%s\t%s\t%s\t%s\n' % (qregion.chrom, qregion.start, qregion.end, region.name, distance))
                    best_dist = distance
                    best_name = region.name
        else:
            out.write('%s\t%s\t%s\t*\t\n' % (qregion.chrom, qregion.start, qregion.end))


def usage(msg=None):
    if msg:
        sys.stdout.write('%s\n\n' % msg)
    sys.stdout.write(__doc__)
    sys.stdout.write('''\
Usage: bedutils nearest {-max val} query.bed reference.bed

Options:
  -max    The maximal distance to look for a nearest region
          (default: 100K)

  -match  Only use regions in the reference that contain the name
          from the query file.

''')

    sys.exit(1)

if __name__ == '__main__':
    qbed_fname = None
    refbed_fname = None
    maxdist = 100000
    match = False
    nostrand = False

    last = None

    for arg in sys.argv[1:]:
        if last == '-max':
            maxdist = int(arg)
            last = None
        elif arg in ['-max']:
            last = arg
        elif arg == '-match':
            match = True
        elif arg == '-nostrand':
            nostrand = True
        elif not qbed_fname and (os.path.exists(arg) or arg == '-'):
            qbed_fname = arg
        elif not refbed_fname and os.path.exists(arg):
            refbed_fname = arg

    if not qbed_fname or not refbed_fname:
        usage()

    qbed = BedFile(qbed_fname)
    refbed = BedFile(refbed_fname)
    find_nearest(qbed, refbed, maxdist, match, nostrand)
