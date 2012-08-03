#!/usr/bin/env python
## category General
## desc Calculates simple stats for a BED file
"""
Calculate a number of summary statistics / distribution for a BED file
"""

import os
import sys

from ngsutils.support.eta import ETA
from ngsutils.support.regions import RegionTagger
import ngsutils.support.ngs_utils


def usage():
    print __doc__
    print """
Usage: bedutils stats in.bed {gene_annotation.gtf}

If a GTF file is given, counts corresponding to exons, introns, promoters,
junctions, intergenic, and mitochondrial regions will be calculated.
"""
    sys.exit(1)


def format_number(n):
    ar = list(str(n))
    for i in range(len(ar))[::-3][1:]:
        ar.insert(i + 1, ',')
    return ''.join(ar)


class Bins(object):
    '''
    Setup simple binning.  Bins are continuous 0->max.  Values are added to
    bins and then means / distributions can be calculated.
    '''
    def __init__(self):
        self.bins = []

    def add(self, val):
        while len(self.bins) <= val:
            self.bins.append(0)
        self.bins[val] += 1

    def mean(self):
        acc = 0
        count = 0

        for i, val in enumerate(self.bins):
            acc += (i * val)
            count += val

        return float(acc) / count

    def max(self):
        return len(self.bins) - 1


def bed_stats(infile, gtf_file=None):
    regiontagger = None
    if gtf_file:
        regiontagger = RegionTagger(gtf_file)

    total = 0
    size = 0
    lengths = Bins()
    refs = {}

    sys.stderr.write('Calculating BED region stats...\n')

    with open(infile) as f:
        eta = ETA(os.stat(infile).st_size, fileobj=f)

        try:
            for line in f:
                if line[0] == '#':
                    continue
                cols = line.strip().split('\t')
                chrom = cols[0]
                start = int(cols[1])
                end = int(cols[2])
                #name = cols[3]
                #score = cols[4]
                strand = cols[5]
                eta.print_status(extra="%s:%s-%s" % (chrom, start, end))
                lengths.add(end - start)
                total += 1

                size += (end - start)

                if not chrom in refs:
                    refs[chrom] = 0
                refs[chrom] += 1

                if regiontagger:
                    regiontagger.add_region(chrom, start, end, strand)

        except KeyboardInterrupt:
            pass

        eta.done()

    print "Regions:\t%s" % format_number(total)
    print "Total coverage:\t%s bases" % format_number(size)
    print "Average size:\t%s bases" % lengths.mean()
    print ""
    print "Reference distribution"
    print "ref\tcount"
    for refname in ngsutils.support.ngs_utils.natural_sort([x for x in refs]):
        print "%s\t%s" % (refname, format_number(refs[refname]))

    if regiontagger:
        print ""
        print "Mapping regions"
        sorted_keys = [x for x in regiontagger.counts]
        sorted_keys.sort()
        for k in sorted_keys:
            print "%s\t%s" % (k, regiontagger.counts[k])


if __name__ == '__main__':
    infile = None
    gtf_file = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif not infile and os.path.exists(arg):
            infile = arg
        elif not gtf_file and os.path.exists(arg):
            gtf_file = arg

    if not infile:
        usage()
    else:
        bed_stats(infile, gtf_file)
