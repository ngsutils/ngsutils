#!/usr/bin/env python
## category General
## desc Calculates simple stats for a BED file
"""
Calculate a number of summary statistics / distribution for a BED file
"""

import os
import sys

from ngsutils.support.regions import RegionTagger
from ngsutils.support.ngs_utils import format_number
from ngsutils.support import Counts
from ngsutils.bed import BedFile
import ngsutils.support.ngs_utils


def usage():
    print __doc__
    print """
Usage: bedutils stats  {-names} in.bed {gene_annotation.gtf}

If a GTF file is given, counts corresponding to exons, introns, promoters,
junctions, intergenic, and mitochondrial regions will be calculated.
"""
    sys.exit(1)


class BedStats(object):
    def __init__(self, bed, gtf_file=None, names=False):
        self.regiontagger = None
        if gtf_file:
            self.regiontagger = RegionTagger(gtf_file)

        self.total = 0
        self.size = 0
        self.lengths = Counts()
        self.refs = {}
        self.names = {}
        for region in bed:
            self.total += 1
            self.size += (region.end - region.start)
            self.lengths.add(region.end - region.start)

            if names:
                if not region.name in self.names:
                    self.names[region.name] = 1
                else:
                    self.names[region.name] += 1

            if not region.chrom in self.refs:
                self.refs[region.chrom] = 0
            self.refs[region.chrom] += 1

            if self.regiontagger:
                self.regiontagger.add_region(region.chrom, region.start, region.end, region.strand)

    def write(self, out=sys.stdout):
        out.write("Regions:\t%s\n" % format_number(self.total))
        out.write("Total coverage:\t%s bases\n" % format_number(self.size))
        out.write("Average size:\t%s bases\n" % self.lengths.mean())
        out.write("\n")
        out.write("Reference distribution\n")
        out.write("ref\tcount\n")
        for refname in ngsutils.support.ngs_utils.natural_sort([x for x in self.refs]):
            out.write("%s\t%s\n" % (refname, format_number(self.refs[refname])))

        if self.names:
            out.write("\nName distribution\n")
            out.write("name\tcount\n")
            for name in ngsutils.support.ngs_utils.natural_sort([x for x in self.names]):
                out.write("%s\t%s\n" % (name, format_number(self.names[name])))

        if self.regiontagger:
            out.write("\n")
            out.write("Mapping regions\n")
            sorted_keys = [x for x in self.regiontagger.counts]
            sorted_keys.sort()
            for k in sorted_keys:
                out.write("%s\t%s\n" % (k, self.regiontagger.counts[k]))


def bed_stats(infile, gtf_file=None, out=sys.stdout, quiet=False, names=False):
    if not quiet:
        sys.stderr.write('Calculating BED region stats...\n')

    stats = BedStats(BedFile(infile), gtf_file, names=names)
    stats.write(out)

if __name__ == '__main__':
    infile = None
    gtf_file = None
    names = False

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif arg == '-names':
            names = True
        elif not infile and (os.path.exists(arg) or arg == '-'):
            infile = arg
        elif not gtf_file and os.path.exists(arg):
            gtf_file = arg

    if not infile:
        usage()
    else:
        bed_stats(infile, gtf_file, names=names)
