#!/usr/bin/env python
"""
Calculate a number of summary statistics / distribution for a BED file
"""

import os
import sys

from support.eta import ETA
from support.refiso import RefIso


def usage():
    print __doc__
    print """
Usage: bedutils stats in.bed {refiso.txt}

If a RefIso file is given, counts corresponding to exons, introns, promoters,
junctions, intergenic, and mitochondrial regions will be calculated.
"""
    sys.exit(1)


def format_number(n):
    ar = list(str(n))
    for i in range(len(ar))[::-3][1:]:
        ar.insert(i + 1, ',')
    return ''.join(ar)


class RangeMatch(object):
    '''
    Simple genomic ranges.  You can define chrom:start-end ranges, then ask if a
    particular genomic coordinate maps to any of those ranges.  This is less-
    efficient than an R-Tree, but easier to code.
    '''
    def __init__(self, name):
        self.ranges = {}
        self.name = name

    def add_range(self, chrom, strand, start, end):
        if not chrom in self.ranges:
            self.ranges[chrom] = {}

        bin = start / 100000
        if not bin in self.ranges[chrom]:
            self.ranges[chrom][bin] = []
        self.ranges[chrom][bin].insert(0, (start, end, strand))

        if (end / 100000) != bin:
            for bin in xrange(bin + 1, (end / 100000) + 1):
                if not bin in self.ranges[chrom]:
                    self.ranges[chrom][bin] = []
                self.ranges[chrom][bin].insert(0, (start, end, strand))

    def get_tag(self, chrom, strand, pos, ignore_strand=False):
        bin = pos / 100000
        rev_match = False
        if not bin in self.ranges[chrom]:
            return None
        for start, end, r_strand in self.ranges[chrom][bin]:
            if pos >= start and pos <= end:
                if ignore_strand or strand == r_strand:
                    return self.name
                rev_match = True
        if rev_match:
            return "%s-rev" % self.name
        return None


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


def bed_stats(infile, ref_file=None):
    regions = []
    counts = {}
    if ref_file:
        sys.stderr.write('Loading gene model: %s\n' % ref_file)
        refiso = RefIso(ref_file)
        exons = RangeMatch('exon')
        introns = RangeMatch('intron')
        promoters = RangeMatch('promoter')
        utr5 = RangeMatch('utr-5')
        utr3 = RangeMatch('utr-3')
        for gene in refiso.genes:
            for transcript in gene.transcripts:
                if gene.strand == '+':
                    promoters.add_range(gene.chrom, gene.strand, transcript.tx_start - 2000, transcript.tx_start)
                    utr5.add_range(gene.chrom, gene.strand, transcript.tx_start, transcript.cds_start)
                    utr3.add_range(gene.chrom, gene.strand, transcript.tx_end, transcript.cds_end)
                else:
                    promoters.add_range(gene.chrom, gene.strand, transcript.tx_end, transcript.tx_end + 2000)
                    utr5.add_range(gene.chrom, gene.strand, transcript.tx_end, transcript.cds_end)
                    utr3.add_range(gene.chrom, gene.strand, transcript.tx_start, transcript.cds_start)
                last_end = None
                for start, end in zip(transcript.exon_starts, transcript.exon_ends):
                    if last_end:
                        introns.add_range(gene.chrom, gene.strand, last_end, start)

                    exons.add_range(gene.chrom, gene.strand, start, end)
                    last_end = end

        regions.append(exons)
        regions.append(introns)
        regions.append(promoters)
        regions.append(utr5)
        regions.append(utr3)

        counts['exon'] = 0
        counts['intron'] = 0
        counts['promoter'] = 0
        counts['exon-rev'] = 0
        counts['intron-rev'] = 0
        counts['promoter-rev'] = 0
        counts['intergenic'] = 0
        counts['mitochondrial'] = 0
        counts['utr-3'] = 0
        counts['utr-5'] = 0

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

                if regions:
                    tag = None

                    if chrom == 'chrM':
                        tag = 'mitochondrial'

                    if not tag:
                        for region in regions:
                            tag = region.get_tag(chrom, strand, start)
                            if tag:
                                break

                    if not tag:
                        tag = 'intergenic'

                    if tag:
                        counts[tag] += 1
        except KeyboardInterrupt:
            pass

        eta.done()

    print "Regions:\t%s" % format_number(total)
    print "Total coverage:\t%s bases" % format_number(size)
    print "Average size:\t%s bases" % lengths.mean()
    print ""
    print "Reference distribution"
    print "ref\tcount"
    refnames = []
    for refname in refs:
        refnames.append(refname)
    refnames.sort()
    for refname in refnames:
        print "%s\t%s" % (refname, format_number(refs[refname]))

    if regions:
        print ""
        print "Mapping regions"
        sorted_matches = [x for x in counts]
        sorted_matches.sort()
        for match in sorted_matches:
            print "%s\t%s" % (match, format_number(counts[match]))


if __name__ == '__main__':
    infile = None
    refiso = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif not infile and os.path.exists(arg):
            infile = arg
        elif not refiso and os.path.exists(arg):
            refiso = arg

    if not infile:
        usage()
    else:
        bed_stats(infile, refiso)
