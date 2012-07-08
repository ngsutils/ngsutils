import sys
from gtf_utils import GTF


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
        if not chrom in self.ranges:
            return None
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


class RegionTagger(object):
    def __init__(self, gtf_file, valid_chroms=None):
        self.regions = []
        self.counts = {}

        sys.stderr.write('Loading gene model: %s\n' % gtf_file)
        gtf = GTF(gtf_file)
        exons = RangeMatch('exon')
        introns = RangeMatch('intron')
        promoters = RangeMatch('promoter')

        for gene in gtf.genes:
            if valid_chroms and not gene.chrom in valid_chroms:
                continue
            if gene.strand == '+':
                promoters.add_range(gene.chrom, gene.strand, gene.start - 2000, gene.start)
            else:
                promoters.add_range(gene.chrom, gene.strand, gene.end, gene.end + 2000)

            for transcript in gene.transcripts:
                last_end = None
                for start, end in transcript.exons:
                    if last_end:
                        introns.add_range(gene.chrom, gene.strand, last_end, start)
                    exons.add_range(gene.chrom, gene.strand, start, end)
                    last_end = end

        self.regions.append(exons)
        self.regions.append(introns)
        self.regions.append(promoters)

        self.counts['exon'] = 0
        self.counts['intron'] = 0
        self.counts['promoter'] = 0
        self.counts['exon-rev'] = 0
        self.counts['intron-rev'] = 0
        self.counts['promoter-rev'] = 0
        self.counts['junction'] = 0
        self.counts['intergenic'] = 0
        self.counts['mitochondrial'] = 0

    def add_read(self, read, chrom):
        if read.is_unmapped:
            return

        tag = None
        strand = '-' if read.is_reverse else '+'

        if chrom == 'chrM':
            tag = 'mitochondrial'

        if not tag:
            for op, length in read.cigar:
                if op == 3:
                    tag = 'junction'
                    break

        if not tag:
            for region in self.regions:
                tag = region.get_tag(chrom, strand, read.pos)
                if tag:
                    break

        if not tag:
            tag = 'intergenic'

        if tag:
            self.counts[tag] += 1

    def add_region(self, chrom, start, end, strand):
        tag = None

        if chrom == 'chrM':
            tag = 'mitochondrial'

        if not tag:
            for op, length in read.cigar:
                if op == 3:
                    tag = 'junction'
                    break

        if not tag:
            for region in self.regions:
                tag = region.get_tag(chrom, strand, read.pos)
                if tag:
                    break

        if not tag:
            tag = 'intergenic'

        if tag:
            self.counts[tag] += 1
