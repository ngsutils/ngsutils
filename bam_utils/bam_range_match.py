#!/usr/bin/env python
'''
For a given BAM file, what is the distribution of reads across 
exons/introns/promoter regions based upon a RefIso annotation file
'''

import sys,gzip,os
sys.path.append(os.path.join(os.path.dirname(__file__),"..","utils"))
sys.path.append(os.path.join(os.path.dirname(__file__),"..","ext"))
from eta import ETA

import pysam

bam_cigar = ['M','I','D','N','S','H','P']

class RangeMatch(object):
    def __init__(self, name):
        self.ranges = {}
        self.name = name

    def add_range(self,chrom,strand,start,end):
        if not chrom in self.ranges:
            self.ranges[chrom] = {}

        bin = start / 100000
        if not bin in self.ranges[chrom]:
            self.ranges[chrom][bin] = []
        self.ranges[chrom][bin].insert(0,(start,end,strand))
        
        if (end / 100000) != bin:
            bin = end / 100000
            if not bin in self.ranges[chrom]:
                self.ranges[chrom][bin] = []
            self.ranges[chrom][bin].insert(0,(start,end,strand))


    def get_tag(self,chrom,strand,pos, ignore_strand = False):
        bin = pos / 100000
        rev_match = False
        if not bin in self.ranges[chrom]:
            return None
        for start,end,r_strand in self.ranges[chrom][bin]:
            if pos >= start and pos <= end:
                if ignore_strand or strand == r_strand:
                    return self.name
                rev_match = True
        if rev_match:
            return "%s-rev" % self.name
        return None

def load_known_gene(exon_ranges, intron_ranges,fname, neargene_ranges=None,neargene_size=0, promoter_2k_ranges=None, promoter_5k_ranges=None):
    if fname[-3:] == '.gz':
        f = gzip.open(fname)
    else:
        f = open(fname)

    eta = ETA(os.stat(fname).st_size,fileobj = f, modulo=100)

    for line in f:
        cols = line.strip().split('\t')
        eta.print_status(extra = cols[1])
        chrom = cols[3]
        strand = cols[4]
        txStart = int(cols[5])
        txEnd = int(cols[6])
        exonStarts = cols[10].split(',')
        exonEnds = cols[11].split(',')
        if promoter_2k_ranges:
            if strand == '+':
                promoter_2k_ranges.add_range(chrom,strand,txStart-2000,txStart)
            else:
                promoter_2k_ranges.add_range(chrom,strand,txEnd,txEnd+2000)
        if promoter_5k_ranges:
            if strand == '+':
                promoter_5k_ranges.add_range(chrom,strand,txStart-5000,txStart)
            else:
                promoter_5k_ranges.add_range(chrom,strand,txEnd,txEnd+5000)
        
        if neargene_ranges and neargene_size > 0:
            neargene_ranges.add_range(chrom,strand,txStart-neargene_size,txEnd+neargene_size)
        if intron_ranges:
            intron_ranges.add_range(chrom,strand,txStart,txEnd)
        for start,end in zip(exonStarts,exonEnds):
            if start and end and exon_ranges:
                exon_ranges.add_range(chrom,strand,int(start),int(end))

    f.close()
    eta.done()

def bam_range_match(fname, ranges):
    matches = {}
    intergenic_chrom = {}
    
    samfile = pysam.Samfile(fname,"rb")
    eta = ETA(0,bamfile=samfile)

    last = None

    for read in samfile.fetch():
        if last == read.qname:
            continue
        last = read.qname
        eta.print_status(extra=read.qname,bam_pos=(read.rname,read.pos))
        chrom = samfile.getrname(read.rname)

        pos = read.pos # zero-based
        
        strand = '-' if read.is_reverse else '+'
        cigar = read.cigar
        tag = None

        if chrom not in intron_ranges.ranges:
            tag='missing'

        if chrom == 'chrM':
            tag = 'mitochondrial'

        if not tag:
            for op,length in cigar:
                if op == 3:
                    tag = 'junction'
                    break

        if not tag:
            for r in ranges:
                if not tag:
                    tag = r.get_tag(chrom,strand,pos)

        if not tag:
            tag = 'intergenic'

        if not tag in matches:
            matches[tag] = 0
        matches[tag] += 1
        
        if tag == 'intergenic':
            if not chrom in intergenic_chrom:
                intergenic_chrom[chrom] = 0
            intergenic_chrom[chrom] += 1
    eta.done()

    sorted_matches = [x for x in matches]
    sorted_matches.sort()
    for match in sorted_matches:
        print "%s\t%s" % (match,matches[match])
    
    print ""
    print "intergenic distribution"
    sorted_matches = [x for x in intergenic_chrom]
    
    for chrom in sorted_matches:
        print "%s\t%s" % (chrom,intergenic_chrom[chrom])


def usage():
    print __doc__
    print "Usage: %s bamfile refiso.txt{.gz}" % os.path.basename(sys.argv[0])

if __name__ == "__main__":
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)

    intron_ranges = RangeMatch('intron')
    exon_ranges = RangeMatch('exon')
    p2k_ranges = RangeMatch('promoter_2k')
    sys.stderr.write('Reading Known Gene definition...\n')
    load_known_gene(exon_ranges,intron_ranges,sys.argv[2], promoter_2k_ranges=p2k_ranges)
    
    sys.stderr.write('Calculating BAM matches...\n')
    bam_range_match(sys.argv[1],[exon_ranges,intron_ranges,p2k_ranges]) 
