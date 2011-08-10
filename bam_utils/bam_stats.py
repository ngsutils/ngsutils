#!/usr/bin/env python
"""
Calculate a number of summary statistics for a BAM file
"""

import os,sys,gzip
import pysam
from support.eta import ETA
from support.refiso import RefIso

def usage():
    print __doc__
    print """
Usage: bamutils stats in.bam {-model refiso.txt} {region}

If a RefIso file is given, counts corresponding to exons, introns, promoters,
junctions, intergenic, and mitochondrial regions will be calculated.

If a region is given, only reads that map to that region will be counted. 
Regions should be be in the format: 'ref:start-end' or 'ref:start' using 
1-based start coordinates.
"""
    sys.exit(1)

flag_descriptions = {
0x1: 'Multiple fragments',
0x2: 'All fragments aligned',
0x4: 'Unmapped',
0x8: 'Next unmapped',
0x10: 'Reverse complimented',
0x20: 'Next reverse complimented',
0x40: 'First fragment',
0x80: 'Last fragment',
0x100: 'Secondary alignment',
0x200: 'QC Fail',
0x400: 'PCR/Optical duplicate'
}

class RangeMatch(object):
    '''
    Simple genomic ranges.  You can define chrom:start-end ranges, then ask if a 
    particular genomic coordinate maps to any of those ranges.  This is less-
    efficient than an R-Tree, but easier to code.
    '''
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
            for bin in xrange(bin+1,(end/100000)+1):
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


class Bins(object):
    '''
    Setup simple binning.  Bins are continuous 0->max.  Values are added to 
    bins and then means / distributions can be calculated.
    '''
    def __init__(self):
        self.bins = []
    def add(self,val):
        while len(self.bins) <= val:
            self.bins.append(0)
        self.bins[val]+=1
    def mean(self):
        acc = 0
        count = 0
        for i,val in enumerate(self.bins):
            acc += (i * val)
            count += val
        
        return float(acc)/ count

    def max(self):
        return len(self.bins)-1
        

def bam_stats(infile,ref_file = None, region = None):
    bamfile = pysam.Samfile(infile,"rb")
    eta = ETA(0,bamfile=bamfile)
    
    regions = []
    counts = {}
    flag_counts = {}
    
    ref = None
    start = None
    end = None
    
    if ref_file:
        sys.stderr.write('Loading gene model: %s\n' % ref_file)
        refiso = RefIso(ref_file)
        exons = RangeMatch('exon')
        introns = RangeMatch('intron')
        promoters = RangeMatch('promoter')
        for gene in refiso.genes:
            if not gene.chrom in bamfile.references:
                continue
            if gene.strand == '+':
                promoters.add_range(gene.chrom,gene.strand,gene.tx_start-2000,gene.tx_start)
            else:
                promoters.add_range(gene.chrom,gene.strand,gene.tx_end,gene.tx_end+2000)
            
            for transcript in gene.transcripts:
                last_end = None
                for start,end in zip(transcript.exon_starts,transcript.exon_ends):
                    if last_end:
                        introns.add_range(gene.chrom,gene.strand,last_end,start)
                    exons.add_range(gene.chrom,gene.strand,start,end)
                    last_end = end
        regions.append(exons)
        regions.append(introns)
        regions.append(promoters)
        counts['exon']=0
        counts['intron']=0
        counts['promoter']=0
        counts['exon-rev']=0
        counts['intron-rev']=0
        counts['promoter-rev']=0
        counts['junction']=0
        counts['intergenic']=0
        counts['mitochondrial']=0
    
    if region:
        ref,startend=region.split(':')
        if '-' in startend:
            start,end = [int(x) for x in startend.split('-')]
            start = start - 1
            sys.stderr.write('Region: %s:%s-%s\n' % (ref,start+1,end))
        else:
            start = int(startend)-1
            end = int(startend)
            sys.stderr.write('Region: %s:%s\n' % (ref,start+1))
    
    total = 0
    mapped = 0
    unmapped = 0
    lengths = Bins()
    alignments = Bins()
    edits = Bins()
    names = set()
    refs = {}

    for rname in bamfile.references:
        refs[rname] = 0

    if region:
        def read_gen():
            for read in bamfile.fetch(ref,start,end):
                yield read
    else:
        def read_gen():
            for read in bamfile:
                yield read
    
    
    sys.stderr.write('Calculating Read stats...\n')
    try:
        for read in read_gen():
            if read.qname in names:
                # reads only count once for this...
                continue
            
            if not read.flag in flag_counts:
                flag_counts[read.flag]=1
            else:
                flag_counts[read.flag]+=1
            
            names.add(read.qname)
            lengths.add(len(read.seq))
            total += 1
            if read.is_unmapped:
                unmapped += 1
                continue

            eta.print_status(extra="%s:%s" % (bamfile.getrname(read.rname),read.pos),bam_pos=(read.rname,read.pos))
        
            mapped += 1
            refs[bamfile.getrname(read.rname)]+=1
            try:
                ih = int(read.opt('IH'))
            except:
                ih = 0
            try:
                nm = int(read.opt('NM'))
            except:
                nm = 0
    
            alignments.add(ih)
            edits.add(nm)
        
            if regions:
                tag = None
                chrom = bamfile.getrname(read.rname)
                strand = '-' if read.is_reverse else '+'

                if chrom == 'chrM':
                    tag = 'mitochondrial'
            
                if not tag:
                    for op,length in read.cigar:
                        if op == 3:
                            tag = 'junction'
                            break
                        
                if not tag:
                    for region in regions:
                        tag = region.get_tag(chrom,strand,read.pos)
                        if tag:
                            break
        
                if not tag:
                    tag = 'intergenic'

                if tag:
                    counts[tag] += 1
    except KeyboardInterrupt:
        pass
        
    eta.done()

    print "Reads:\t%s" % total
    print "Mapped:\t%s" % mapped
    print "Unmapped:\t%s" % unmapped
    
    if total > 0:
        print ""
        print "Flag distribution"
        
        tmp = []
        maxsize = 0
        for flag in flag_descriptions:
            if flag in flag_counts:
                tmp.append(flag)
                maxsize = max(maxsize,len(flag_descriptions[flag]))
        tmp.sort()
        
        for flag in tmp:
            count = 0
            for f in flag_counts:
                if (f & flag) > 0:
                    count += flag_counts[f]
                    
            if count > 0:
                print "[0x%03x] %-*s:\t%s (%.1f%%)" % (flag,maxsize,flag_descriptions[flag],count,(float(count)*100/total))
                
        print ""
        print ""
        print "Ave length:\t%s" % lengths.mean()
        print ""
        print "Ave # alignments (IH):\t%s" % alignments.mean()
        print "Max # alignments (IH):\t%s" % alignments.max()
        print
        print "Ave edit distance (NM):\t%s" % edits.mean()
        print "Max edit distance (NM):\t%s" % edits.max()
        print ""
        print "Read lengths"
        for i,v in enumerate(lengths.bins[::-1]):
            if v:
                print "%s\t%s" % (lengths.max()-i,v)
        print ""
        print "# of alignments (IH)"
        for i,v in enumerate(alignments.bins):
            if v:
                print "%s\t%s" % (i,v)

        print ""
        print "Edit distances (NM)"
        for i,v in enumerate(edits.bins):
            if v:
                print "%s\t%s" % (i,v)
        print ""
    
        print "Reference distribution"
        print "ref\tlength\tcount\tcount per million bases"
        for refname,reflen in zip(bamfile.references,bamfile.lengths):
            print "%s\t%s\t%s\t%s" % (refname,reflen,refs[refname],refs[refname]/(float(reflen)/1000000))

        if regions:
            print ""
            print "Mapping regions"
            sorted_matches = [x for x in counts]
            sorted_matches.sort()
            for match in sorted_matches:
                print "%s\t%s" % (match,counts[match])
        
        
    
    bamfile.close()
    

if __name__ == '__main__':
    infile = None
    refiso = None
    region = None
    
    last = None
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif not infile and os.path.exists(arg):
            infile = arg
        elif last == '-model':
            refiso = arg
            last = None
        elif arg in ['-model']:
            last = arg
        else:
            region = arg

    if not infile:
        usage()
    else:
        bam_stats(infile,refiso,region)
        
