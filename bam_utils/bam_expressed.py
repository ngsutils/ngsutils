#!/usr/bin/env python
'''
Scans a BAM file, looking for regions of mapped reads. Reads in a BAM file, 
and outputs regions of expression in BED format.
'''
import sys,gzip,os
from support.eta import ETA

import pysam

bam_cigar = ['M','I','D','N','S','H','P']


class ExpressedRegion(object):
    _count = 0
    def __init__(self,chrom,only_uniq_starts = False):
        ExpressedRegion._count += 1
        self.name = 'region_%s' % ExpressedRegion._count
        
        self.chrom = chrom
        self.start = None
        self.end = None
        self.fwd_count = 0
        self.rev_count = 0
        
        self.reads = set()
        self.read_count = 0
        
        self.only_uniq_starts = only_uniq_starts
        self.uniq_starts = set()
        
    def add_column(self,read,pos):
        if not self.start:
            self.start = pos
        if not self.end or pos >= self.end:
            self.end = pos+1

        if not read.alignment.qname in self.reads:
            self.reads.add(read.alignment.qname)
            self.read_count += 1

            if self.only_uniq_starts:
                if not read.alignment.is_reverse:
                    if not read.alignment.pos in self.uniq_starts:
                        self.uniq_starts.add(read.alignment.pos)
                        self.fwd_count += 1
                else:
                    if not read.alignment.aend in self.uniq_starts:
                        self.uniq_starts.add(read.alignment.aend)
                        self.rev_count += 1
                    
            else:
                if read.alignment.is_reverse:
                    self.rev_count += 1
                else:
                    self.fwd_count += 1
                    
    
    def write(self,fs):
        cols = []
        cols.append(self.chrom)
        cols.append(str(self.start))
        cols.append(str(self.end))
        cols.append(self.name)
        if self.only_uniq_starts:
            cols.append(str(len(self.uniq_starts)))
        else:
            cols.append(str(self.read_count))

        if self.fwd_count > self.rev_count:
            cols.append('+')
        else:
            cols.append('-')
        
        fs.write('\t'.join(cols))
        fs.write('\n')

def bam_find_regions(bam_name, merge_distance = 10, min_read_count = 2, only_uniq_starts = False, nostrand = False):

    bamfile = pysam.Samfile(bam_name,"rb")
    eta = ETA(0,bamfile=bamfile)

    region_plus = None
    region_minus = None

    for pileup in bamfile.pileup(stepper='all', mask=1540):
        chrom = bamfile.getrname(pileup.tid)
        eta.print_status(extra='%s:%s' % (chrom,pileup.pos),bam_pos=(pileup.tid,pileup.pos))
        
        for read in pileup.pileups:
            if nostrand or not read.alignment.is_reverse:
                if not region_plus or region_plus.chrom != chrom or (region_plus.end+merge_distance) < pileup.pos:
                    if region_plus and region_plus.read_count >= min_read_count:
                        region_plus.write(sys.stdout)
                
                    region_plus = ExpressedRegion(chrom,only_uniq_starts)
                region_plus.add_column(read,pileup.pos)
            else:
                if not region_minus or region_minus.chrom != chrom or (region_minus.end+merge_distance) < pileup.pos:
                    if region_minus and region_minus.read_count >= min_read_count:
                        region_minus.write(sys.stdout)
                
                    region_minus = ExpressedRegion(chrom,only_uniq_starts)
                region_minus.add_column(read,pileup.pos)
        
    if region_plus and region_plus.read_count >= min_read_count:
        region_plus.write(sys.stdout)
    if region_minus and region_minus.read_count >= min_read_count:
        region_minus.write(sys.stdout)
    
    eta.done()
    bamfile.close()


def usage():
    print __doc__
    print """\
Usage: %s {options} bamfile

Options:
-ns             Ignore strandedness when creating regions
                (default: false)

-uniq           Only use unique starting positions when performing counts
                (default: false)

-dist N         minimum distance required between regions - they will be 
                merged if w/in this distance
                (default: 10)
           
-mincount N     The minimum number of reads required in a region
                (default: 2)

""" % os.path.basename(sys.argv[0])
    sys.exit(1)

if __name__ == "__main__":
    uniq = False
    dist = 10
    mincount = 2
    fname = None
    last = None
    nostrand = False
    
    for arg in sys.argv[1:]:
        if last == '-dist':
            dist = int(arg)
            last = None
        elif last == '-mincount':
            mincount = int(arg)
            last = None
        elif arg == '-h':
            usage()
        elif arg == '-uniq':
            uniq = True
        elif arg == '-ns':
            nostrand = True
        elif arg in ['-dist','-mincount']:
            last = arg
        elif not fname and os.path.exists(arg):
            fname = arg
        else:
            print 'Unknown argument: %s' % arg
            usage()
    if not fname:
        usage()
        
    bam_find_regions(fname,dist,mincount,uniq,nostrand)
