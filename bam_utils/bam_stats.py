#!/usr/bin/env python
"""
Calculate a number of summary statistics for a BAM file
"""

import os,sys,gzip
import pysam
from eta import ETA

def usage():
    base = os.path.basename(sys.argv[0])
    print __doc__
    print """
Usage: %s in.bam
""" % (base)
    sys.exit(1)

class Bins(object):
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
        

def bam_stats(infile):
    bamfile = pysam.Samfile(infile,"rb")
    eta = ETA(0,bamfile=bamfile)
    
    total = 0
    mapped = 0
    unmapped = 0
    alignments = Bins()
    edits = Bins()
    names = set()
    
    refs = {}
    for rname in bamfile.references:
        refs[rname] = 0
    
    sys.stderr.write('Calculating Read stats...\n')
    for read in bamfile:
        if read.qname in names:
            # reads only count once for this...
            continue
        
        eta.print_status(extra="%s:%s" % (bamfile.getrname(read.rname),read.pos),bam_pos=(read.rname,read.pos))
        names.add(read.qname)
        
        total += 1
        if read.is_unmapped:
            unmapped += 1
        else:
            mapped += 1
            refs[bamfile.getrname(read.rname)]+=1
            ih = int(read.opt('IH'))
            nm = int(read.opt('NM'))
        
            alignments.add(ih)
            edits.add(nm)

    eta.done()
    
    print "Reads:\t%s" % total
    print "Mapped:\t%s" % mapped
    print "Unmapped:\t%s" % unmapped
    print ""
    print "Ave # alignments (IH):\t%s" % alignments.mean()
    print "Max # alignments (IH):\t%s" % alignments.max()
    print
    print "Ave edit distance (NM):\t%s" % edits.mean()
    print "Max edit distance (NM):\t%s" % edits.max()
    print ""
    print "# of alignments"
    for i,v in enumerate(alignments.bins):
        if v:
            print "%s\t%s" % (i,v)

    print ""
    print "Edit distances"
    for i,v in enumerate(edits.bins):
        if v:
            print "%s\t%s" % (i,v)
    print ""
    
    print "Reference distribution"
    print "ref\tlength\tcount\tcount per million bases"
    for refname,reflen in zip(bamfile.references,bamfile.lengths):
        print "%s\t%s\t%s\t%s" % (refname,reflen,refs[refname],refs[refname]/(float(reflen)/1000000))

    bamfile.close()
    

if __name__ == '__main__':
    infile = None
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif not infile and os.path.exists(arg):
            infile = arg

    if not infile:
        usage()
    else:
        bam_stats(infile)
        
