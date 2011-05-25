#!/usr/bin/env python
"""
Given a number of BAM files, this script will merge them together, taking 
only the best matches.  There can be any number of files, but the BAM header
will be taken from the first one.  The input files should be sorted by read 
name, or at least have reads in the same order.

The first input file should have a record for every read in the other files.
However, the secondary files *may* have missing lines, so long as they are in 
the same order as the first file.
"""

import os,sys,gzip
from support.eta import ETA
import pysam

def usage():
    base = os.path.basename(sys.argv[0])
    print __doc__
    print """
Usage: %s out.bam in1.bam in2.bam ...

""" % (base,)
    sys.exit(1)

def bam_reads_batch(bam):
    reads = []
    last = None
    for read in bam:
        if last and read.qname != last:
            yield reads
            reads = []
        last = read.qname
        reads.append(read)
    
    if reads:
        yield reads

def bam_merge(fname,infiles,tag='AS'):
    bams = []
    last_reads = []
    bamgens = []
    
    for infile in infiles:
        bam=pysam.Samfile(infile,"rb")
        last_reads.append(None)
        bams.append(bam)
        bamgens.append(bam_reads_batch(bam))
        
    outfile = pysam.Samfile(fname,"wb",template=bams[0])
    i = 0
    
    while True:
        found = False
        for i,bamgen in enumerate(bamgens):
            if last_reads[i] == None:
                try:
                    last_reads[i] = bamgen.next()
                    if last_reads[i]:
                        found = True
                except:
                    pass
            else:
                found = True
        if not found:
            break
        
        best_val = None
        best_reads = None
        
        for fn,reads in zip(infiles,last_reads):
            print os.path.basename(fn),reads[0].qname,reads[0].is_unmapped,reads[0].opt(tag)
        print
        
        for i in xrange(last_reads):
            if not last_reads[i]:
                continue
            for read in last_reads[i].reads:
                if read.qname == last_reads[0][0].qname:
                    if not read.is_unmapped:
                        tag_val = int(read.opt(tag))
                        if not best_val or tag_val > best_val:
                            best_val = tag_val
                            best_reads = last_reads[i]
                            break
                    last_reads[i]=None
    
        if best_reads:
            for read in best_reads:
                outfile.write(read)
        i += 1
        if i > 100:
            return

    eta.done()
    outfile.close()
    for bam in bams:
        bam.close()


if __name__ == '__main__':
    infiles = []
    outfile = None
    num=1000000
    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif not outfile:
            outfile = arg
        elif os.path.exists(arg):
            infiles.append(arg)

    if not infiles or not outfile:
        usage()
    else:
        bam_merge(outfile,infiles)
        
