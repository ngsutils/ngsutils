#!/usr/bin/env python
"""
Given a BAM file and a BED file, this script will extract only reads that
map to the given BED regions.

"""

import os,sys,gzip
from support.eta import ETA
import pysam

def usage():
    base = os.path.basename(sys.argv[0])
    print __doc__
    print """
Usage: %s {opts} in.bam out.bam regions.bed

Options:
  -ns    Ignore strandedness of reads and regions
""" % (base,base)
    sys.exit(1)


def bam_extract(infile,outfile,bedfile,nostrand=False):
    with open(bedfile) as f:
        bamfile = pysam.Samfile(infile,"rb")
        outfile = pysam.Samfile(outfile,"wb",template=bamfile)
        if failedfile:
            failed_out = open(failedfile,'w')
        else:
            failed_out = None
        eta = ETA(os.stat(bedfile).st_size,fileobj=f)
    
        passed = 0
    
        for line in bedfile:
            eta.print_status(extra="extracted:%s" % (passed,))
            if line[0] == '#':
                continue
            cols = line.strip().split('\t')
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            strand = cols[5]
            
            for read in bamfile.fetch(chrom,start,end):
                if nostrand or (read.is_reverse and strand == '-') or (not read.is_reverse and strand == '+'):
                    if start <= read.pos <= end or start <= read.apos <= end:
                        outfile.write(read)
                        passed += 1

        eta.done()
        bamfile.close()
        outfile.close()

    sys.stderr.write("%s extracted\n" % (passed,))

if __name__ == '__main__':
    infile = None
    outfile = None
    bedfile = None
    nostrand = False
    
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif arg == '-ns':
            nostrand = True
        elif not infile and os.path.exists(arg):
            infile = arg
        elif not outfile:
            outfile = arg
        elif not bedfile and os.path.exists(arg):
            bedfile = arg
        else:
            print "Unknown argument: %s" % arg
            usage()

    if not infile or not outfile or not bedfile:
        usage()
    else:
        bam_extract(infile,outfile,bedfile,nostrand)
        
