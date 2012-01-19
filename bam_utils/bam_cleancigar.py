#!/usr/bin/env python
'''
Fixes BAM files where the CIGAR alignment has a zero length element
'''
import sys
import os
import pysam

def bam_cleancigar(infile, outfile):
    bam = pysam.Samfile(infile,"rb")
    out = pysam.Samfile(outfile,"wb", template=bam)
    for read in bam:
        newcigar = []
        changed = False
        for op, length in read.cigar:
            if length > 0:
                newcigar.append((op,length))
            else:
                changed = True
        
        if changed:
            read.cigar = newcigar
        
        out.write(read)
        
    bam.close()
    out.close()
    
    

def usage():
    print __doc__
    print "Usage: bamutils cleancigar inbamfile outbamfile"
    sys.exit(-1)

if __name__ == "__main__":
    infile = None
    outfile = None
    
    for arg in sys.argv[1:]:
        if arg == "-h":
            usage()
        elif os.path.exists(arg):
            if not infile:
                infile = arg
            elif not outfile:
                outfile = arg
            else:
                usage()
        else:
            sys.stderr.write("File: %s not found!\n" % arg)
            usage()
            
    if not infile or not outfile:
        usage()

    bam_cleancigar(infile,outfile)
