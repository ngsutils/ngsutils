#!/usr/bin/env python

import sys,gzip,os
sys.path.append(os.path.join(os.path.dirname(__file__),"..","utils"))
sys.path.append(os.path.join(os.path.dirname(__file__),"..","ext"))
from eta import ETA
import pysam

bam_cigar = ['M','I','D','N','S','H','P']

def bam_read_pos(fname,all=False):
    bamfile = pysam.Samfile(fname,"rb")
    eta = ETA(0,bamfile=bamfile)

    for read in bamfile.fetch():
        eta.print_status(extra=read.qname,bam_pos=(read.rname,read.pos))
        if not read.is_unmapped:
            print '%s\t%s\t%s\t%s' % (read.qname, bamfile.getrname(read.rname), read.pos, ''.join(['%s%s' % (length,bam_cigar[op]) for op,length in read.cigar]))
        elif all:
            print '%s\t*\t0\t\n' % (read.qname)
    
    eta.done()
    bamfile.close()

def usage():
    print """\
Usage: %s bamfile

Ouputs the read positions of all mapped reads in a tab-delimited format.

Format:
name    ref    pos    cigar


""" % os.path.basename(sys.argv[0])

if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[-1]):
        usage()
        sys.exit(1)

    if sys.argv[1] == '-all':
        bam_read_pos(sys.argv[2],True)
    else:
        bam_read_pos(sys.argv[1])

