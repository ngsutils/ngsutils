#!/usr/bin/env python

import sys,gzip,os
sys.path.append(os.path.join(os.path.dirname(__file__),"..","utils"))
sys.path.append(os.path.join(os.path.dirname(__file__),"..","ext"))
from eta import ETA
import pysam

def bam_read_pos(fname):
    bamfile = pysam.Samfile(fname,"rb")
    eta = ETA(0,bamfile=bamfile)

    for read in bamfile.fetch():
        eta.print_status(extra=read.qname,bam_pos=(read.rname,read.pos))
        if not read.is_unmapped:
            print '%s\t%s\t%s\t%s\t0\t%s' % (bamfile.getrname(read.rname), read.pos, read.aend, read.qname, '-' if read.is_reverse else '+')
    
    eta.done()
    bamfile.close()

def usage():
    print """\
Usage: %s bamfile

Ouputs the read positions of all mapped reads in BED format.


""" % os.path.basename(sys.argv[0])

if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[-1]):
        usage()
        sys.exit(1)

    bam_read_pos(sys.argv[1])
