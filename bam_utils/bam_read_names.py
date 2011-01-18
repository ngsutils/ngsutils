#!/usr/bin/env python

import sys,gzip,os
from support.eta import ETA
import pysam

def bam_read_names(fname,show_all=False):
    bamfile = pysam.Samfile(fname,"rb")
    eta = ETA(0,bamfile=bamfile)

    for read in bamfile.fetch():
        eta.print_status(extra=read.qname,bam_pos=(read.rname,read.pos))
        if show_all or not read.is_unmapped:
            sys.stdout.write(read.qname)
            sys.stdout.write('\n')
    
    eta.done()
    bamfile.close()

def usage():
    print """\
Usage: %s {-all} bamfile

Ouputs the names of all mapped reads.

Options
-all   Output all reads, even if they are unmapped

""" % os.path.basename(sys.argv[0])

if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[-1]):
        usage()
        sys.exit(1)

    if sys.argv[1] == '-all':
        bam_read_names(sys.argv[-1],True)
    else:
        bam_read_names(sys.argv[-1])
