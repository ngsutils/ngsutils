#!/usr/bin/env python

import sys,gzip,os
sys.path.append(os.path.join(os.path.dirname(__file__),"..","utils"))
sys.path.append(os.path.join(os.path.dirname(__file__),"..","ext"))
from eta import ETA
import pysam

def bam2fasta(fname,show_all=False):
    if sam_fname[-4:].lower == '.bam':
        samfile = pysam.Samfile(sam_fname,'rb')
    else:
        samfile = pysam.Samfile(sam_fname,'r')

    eta = ETA(0,bamfile=samfile)

    for read in samfile:
        eta.print_status(extra=read.qname,bam_pos=(read.rname,read.pos))
        if show_all or not read.is_unmapped:
            print '>%s\n%s' % (read.qname, read.seq)
    
    eta.done()
    samfile.close()

def usage():
    print """\
Usage: %s {-all} file.[sam,bam]

Ouputs the sequences of all mapped reads to FASTA format.

Options
-all   Output all reads, even if they are unmapped

""" % os.path.basename(sys.argv[0])

if __name__ == "__main__":
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[-1]):
        usage()
        sys.exit(1)

    if sys.argv[1] == '-all':
        bam2fasta(sys.argv[-1],True)
    else:
        bam2fasta(sys.argv[-1])
