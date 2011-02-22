#!/usr/bin/env python
'''
Ouputs the names and positions of reads in the BAM file.
'''
import sys,os
from support.eta import ETA
import pysam

def bam_read_names(fname,mapped=False,unmapped=False):
    bamfile = pysam.Samfile(fname,"rb")
    eta = ETA(0,bamfile=bamfile)

    for read in bamfile:
        eta.print_status(extra=read.qname,bam_pos=(read.rname,read.pos))
        if mapped and not read.is_unmapped:
            print '%s\t%s\t%s\t%s' % (read.qname, bamfile.getrname(read.rname), read.pos, ''.join(['%s%s' % (length,bam_cigar[op]) for op,length in read.cigar]))
        elif unmapped and read.is_unmapped:
            print '%s\t*\t0\t\n' % (read.qname)
    
    eta.done()
    bamfile.close()

def usage():
    print __doc__
    print """\
Usage: %s {-mapped} bamfile

Options
-mapped     Output only mapped reads
-unmapped   Output only unmapped reads

""" % os.path.basename(sys.argv[0])

if __name__ == "__main__":
    mapped = False
    unmapped = False
    fname = None
    for arg in sys.argv[1:]:
        if arg == '-unmapped':
            unmapped = True
        elif arg == '-mapped':
            mapped = True
        elif os.path.exists(arg):
            fname = arg
            
    if not fname
        usage()
        sys.exit(1)

    if not unmapped and not mapped:
        bam_read_names(fname,True,True)
    else:
        bam_read_names(fname,mapped,unmapped)

