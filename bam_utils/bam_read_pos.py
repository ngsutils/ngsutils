#!/usr/bin/env python

import sys,os
import pysam

def bam_read_pos(fname,readname):
    bamfile = pysam.Samfile(fname,"rb")
    count = 0
    pos = 0
    
    for read in bamfile:
        if read.qname == readname:
            pos = count
        count += 1
        
    bamfile.close()
    print "%s position: %s of %s" % (readname,pos,count)

def usage():
    print """\
Usage: %s bamfile read_name

""" % os.path.basename(sys.argv[0])

if __name__ == "__main__":
    if len(sys.argv) != 3 or not os.path.exists(sys.argv[1]):
        usage()
        sys.exit(1)

    bam_read_pos(sys.argv[1], sys.argv[2])

