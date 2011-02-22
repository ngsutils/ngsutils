#!/usr/bin/env python

import sys,os
import pysam

def bam_read_pos(fname,readname=None,ref=None,pos=None):
    bamfile = pysam.Samfile(fname,"rb")
    count = 0
    i = 0
    
    for read in bamfile:
        if readname and read.qname == readname:
            i = count
        elif bamfile.getrname(read.rname) == ref:
            if not i and read.pos >= pos:
                i = count
        count += 1
        
    bamfile.close()

    if readname:
        print "%s position: %s of %s (%s%%)" % (readname,i,count,float(i)*100/count)
    else:
        print "%s:%s position: %s of %s (%s%%)" % (ref,pos,i,count,float(i)*100/count)
        

def usage():
    print """\
Usage: %s {-read read_name} {-pos chr:pos}  bamfile

""" % os.path.basename(sys.argv[0])

if __name__ == "__main__":
    readname = None
    ref = None
    pos = None
    fname = None
    
    last = None
    for arg in sys.argv[1:]:
        if last == '-read':
            readname = arg
            last = None
        elif last == '-pos':
            ref = arg.split(':')[0]
            pos = int(arg.split(':')[1])
            last = None
        elif arg in ['-pos','-read']:
            last = arg
        elif os.path.exists(arg):
            fname = arg
    
    if readname and fname:
        bam_read_pos(fname,read=readname)
    elif ref and pos and fname:
        bam_read_pos(fname,ref=ref,pos=pos)
    else:
        usage()
        sys.exit(1)
