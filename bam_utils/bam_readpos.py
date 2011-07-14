#!/usr/bin/env python
'''
Finds the realtive position of a read or reference:position in a BAM file.
This is returned as a line number out of the total number of lines in the 
file.
'''



import sys,os
import pysam

def bam_readpos(fname,readname=None,ref=None,pos=None):
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
        print "%s position: %s of %s (%.2f%%)" % (readname,i,count,float(i)*100/count)
    else:
        print "%s:%s position: %s of %s (%.2f%%)" % (ref,pos,i,count,float(i)*100/count)
        

def usage():
    print __doc__
    print """
Usage: bamutils readpos {-read read_name} {-pos chr:pos}  bamfile
"""

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
        bam_readpos(fname,read=readname)
    elif ref and pos and fname:
        bam_readpos(fname,ref=ref,pos=pos)
    else:
        usage()
        sys.exit(1)
