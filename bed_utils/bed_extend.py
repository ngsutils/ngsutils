#!/usr/bin/env python
'''
Extends BED regions to be a fixed size
'''

import os,sys
sys.path.append(os.path.join(os.path.dirname(__file__),"..","utils")) #eta
import ngs_utils

def usage():
    print __doc__
    print """\
Usage: %s SIZE bedfile

SIZE is what the total size of the region should be.  The size of the region 
will be extended or reduced to make the total length exactly SIZE.

However, if the first character of SIZE is '+', then the region is extended 
SIZE bases, regardless of how long it is to start with.

""" % os.path.basename(sys.argv[0])
    sys.exit(1)

def bed_extend(fname,size,relative=False):
    with ngs_utils.gzip_opener(fname) as f:
        for line in f:
            chrom,start,end,name,score,strand = line.strip().split('\t')

            if strand == '+':
                if relative:
                    end = int(end) + size
                else:
                    end = int(start) + size
            else:
                if relative:
                    start = int(start) - size
                else:
                    start = int(end) - size
        
            if start < 0:
                start = 0
        
            sys.stdout.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom,start,end,name,score,strand))

if __name__ == '__main__':
    fname = None
    size = None
    relative = False
    last = None
    
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if not size:
            if arg[0] == '+':
                relative = True
                size = int(arg[1:])
            else:
                size = int(arg)
        elif not fname and os.path.exists(arg):
            fname = arg
        else:
            print "Unknown option: %s" % arg
            usage()
    
    if not fname or not size:
        usage()
    
    bed_extend(fname,size,relative)