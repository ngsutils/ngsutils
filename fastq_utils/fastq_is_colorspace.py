#!/usr/bin/env python
'''
Scans reads in a FASTQ file to see if they are in colorspace.
'''

import os,sys

from fastq_utils import read_fastq

def fastq_is_colorspace(fname,quiet=False):
    '''
    This works by scanning the first 10 reads that have sequences (aren't Ns 
    or 4s). If there are any colorspace values, the entire file is called as 
    colorspace.
    
    It's a bit overkill...
    '''
    seqs_checked=0
    is_colorspace = False

    valid_basespace="atcgATCG"
    valid_colorspace="0123456"

    for name,seq,qual in read_fastq(fname,quiet=True):
        if seqs_checked > 10:
            break
        checked = False
        for base in seq:
            if base in valid_colorspace:
                is_colorspace = True
                checked = True
            elif base in valid_basespace:
                checked = True
        if checked:
            seqs_checked += 1

    if not quiet:
        if is_colorspace:
            print "colorspace"
        else:
            print "basespace"
    
    return is_colorspace

def usage():
    print """Usage: %s filename.fastq{.gz}""" % os.path.basename(sys.argv[0])
    sys.exit(1)

if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if os.path.exists(arg):
            fname = arg

    if not fname:
        usage()
        
    fastq_is_colorspace(fname)
