#!/usr/bin/env python
'''
Writes out the read names present in the FASTQ file.
'''

import os,sys

from fastq_utils import read_fastq

def usage():
    print __doc__
    print """Usage: %s filename.fastq{.gz}""" % os.path.basename(sys.argv[0])
    sys.exit(1)

if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if os.path.exists(arg):
            fname = arg

    if not fname:
        usage()
        
    for name,seq,qual in read_fastq(fname):
        sys.stdout.write('%s\n' % name.split()[0][1:])
        
