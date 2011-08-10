#!/usr/bin/env python
'''
Converts a file in basecall format to BED3 format.
'''
import sys,os

def bed_frombasecall(fname):
    if fname == '-':
        f = sys.stdin
    else:
        f = open(fname)
    for line in f:
        if line[0] == '#':
            continue
        cols = line.strip('\n').split('\t')
        chrom = cols[0]
        pos = int(cols[1])-1
        sys.stdout.write('%s\t%s\t%s\n' % (chrom,pos,pos+1,))
    if fname != '-':
        f.close()

def usage():
    print __doc__
    print 'Usage: bedutils frombasecall basecall.txt (- for stdin)'
    sys.exit(1)
if __name__ == '__main__':
    fname = None
    for arg in sys.argv[1:]:
        if not fname and os.path.exists(arg):
            fname = arg
        elif arg == '-h':
            usage()
        else:
            usage()
            
    if not fname:
        usage()
        
    bed_frombasecall(fname)