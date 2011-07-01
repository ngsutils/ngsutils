#!/usr/bin/env python
'''
Converts a file in basecall format to BED3 format.
'''
import sys,os

def basecall_to_bed(fname):
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
    print 'Usage: %s basecall.txt (- for stdin)' % os.path.basename(sys.argv[0])
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
        
    basecall_to_bed(fname)