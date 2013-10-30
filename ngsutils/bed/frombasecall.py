#!/usr/bin/env python
## category Conversion
## desc Converts a file in basecall format to BED3 format
'''
Converts a file in basecall format to BED3 format
'''
import sys
import os


def bed_frombasecall(fname):
    if fname == '-':
        f = sys.stdin
    else:
        f = open(fname)

    _bed_frombasecall(f)

    if fname != '-':
        f.close()


def _bed_frombasecall(fileobj, out=sys.stdout):
    for line in fileobj:
        line = line.strip()
        if not line or line[0] == '#':
            continue
        cols = line.split('\t')
        chrom = cols[0]
        pos = int(cols[1]) - 1
        out.write('%s\t%s\t%s\n' % (chrom, pos, pos + 1,))


def usage():
    print __doc__
    print 'Usage: bedutils frombasecall basecall.txt (- for stdin)'
    sys.exit(1)
if __name__ == '__main__':
    fname = None
    for arg in sys.argv[1:]:
        if not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg
        elif arg == '-h':
            usage()
        else:
            usage()

    if not fname:
        usage()

    bed_frombasecall(fname)
