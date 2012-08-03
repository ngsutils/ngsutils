#!/usr/bin/env python
## category Conversion
## desc Converts a file in VCF format to BED6
'''
Converts a file in VCF format to BED6
'''
import sys
import os
import gzip


def bed_fromvcf(fname):
    if fname == '-':
        f = sys.stdin
    elif fname[-4:] == '.bgz' or fname[-3:] == '.gz':
        f = gzip.open(fname)
    else:
        f = open(fname)
    for line in f:
        if line[0] == '#':
            continue
        cols = line.strip('\n').split('\t')
        start = int(cols[1]) - 1
        end = start + len(cols[3])

        sys.stdout.write('%s\t%s\t%s\t%s/%s\t%s\t+\n' % (cols[0], start, end, cols[3], cols[4], cols[5]))
    if fname != '-':
        f.close()


def usage():
    print __doc__
    print 'Usage: bedutils fromvcf file.vcf{.bgz} (- for stdin)'
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

    bed_fromvcf(fname)
