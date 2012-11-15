#!/usr/bin/env python
## category Conversion
## desc Converts a file in VCF format to BED6
'''
Converts a file in VCF format to BED6

See: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
'''
import sys
import os
from ngsutils.support.ngs_utils import gzip_opener


def bed_fromvcf(fname, out=sys.stdout):
    with gzip_opener(fname) as f:
        _bed_fromvcf(f, out)


def _bed_fromvcf(fileobj, out=sys.stdout):
    for line in fileobj:
        if line[0] == '#':
            continue
        cols = line.strip('\n').split('\t')
        start = int(cols[1]) - 1
        end = start + len(cols[3])

        out.write('%s\t%s\t%s\t%s/%s\t%s\t+\n' % (cols[0], start, end, cols[3], cols[4], cols[5]))


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
