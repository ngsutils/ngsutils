#!/usr/bin/env python
## category General
## desc Filter annotations from a GTF file
'''
Filter annotations from a GTF file.
'''

import sys
import os

from ngsutils.support import gzip_reader


def usage():
    print __doc__
    print '''\
Usage: gtfutils filter {filters} filename.gtf

Possible filters:
    -chr str    Remove annotations from chromosomes with 'str' in the name
    
    -to-ucsc    Rename Ensembl-style chromosome names (1, 2, etc) to 
                UCSC/NCBI-style names (chr1, chr2, etc.) 
                
                (This will keep only chromosomes 1-22, X, Y, and MT for Human)

'''
    sys.exit(1)


def gtf_filter(fname, filters, out=sys.stdout):
    for line in gzip_reader(fname):
        cols = line.strip('\n').split('\t')
        good = True
        for filt in filters:
            cols = filt.process(cols)
            if not cols:
                good = False
                break

        if good:
            out.write('%s\n' % '\t'.join([str(x) for x in cols]))


class GTFFilter(object):
    def process(self, cols):
        raise NotImplementedError


class ChrSubstr(GTFFilter):
    def __init__(self, substr):
        self.substr = substr

    def process(self, cols):
        if self.substr in cols[0]:
            return None
        return cols


class ToUCSCChrom(GTFFilter):
    def __init__(self):
        self.valid = set('1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT'.split())

    def process(self, cols):
        if cols[0] in self.valid:
            out = cols[:]
            if cols[0] == 'MT':
                out[0] = 'chrM'
            else:
                out[0] = 'chr%s' % cols[0]
            return out
        return None


if __name__ == '__main__':
    fname = None
    last = None
    filters = []
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-chr':
            filters.append(ChrSubstr(arg))
            last = None
        elif arg == '-to-ucsc':
            filters.append(ToUCSCChrom())
        elif arg in ['-chr']:
            last = arg
        elif not fname and os.path.exists(arg):
            fname = arg
        else:
            usage()

    if not fname or not filters:
        usage()



    gtf_filter(fname, filters)
