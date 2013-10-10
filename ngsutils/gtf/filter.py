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

'''
    sys.exit(1)


def gtf_filter(fname, filters, out=sys.stdout):
    for line in gzip_reader(fname):
        cols = line.strip('\n').split('\t')
        good = True
        for filt in filters:
            if not filt.process(cols):
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
            return False
        return True


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
        elif arg in ['-chr']:
            last = arg
        elif not fname and os.path.exists(arg):
            fname = arg
        else:
            usage()

    if not fname or not filters:
        usage()



    gtf_filter(fname, filters)
