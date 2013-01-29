#!/usr/bin/env python
## category General
## desc Adds a prefix or suffix to the read names in a FASTQ file
'''
Adds a prefix or suffix to the read names for all reads in the file

Writes output to stdout
'''

import os
import sys

from ngsutils.fastq import FASTQ


def fastq_tag(fastq, prefix='', suffix='', out=sys.stdout, quiet=False):
    if not prefix and not suffix:
        raise ValueError('Must pass at least one of: prefix, suffix.')

    for read in fastq.fetch(quiet):
        if prefix and suffix:
            nname = '%s%s%s' % (prefix, read.name, suffix)
        elif prefix:
            nname = '%s%s' % (prefix, read.name)
        elif suffix:
            nname = '%s%s' % (read.name, suffix)

        read.clone(name=nname).write(out)


def usage():
    print __doc__
    print "fastqutils tag {-suf suffix} {-pre prefix} filename.fastq"
    sys.exit(1)

if __name__ == '__main__':
    prefix = None
    suffix = None
    fname = None
    last = None
    for arg in sys.argv[1:]:
        if last == '-pre':
            prefix = arg
            last = None
        elif last == '-suf':
            suffix = arg
            last = None
        elif arg in ['-pre', '-suf']:
            last = arg
        elif not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg

    if not (prefix or suffix) or not fname:
        usage()

    fq = FASTQ(fname)
    fastq_tag(fq, prefix, suffix)
    fq.close()
