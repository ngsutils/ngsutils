#!/usr/bin/env python
## category General
## desc Adds a prefix or suffix to the read names in a FASTQ file
'''
Adds a prefix or suffix to the read names for all reads in the file

Writes output to stdout
'''

import os
import sys

from fastq_utils import read_fastq


def fastq_tag(fname, prefix, suffix):
    for name, seq, qual in read_fastq(fname):
        spl = name[1:].split(None, 1)
        nname = ''
        if len(spl) > 1:
            desc = spl[1]
        else:
            desc = None

        if prefix and suffix:
            nname = '%s%s%s' % (prefix, spl[0], suffix)
        elif prefix:
            nname = '%s%s' % (prefix, spl[0])
        elif suffix:
            nname = '%s%s' % (spl[0], suffix)

        if desc:
            nname = '%s %s' % (nname, desc)

        sys.stdout.write('@%s\n%s\n+\n%s\n' % (nname, seq, qual))


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

    fastq_tag(fname, prefix, suffix)
