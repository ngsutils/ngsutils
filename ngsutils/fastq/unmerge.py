#!/usr/bin/env python
## category General
## desc Unmerged paired FASTQ files into two (or more) files
'''
Splits a combined FASTQ file into two (or more) FASTQ files.

Files will be written as "out_template.N.fastq" where N is the
pair number.
'''

import os
import sys
import gzip

from ngsutils.fastq import FASTQ


def fastq_unmerge(combined_fname, out_template, gz=False):
    outs = []
    if gz:
        outs.append(gzip.open('%s.1.fastq.gz' % out_template, 'w'))
    else:
        outs.append(open('%s.1.fastq' % out_template, 'w'))

    outidx = 1

    last_read = None
    fq = FASTQ(combined_fname)
    for read in fq.fetch():
        if last_read and last_read.name == read.name:
            outidx += 1
            if len(outs) < outidx:
                if gz:
                    outs.append(gzip.open('%s.%s.fastq.gz' % (out_template, outidx), 'w'))
                else:
                    outs.append(open('%s.%s.fastq' % (out_template, outidx), 'w'))
            read.write(outs[outidx - 1])
        else:
            outidx = 1
            read.write(outs[0])

        last_read = read

    fq.close()
    for out in outs:
        out.close()


def usage():
    print __doc__
    print """Usage: fastqutils unmerge {options} combined.fastq out_template

Options:
  -gz    gzip compress the output files
"""
    sys.exit(1)

if __name__ == '__main__':
    combined_fname = None
    out_template = None
    gz = False

    for arg in sys.argv[1:]:
        if arg == '-gz':
            gz = True
        elif not combined_fname and os.path.exists(arg):
            combined_fname = arg
        elif not out_template:
            out_template = arg

    if not combined_fname or not out_template:
        usage()

    fastq_unmerge(combined_fname, out_template, gz=gz)
