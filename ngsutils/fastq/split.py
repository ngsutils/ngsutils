#!/usr/bin/env python
## category General
## desc Splits a FASTQ file into N chunks
'''
Splits a FASTQ file into multiple smaller files

Output is a set of gzip compressed FASTQ files
'''

import os
import sys
import gzip

from ngsutils.fastq import FASTQ


def fastq_split(fname, outbase, chunks, ignore_pairs=False, gz=False, count_fname=None, quiet=False):
    fastq = FASTQ(fname)

    if ignore_pairs:
        is_paired = False
    else:
        is_paired = fastq.is_paired

    outs = []
    fnames = []
    for i in xrange(chunks):
        if gz:
            fn = '%s.%s.fastq.gz' % (outbase, i + 1)
            tmp = os.path.join(os.path.dirname(fn), '.tmp.%s' % os.path.basename(fn))
            fnames.append((tmp, fn))

            if not quiet:
                sys.stderr.write('Output file: %s\n' % fn)
            outs.append(gzip.open(tmp, 'w'))
        else:
            fn = '%s.%s.fastq' % (outbase, i + 1)
            tmp = os.path.join(os.path.dirname(fn), '.tmp.%s' % os.path.basename(fn))
            fnames.append((tmp, fn))

            if not quiet:
                sys.stderr.write('Output file: %s\n' % fn)
            outs.append(open(tmp, 'w'))

    i = chunks
    last_name = None

    for read in fastq.fetch(quiet=quiet):
        if not is_paired:
            i += 1
        elif read.name != last_name:
            i += 1

        if i >= len(outs):
            i = 0

        last_name = read.name

        read.write(outs[i])

    for out in outs:
        out.close()

    fastq.close()

    for tmp, fname in fnames:
        os.rename(tmp, fname)


def usage(msg=None):
    if msg:
        print msg

    print __doc__
    print """\
Usage: fastqutils split {opts} filename.fastq{.gz} out_template num_chunks

Options:
  -ignorepaired    Normally for paired-end samples, each read of the pair is
                   written to the same file. With this option, paired reads
                   can be written to any file. This is useful for splitting a
                   paired FASTQ file back into separate files for each
                   fragment.

  -gz              gzip compress the output FASTQ files

"""
    sys.exit(1)

if __name__ == '__main__':
    fname = None
    outtemplate = None
    chunks = 0
    ignore_pairs = False
    gz = False
    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if arg == '-ignorepaired':
            ignore_pairs = True
        elif arg == '-gz':
            gz = True
        elif not fname:
            if not os.path.exists(arg):
                usage("Missing file: %s" % arg)
            fname = arg
        elif not outtemplate:
            outtemplate = arg
        else:
            chunks = int(arg)

    if not fname or not chunks or not outtemplate:
        usage()

    fastq_split(fname, outtemplate, chunks, ignore_pairs, gz)
