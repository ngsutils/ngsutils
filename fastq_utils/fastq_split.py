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

from fastq_utils import read_fastq, is_paired_fastq


def fastq_split(fname, outbase, chunks, ignore_pairs=False, gz=False, count_fname=None):
    i = 0

    if ignore_pairs:
        is_paired = False
    else:
        is_paired = is_paired_fastq(fname)

    outs = []
    fnames = []
    for i in xrange(chunks):
        if gz:
            fn = '%s.%s.fastq.gz' % (outbase, i + 1)
            tmp = os.path.join(os.path.dirname(fn), '.tmp.%s' % os.path.basename(fn))
            fnames.append((tmp, fn))

            sys.stderr.write('Output file: %s\n' % fn)
            outs.append(gzip.open(tmp, 'w'))
        else:
            fn = '%s.%s.fastq' % (outbase, i + 1)
            tmp = os.path.join(os.path.dirname(fn), '.tmp.%s' % os.path.basename(fn))
            fnames.append((tmp, fn))

            sys.stderr.write('Output file: %s\n' % fn)
            outs.append(open(tmp, 'w'))

    i = 0
    last_name = None
    count = 0

    for name, seq, qual in read_fastq(fname):
        count += 1
        sn = name.split()[0]
        if not is_paired:
            i += 1
        elif sn != last_name:
            i += 1

        if i >= len(outs):
            i = 0

        last_name = sn

        outs[i].write('%s\n%s\n+\n%s\n' % (name, seq, qual))

    for out in outs:
        out.close()

    for tmp, fname in fnames:
        os.rename(tmp, fname)

    if count_fname:
        with open(count_fname, 'w') as f:
            f.write(count)


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

  -count fname     Write the number of reads in the FASTQ file to a text file

  -gz              gzip compress the output FASTQ files

"""
    sys.exit(1)

if __name__ == '__main__':
    fname = None
    outtemplate = None
    chunks = 0
    ignore_pairs = False
    gz = False
    countfname = None
    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if arg == '-ignorepaired':
            ignore_pairs = True
        elif arg == '-gz':
            gz = True
        elif last == '-count':
            countfname = arg
            last = None
        elif arg in ['-count']:
            last = arg
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

    fastq_split(fname, outtemplate, chunks, ignore_pairs, gz, countfname)
