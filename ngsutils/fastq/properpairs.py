#!/usr/bin/env python
## category General
## desc Find properly paired reads (when fragments are filtered separately)
'''
When you filter paired FASTQ files separately, fragments my be filtered or
removed from either file. This program will then find all of the properly
paired reads. Valid pairs will be written to new output files that may be
optionally gzip compressed.

Caution: this can use up to 2X the disk space of each FASTQ file!
'''

import os
import sys
import gzip

import tempfile

from ngsutils.fastq import FASTQ

import ngsutils.fastq.sort


def find_fastq_pairs(fq1, fq2, out1, out2, tmpdir=None, quiet=False):
    tmp1 = tempfile.NamedTemporaryFile(delete=False, prefix='.tmp', suffix='.gz', dir=tmpdir if tmpdir else os.path.dirname(fq1.fname))
    tmp1_fname = tmp1.name
    tmp1_out = gzip.GzipFile(fileobj=tmp1)

    ngsutils.fastq.sort.fastq_sort(fq1, out=tmp1_out, tmpdir=tmpdir if tmpdir else os.path.dirname(fq1.fname))
    tmp1_out.close()
    tmp1.close()

    tmp2 = tempfile.NamedTemporaryFile(delete=False, prefix='.tmp', suffix='.gz', dir=tmpdir if tmpdir else os.path.dirname(fq2.fname))
    tmp2_fname = tmp2.name
    tmp2_out = gzip.GzipFile(fileobj=tmp2)

    ngsutils.fastq.sort.fastq_sort(fq2, out=tmp2_out, tmpdir=tmpdir if tmpdir else os.path.dirname(fq2.fname))
    tmp2_out.close()
    tmp2.close()

    sys.stderr.write('Finding properly paired FASTQ reads...\n')

    fq_tmp1 = FASTQ(tmp1_fname)
    fq_tmp2 = FASTQ(tmp2_fname)

    reader1 = fq_tmp1.fetch(quiet=quiet)
    reader2 = fq_tmp2.fetch(quiet=True)

    read1 = reader1.next()
    read2 = reader2.next()

    pairs = 0
    discarded_1 = 0
    discarded_2 = 0

    while read1 and read2:
        if read1.name == read2.name:
            read1.write(out1)
            read2.write(out2)

            try:
                read1 = reader1.next()
                read2 = reader2.next()
            except StopIteration:
                break

            pairs += 1
        elif read1.name < read2.name:
            discarded_1 += 1
            try:
                read1 = reader1.next()
            except StopIteration:
                break
        else:
            discarded_2 += 1
            try:
                read2 = reader2.next()
            except StopIteration:
                break

    fq_tmp1.close()
    fq_tmp2.close()

    os.unlink(tmp1_fname)
    os.unlink(tmp2_fname)

    return pairs, discarded_1, discarded_2

def usage(msg=""):
    if msg:
        print '%s\n' % msg
    print """Usage: fastqutils properpairs {opts} filename1.fastq{.gz} filename2.fastq{.gz} output1 output2

Options:
  -f       Force overwriting output file (if it exists)
  -z       Output files should be gzip compressed
  -t dir   Use {dir} for temporary files
"""
    sys.exit(1)

if __name__ == '__main__':
    fqname1 = None
    fqname2 = None
    outname1 = None
    outname2 = None
    tmpdir = None
    force = False

    gz = False

    last = None

    for arg in sys.argv[1:]:
        if arg == '-z':
            gz = True
        elif arg == '-f':
            force = True
        elif last == '-t':
            if os.path.exists(arg) and os.path.isdir(arg):
                tmpdir = arg
            else:
                usage('%s is not a valid temp-directory!' % arg)

            last = None
        elif arg in ['-t']:
            last = arg
        elif not fqname1:
            if not os.path.exists(arg):
                usage("File %s doesn't exist!" % arg)
            fqname1 = arg
        elif not fqname2:
            if not os.path.exists(arg):
                usage("File %s doesn't exist!" % arg)
            fqname2 = arg
        elif not outname1:
            outname1 = arg
        elif not outname2:
            outname2 = arg

    if not fqname1 or not fqname2 or not outname1 or not outname2:
        usage()

    if not force:
        for fname in [outname1, outname2]:
            if os.path.exists(fname):
                usage("File %s exists!" % fname)

    fq1 = FASTQ(fqname1)
    fq2 = FASTQ(fqname2)

    if gz:
        out1 = gzip.open(outname1, 'w')
        out2 = gzip.open(outname2, 'w')
    else:
        out1 = open(outname1, 'w')
        out2 = open(outname2, 'w')

    paired, discard_1, discard_2 = find_fastq_pairs(fq1, fq2, out1, out2, tmpdir)

    print "Proper pairs: %s" % paired
    print "Discarded 1 : %s" % discard_1
    print "Discarded 2 : %s" % discard_2

    fq1.close()
    fq2.close()
    out1.close()
    out2.close()
