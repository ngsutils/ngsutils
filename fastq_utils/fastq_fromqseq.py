#!/usr/bin/env python
## category Conversion
## desc Converts Illumina qseq (export/sorted) files to FASTQ
'''
Converts an Illumina export.txt file to a FASTQ formatted file.
This mainly consists of converting the quality scores to Phred format.
In addtion, this can also remove reads which have failed the QC filter and
trim away sequence that have a trailing 'B' for quality (Read Segment Quality
Control Indicator).
Based on: http://www.cassj.co.uk/blog/?p=490
          http://en.wikipedia.org/wiki/FASTQ_format
'''

import sys
import math
import os
import gzip

from support.eta import ETA


def read_illumina_export(filename, solexa_quals=True, min_length=0, trim=False, tag=None, qc_remove=True):
    if filename[-3:] == '.gz':
        f = gzip.open(filename)
    else:
        f = open(filename)

    fsize = os.stat(filename).st_size
    eta = ETA(fsize, fileobj=f, modulo=5000)

    lengths = {}
    reads = 0
    passed = 0
    lenfailed = 0
    qcfailed = 0
    for line in f:
        reads += 1
        cols = line.strip().split('\t')
        name = '%s #%s/%s' % (':'.join(cols[0:6]), cols[6], cols[7])
        if cols[10] == 'QC':
            if qc_remove:
                qcfailed += 1
                continue
            else:
                name = '%s QCFAIL' % name

        if tag:
            name = "%s_%s" % (tag, name)

        seq = cols[8]
        il_qual = cols[9]

        if trim:
            # Trim trailing 'B's
            while il_qual and il_qual[-1] == 'B':
                il_qual = il_qual[:-1]
                seq = seq[:-1]

            if not il_qual or len(seq) < min_length:
                lenfailed += 1
                continue

        if solexa_quals:
            qual = convert_solexa_qual(il_qual)
        else:
            qual = convert_illumina_qual(il_qual)

        sys.stdout.write('@%s\n%s\n+\n%s\n' % (name, seq, qual))

        passed += 1
        if not len(seq) in lengths:
            lengths[len(seq)] = 1
        else:
            lengths[len(seq)] += 1

        if eta:
            eta.print_status()

    eta.done()
    f.close()

    sys.stderr.write('Reads processed     : %s\n' % reads)
    sys.stderr.write('Failed QC filter    : %s\n' % qcfailed)
    if trim:
        sys.stderr.write('Failed length filter: %s\n' % lenfailed)
        sys.stderr.write('Passed length filter: %s\n' % passed)
    else:
        sys.stderr.write('Passed              : %s\n' % passed)

    sys.stderr.write('Summary of read lengths\n')
    sortedar = []
    for length in lengths:
        sortedar.append((lengths[length], length))
    sortedar.sort()

    for count, length in sortedar:
        sys.stderr.write('Length %s: %s\n' % (length, count))


def convert_solexa_qual(qual):
    '''
    Illumina char: QSolexa + 64  (note: this is for very old samples)
    Phred char: QPhred + 33

    QPhred = -10 * log10 (1/error)
    QSolexa = -10 * log10 (error/(1-error))

    QPhred = 10 * log10 (10 ^ (QSolexa/10) + 1)

    '''

    rv = []
    for q in qual:
        val = ord(q) - 64
        qp = 10 * math.log10(10 ** (val / 10) + 1)
        rv.append(chr(qp + 33))
    return ''.join(rv)


def convert_illumina_qual(qual):
    '''
    Illumina char: QPhred + 64
    Phred char: QPhred + 33
    '''

    return ''.join([chr(ord(q) - 31) for q in qual])


def usage():
    print __doc__
    print """\
Usage: fastqutils fromqseq {opts} export.txt

Options:
  -solexa     file is pre-1.3 format (uses Solexa calculation for quality)
  -tag tag    Prefix the read names with this tag (such as the sample name)
  -trim       perform quality control indicator (trailing B) trimming
  -min N      the minimum allowed length for a read, post B-trimming
  -noqc       Don't remove reads that failed QC (for matching paired end data)
"""
    sys.exit(1)

if __name__ == '__main__':
    solexa_quals = False
    fname = None
    min_length = 0
    trim = False
    last = None
    tag = None
    qc_remove = True

    for arg in sys.argv[1:]:
        if arg in ['-h', '--help']:
            usage()
            sys.exit(1)

        if last == '-min':
            min_length = int(arg)
            last = None
        elif arg == '-min':
            last = arg
        elif last == '-tag':
            tag = arg
            last = None
        elif arg == '-h':
            usage()
        elif arg == '-tag':
            last = arg
        elif arg == '-trim':
            trim = True
        elif arg == '-solexa':
            solexa_quals = True
        elif arg == '-illumina':
            solexa_quals = False
        elif arg == '-noqc':
            qc_remove = False
        elif os.path.exists(arg):
            fname = arg

    if not fname or solexa_quals is None:
        usage()

    sys.stderr.write("Converting file: %s\n(using %s scaling)\n%s" % (fname, 'Solexa' if solexa_quals else 'Illumina', '(min-length %d)\n' % min_length if min_length else ''))
    read_illumina_export(fname, solexa_quals, min_length, trim, tag, qc_remove)
