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
import os
import collections
from ngsutils.support import gzip_reader
from ngsutils.fastq import convert_solexa_qual, convert_illumina_qual

QseqConversionResults = collections.namedtuple('QseqConversionResults', 'reads qcfailed lenfailed passed lengths')


class QseqRecord(collections.namedtuple('QseqRecord', 'machine run lane tile x y index read_num seq qual hidden_qcpass')):
    @property
    def qcpass(self):
        if self.hidden_qcpass == '0' or self.hidden_qcpass == 'QC':
            return False
        return True

    @property
    def name(self):
        name = '%s #%s/%s' % (':'.join([self.machine, self.run, self.lane, self.tile, self.x, self.y]), self.index, self.read_num)
        if not self.qcpass:
            name += ' QCFAIL'
        return name


def qseq_reader(fname=None, fileobj=None, quiet=False):
    if not fileobj:
        if not fname:
            raise ValueError('Must pass fname or fileobj!')

        for line in gzip_reader(fname, quiet=quiet):
            yield QseqRecord(*line.strip().split('\t')[:11])
    else:
        for line in fileobj:
            yield QseqRecord(*line.strip().split('\t')[:11])


def read_illumina_export(qseqreader, solexa_quals=False, min_length=0, trim=False, tag=None, qc_remove=True, out=sys.stdout):
    lengths = {}
    reads = 0
    passed = 0
    lenfailed = 0
    qcfailed = 0
    for line in qseqreader:
        reads += 1

        if qc_remove and not line.qcpass:
            qcfailed += 1
            continue

        if tag:
            name = "%s%s" % (tag, line.name)
        else:
            name = line.name

        seq = line.seq
        il_qual = line.qual

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

        out.write('@%s\n%s\n+\n%s\n' % (name, seq, qual))

        passed += 1
        if not len(seq) in lengths:
            lengths[len(seq)] = 1
        else:
            lengths[len(seq)] += 1

    return QseqConversionResults(reads, qcfailed, lenfailed, passed, lengths)


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
    results = read_illumina_export(qseq_reader(fname), solexa_quals, min_length, trim, tag, qc_remove)

    sys.stderr.write('Reads processed     : %s\n' % results.reads)
    sys.stderr.write('Failed QC filter    : %s\n' % results.qcfailed)
    sys.stderr.write('Failed length filter: %s\n' % results.lenfailed)
    sys.stderr.write('Passed              : %s\n' % results.passed)
    sys.stderr.write('\n')
    sys.stderr.write('Summary of read lengths\n')
    sortedar = []
    for length in results.lengths:
        sortedar.append((results.lengths[length], length))
    sortedar.sort()

    for count, length in sortedar:
        sys.stderr.write('Length %s: %s\n' % (length, count))
