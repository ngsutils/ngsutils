#!/usr/bin/env python
## category Conversion
## desc Convert BAM coverage to bedGraph (for visualization)
'''
Convert BAM coverage to bedGraph based on read depth. This will take into
account gaps in RNAseq alignments and not display any coverage across introns.

This can optionally normalize the counts by a given factor or display only
coverage for a specific strand.

 See: http://genome.ucsc.edu/goldenPath/help/bedgraph.html
      http://genome.ucsc.edu/goldenPath/help/bigWig.html
'''

import sys
import os
from array import array
from ngsutils.bam import bam_iter
import pysam


def write_bedgraph(chrom, start, end, count, normalize=None, out=sys.stdout):
    if start and end and count:
        if normalize:
            out.write('%s\t%s\t%s\t%s\n' % (chrom, start, end, normalize * count))
        else:
            out.write('%s\t%s\t%s\t%s\n' % (chrom, start, end, count))


class BamCounter(object):
    def __init__(self, normalization_factor=1, strand=None, out=sys.stdout):
        self.normalization_factor = normalization_factor

        self.pos_counts = None

        self.cur_tid = None
        self.cur_chrom = None

        self._last_val = 0
        self._last_start = None
        self._last_pos = None

        self.out = out
        self.strand = strand


    def _add_read(self, read):
        refpos = read.pos
        for op, size in read.cigar:
            if op == 0:
                self.incr_count(refpos, refpos + size)
                refpos += size
            elif op == 1:
                pass
            elif op == 2:
                refpos += size
            elif op == 3:
                refpos += size
            pass


    def get_counts(self, bam, ref=None, start=None, end=None, quiet=False):
        for read in bam_iter(bam, ref=ref, start=start, end=end, quiet=quiet, callback=lambda x: '%s:%s (%s)' % (self.cur_chrom, x.pos, len(self.pos_counts) if self.pos_counts else 0)):
            if read.is_unmapped:
                continue
            if self.strand:
                if self.strand == '+' and read.is_reverse:
                    continue
                elif self.strand == '-' and not read.is_reverse:
                    continue

            if self.cur_tid is None or read.tid != self.cur_tid:
                if self.pos_counts:
                    self.flush()

                self.pos_counts = array('I', [0,] * bam.lengths[read.tid])

                self.cur_tid = read.tid
                self.cur_chrom = bam.references[read.tid]

            self._add_read(read)

        self.flush()

    def incr_count(self, start, end=None):
        if not end:
            end = start

        for pos in xrange(start, end):
            self.pos_counts[pos] += 1

    def write(self, pos, val):
        if val == self._last_val:
            self._last_pos = pos
        else:
            if self._last_val:
                self.out.write('%s\t%s\t%s\t%s\n' % (self.cur_chrom, self._last_start, self._last_pos+1, (self._last_val * self.normalization_factor)))
            self._last_val = val
            self._last_start = pos
            self._last_pos = pos

    def flush(self):
        for i, count in enumerate(self.pos_counts):
            self.write(i, count)

        if self._last_val:
            self.out.write('%s\t%s\t%s\t%s\n' % (self.cur_chrom, self._last_start, self._last_pos+1, (self._last_val * self.normalization_factor)))

        self._last_val = 0
        self._last_start = None
        self._last_pos = None


def bam_tobedgraph(bamfile, strand=None, normalize=None, ref=None, start=None, end=None, out=sys.stdout):
    if normalize is None:
        normalize = 1
    counter = BamCounter(normalize, strand, out)
    counter.get_counts(bamfile, ref, start, end)


def usage():
    print __doc__
    print """\
Usage: bamutils tobedgraph {opts} bamfile

Options:
    -plus             Only count reads on the plus strand
                      (default: count all reads)
    -minus            Only count reads on the minus strand

    -norm VAL         The count at every position is calculated as:
                      floor(count * VAL).

    -ref name         Only count reads mapping to this reference (chrom)

    -region chr:start-end    Count reads mapping to this genome region
                             (start is 1-based)
"""
    sys.exit(1)

if __name__ == "__main__":
    bam = None
    strand = None
    norm = 1
    ref = None
    start = None
    end = None

    last = None
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if last == '-norm':
            norm = float(arg)
            last = None
        elif last == '-ref':
            ref = arg
            last = None
        elif last == '-region':
            ref, se = arg.split(':')
            start, end = [int(x) for x in se.split('-')]
            start = start - 1
            last = None
        elif arg in ['-norm', '-ref', '-region']:
            last = arg
        elif arg == '-plus':
            strand = '+'
        elif arg == '-minus':
            strand = '-'
        elif not bam and os.path.exists(arg) and os.path.exists('%s.bai' % arg):
            bam = arg
        else:
            print "Unknown option or missing index: %s" % arg
            usage()

    if not bam:
        usage()

    bamfile = pysam.Samfile(bam, "rb")
    bam_tobedgraph(bamfile, strand, norm, ref, start, end, out=sys.stdout)
    bamfile.close()
