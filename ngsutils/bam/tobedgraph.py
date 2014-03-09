#!/usr/bin/env python
## category Conversion
## desc Convert BAM coverage to bedGraph (for visualization)
'''
Convert BAM coverage to bedGraph based on pileup counts

Takes a BAM file and produces a bedGraph file.  This can optionally
normalize the counts by a given factor.

 See: http://genome.ucsc.edu/goldenPath/help/bedgraph.html
      http://genome.ucsc.edu/goldenPath/help/bigWig.html
'''

import sys
import os
from array import array
from ngsutils.bam import bam_pileup_iter, cigar_tostr, read_cigar_at_pos, bam_iter
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

        self.pos_counts = array('L')

        self.cur_tid = None
        self.cur_pos = None
        self.cur_chrom = None

        self._last_val = 0
        self._last_start = None
        self._last_pos = None

        self.out = out
        self.strand = strand


    def _clear_pos_counts(self, readpos):
        if self.cur_pos == readpos:
            return

        i = 0
        while i < len(self.pos_counts) and i < readpos - self.cur_pos:
            self.write(self.cur_pos + i, self.pos_counts[i])
            i += 1

        self.pos_counts = self.pos_counts[i:]
        self.cur_pos = readpos
        


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

    def get_counts(self, bam, quiet=False):
        for read in bam_iter(bam, quiet=quiet):
            if read.is_unmapped:
                continue
            if self.strand:
                if self.strand == '+' and read.is_reverse:
                    continue
                elif self.strand == '-' and not read.is_reverse:
                    continue

            if self.cur_tid is None or read.tid != self.cur_tid:
                if self.pos_counts:
                    for i, count in enumerate(self.pos_counts):
                        self.write(self.cur_pos + i, count)
                    self.flush()

                self.pos_counts = array('L')
                self.read_cache = set()

                self.cur_pos = 0
                self.cur_tid = read.tid
                self.cur_chrom = bam.references[read.tid]

            self._clear_pos_counts(read.pos)
            self._add_read(read)

        for i, count in enumerate(self.pos_counts):
            self.write(self.cur_pos + i, count)
        self.flush()


    def incr_count(self, start, end=None):
        if not end:
            end = start
       
        if len(self.pos_counts) < end:
            self.pos_counts.extend([0,] * (end+1-self.cur_pos+len(self.pos_counts)))

        for pos in xrange(start, end):
            self.pos_counts[pos - self.cur_pos] += 1

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
        if self._last_val:
            self.out.write('%s\t%s\t%s\t%s\n' % (self.cur_chrom, self._last_start, self._last_pos+1, (self._last_val * self.normalization_factor)))

        self._last_val = 0
        self._last_start = None
        self._last_pos = None


def bam_tobedgraph(bamfile, strand=None, normalize=None, out=sys.stdout):
    if normalize is None:
        normalize = 1
    counter = BamCounter(normalize, strand, out)
    counter.get_counts(bamfile)


def usage():
    print __doc__
    print """\
Usage: bamutils tobedgraph [-plus | -minus] {-norm N} bamfile

Options:
    -plus             only count reads on the plus strand
                      (default: count all reads)
    -minus            only count reads on the minus strand

    -norm VAL         the count at every position is calculated as:
                      floor(count * VAL).
"""
    sys.exit(1)

if __name__ == "__main__":
    bam = None
    strand = None
    norm = 1 

    last = None
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if last == '-norm':
            norm = float(arg)
            last = None
        elif arg in ['-norm']:
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
    bam_tobedgraph(bamfile, strand, norm, out=sys.stdout)
    bamfile.close()
