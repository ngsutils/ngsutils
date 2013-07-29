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
from ngsutils.bam import bam_pileup_iter, cigar_tostr, read_cigar_at_pos
import pysam


def write_bedgraph(chrom, start, end, count, normalize=None, out=sys.stdout):
    if start and end and count:
        if normalize:
            out.write('%s\t%s\t%s\t%s\n' % (chrom, start, end, normalize * count))
        else:
            out.write('%s\t%s\t%s\t%s\n' % (chrom, start, end, count))


def bam_tobedgraph(bamfile, strand=None, normalize=None, nogaps=False, out=sys.stdout):
    last_chrom = None
    last_count = 0
    last_start = None
    last_end = None

    for pileup in bam_pileup_iter(bamfile):
        # sys.stdin.readline()
        chrom = bamfile.getrname(pileup.tid)

        if chrom != last_chrom and last_count > 0 and last_end:
            write_bedgraph(last_chrom, last_start, last_end + 1, last_count, normalize, out)
            last_count = 0
            last_start = None
            last_end = None

        count = 0
        if strand is None:
            if not nogaps:
                count = pileup.n
            else:
                for read in pileup.pileups:
                    op = read_cigar_at_pos(read.alignment.cigar, read.qpos, read.is_del)
                    if op != 3:
                        count += 1
        else:
            #
            #  TODO - add rev_read2 option
            #
            for read in pileup.pileups:
                if not read.alignment.is_reverse and strand == '+':
                    if not nogaps:
                        count += 1
                    else:
                        op = read_cigar_at_pos(read.alignment.cigar, read.qpos, read.is_del)
                        if op != 3:
                            count += 1

                elif read.alignment.is_reverse and strand == '-':
                    if not nogaps:
                        count += 1
                    else:
                        op = read_cigar_at_pos(read.alignment.cigar, read.qpos, read.is_del)
                        if op != 3:
                            count += 1

            # print pileup.pos,count,last_start,last_end
        if count != last_count or not last_end or (pileup.pos - last_end) > 1:
            if last_count > 0:
                write_bedgraph(last_chrom, last_start, last_end + 1, last_count, normalize, out)

            if count == 0:
                last_start = None
            # elif not last_end or (pileup.pos-last_end) > 1:
            #     last_start = pileup.pos
            else:
                last_start = pileup.pos

        last_end = pileup.pos
        last_count = count
        last_chrom = chrom

    write_bedgraph(last_chrom, last_start, last_end + 1, last_count, normalize, out)


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

    -nogaps           Don't include gaps across splice-junctions (RNA-seq)
                      Warning: adds significant processing time!
"""
    sys.exit(1)

if __name__ == "__main__":
    bam = None
    strand = None
    norm = None
    nogaps = False

    last = None
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if last == '-norm':
            norm = float(arg)
            last = None
        elif arg in ['-norm']:
            last = arg
        elif arg == '-nogaps':
            nogaps = True
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
    bam_tobedgraph(bamfile, strand, norm, nogaps, out=sys.stdout)
    bamfile.close()
