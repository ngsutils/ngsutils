#!/usr/bin/env python
## category RNA-seq
## desc Finds regions of unusual deletions (CLIP-seq)
## experimental
"""
Finds regions of unusual deletions (CLIP-seq)

Given a set of BAM files, we search for areas where there are an unusual
amount of deletions.  For CLIP-Seq, this can be an indicator of the location
of protein-RNA interaction.

Output is either a BED file or a FASTA format file containing these hotspots.
Only unique regions are returned across all files.

See: Zhang and Darnell, Nature Biotechnology (2011)
     doi:10.1038/nbt.1873
     pmid:21633356
"""

import os
import sys
from ngsutils.bam import bam_pileup_iter
import pysam


def usage():
    print __doc__
    print """
Usage: bamutils cims {opts} in.bam {in.bam...}

Options:
    -fasta ref.fa    Ouput in FASTA format (requires reference genome.fa)
                     [default: BED output]

    -flanking N      The number of flanking bases on either side to report
                     (FASTA output only) [default: 12]

    -cutoff N        Cut-off % for deletions - if the % of reads that
                     include a deletion at a position is higher than this
                     number, the fragment is reported (0->1.0)
                     [default: 0.1]

    -ns              Don't take the strand of the read into account

    -window N        The maximum length of a deletion window
                     [default: 20]
"""
    sys.exit(1)


class BEDEmitter(object):
    def __init__(self, out=None):
        self.num = 1
        if out:
            self.out = out
        else:
            self.out = sys.stdout
        pass

    def emit(self, chrom, start, end, strand):
        if not strand:
            strand = '+'
        self.out.write('%s\t%s\t%s\tregion_%s\t%s\t%s\n' % (chrom, start, end, self.num, 0, strand))
        self.num += 1

    def close(self):
        if self.out != sys.stdout:
            self.out.close()


class FASTAEmitter(object):
    def __init__(self, ref_fname, flanking=12, out=None):
        self.num = 1
        self.ref = pysam.Fastafile(ref_fname)
        assert flanking > 0
        self.flanking = flanking

        if out:
            self.out = out
        else:
            self.out = sys.stdout

    def close(self):
        self.ref.close()
        if self.out != sys.stdout:
            self.out.close()

    def emit(self, chrom, start, end, strand):
        seq = self.ref.fetch(chrom, start - self.flanking, end + self.flanking)
        seq = '%s%s%s' % (seq[:self.flanking].upper(), seq[self.flanking:end - start + self.flanking].lower(), seq[-self.flanking:].upper())

        if strand == '-':
            rc = []
            for base in seq[::-1]:
                if base == 'A':
                    rc.append('T')
                elif base == 'T':
                    rc.append('A')
                elif base == 'C':
                    rc.append('G')
                elif base == 'G':
                    rc.append('C')
                elif base == 'a':
                    rc.append('t')
                elif base == 't':
                    rc.append('a')
                elif base == 'c':
                    rc.append('g')
                elif base == 'g':
                    rc.append('c')
            seq = ''.join(rc)

        self.out.write('>%s:%s%s%s\n%s\n' % (chrom, start - self.flanking, strand, end + self.flanking, seq))


class RegionManager(object):
    def __init__(self, emitter, strand='', max_window=20):
        self.emitter = emitter

        self.strand = strand
        self.max_window = max_window

        self.last_chrom = None
        self.start = 0
        self.end = 0
        self.del_reads = set()
        self.total_reads = set()

    def emit(self):
        if self.last_chrom:
            self.emitter.emit(self.last_chrom, self.start, self.end, self.strand)

    def reset(self, new_chrom, new_pos):
        self.last_chrom = new_chrom
        self.start = new_pos
        self.end = new_pos
        self.del_reads = set()
        self.total_reads = set()

    def add(self, chrom, pos, strand, del_reads, total_reads):
        if self.strand and strand != self.strand:
            # ignore this if the strand doesn't match
            return

        if chrom != self.last_chrom:
            self.emit()
            self.reset(chrom, pos)
        elif pos - self.start >= self.max_window:
            self.emit()
            self.reset(chrom, pos)

        self.end = pos
        self.del_reads |= del_reads
        self.total_reads |= total_reads

    def close(self):
        self.emit()


def is_read_del_at_pos(read, pos, ppos=0):
    last_op = None
    idx = 0
    for op, length in read.cigar:
        if op in [0, 1]:
            idx += length

        if pos < idx:
            if op == 2 and last_op != 3:
                return True
        last_op = op

    return False


#TODO: Check this...
def is_read_match_at_pos(read, pos):
    idx = 0
    for op, length in read.cigar:
        if op in [0, 1]:
            idx += length

        if pos < idx:
            if op == 2 or op == 0:
                return True

    return False


def bam_cims_finder(bam_fnames, output='bed', ref_fname=None, flanking=12, cutoff=0.1, stranded=True, window_size=20):
    for bam_fname in bam_fnames:
        sys.stderr.write('%s\n' % bam_fname)
        bam = pysam.Samfile(bam_fname, "rb")

        if output == 'fasta':
            emitter = FASTAEmitter(ref_fname, flanking)
        else:
            emitter = BEDEmitter()

        if stranded:
            strands = ['+', '-']
        else:
            strands = ['']

        for strand in strands:
            manager = RegionManager(emitter, strand, window_size)
            for pileup in bam_pileup_iter(bam, mask=1540):
                chrom = bam.getrname(pileup.tid)

                deletions = 0.0
                total = 0.0

                del_reads = set()
                total_reads = set()

                for pileupread in pileup.pileups:
                    if not strand or (strand == '+' and not pileupread.alignment.is_reverse) or (strand == '-' and pileupread.alignment.is_reverse):
                        if is_read_match_at_pos(pileupread.alignment, pileupread.qpos):
                            total += 1
                            total_reads.add(pileupread.alignment.qname)

                        if is_read_del_at_pos(pileupread.alignment, pileupread.qpos):
                            deletions += 1
                            del_reads.add(pileupread.alignment.qname)
                            # print ""
                            # print chrom
                            # print pileup.pos
                            # print pileupread.alignment.qname
                            # print pileupread.alignment.pos
                            # print pileupread.alignment.cigar
                            # print pileupread.qpos

                if total > 0:
                    pct = deletions / total

                    if pct > cutoff:
                        manager.add(chrom, pileup.pos, strand, del_reads, total_reads)

            manager.close()
        bam.close()
        emitter.close()

if __name__ == '__main__':
    bams = []
    ref = None
    output = 'bed'
    cutoff = 0.1
    flanking = 12
    stranded = True
    window = 20

    last = None
    for arg in sys.argv[1:]:
        if last == '-flanking':
            flanking = int(arg)
            last = None
        elif last == '-cutoff':
            cutoff = float(arg)
            last = None
        elif last == '-window':
            window = float(arg)
            last = None
        elif last == '-fasta' and not ref and os.path.exists(arg) and os.path.exists('%s.fai' % arg):
            output = 'fasta'
            ref = arg
            last = None
        elif arg == '-h':
            usage()
        elif arg in ['-flanking', '-fasta', '-cutoff', '-window']:
            last = arg
        elif arg == '-ns':
            stranded = False
        elif os.path.exists(arg) and os.path.exists('%s.bai' % arg):
            bams.append(arg)
        else:
            print "Unknown option or missing index: %s" % arg
            usage()

    if not bams:
        usage()
    else:
        bam_cims_finder(bams, output, ref, flanking, cutoff, stranded, window)
