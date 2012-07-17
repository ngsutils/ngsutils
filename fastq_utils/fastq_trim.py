#!/usr/bin/env python
## category General
## desc Remove 5' and 3' linker sequences (S/W aligned)
'''
Removes linkers from 5' and 3' ends by performing local
s/w alignments with the linker sequences
'''

import os
import sys

from fastq_utils import read_fastq, is_colorspace_fastq
import support.localalign


def fastq_trim(fname, linker_5=None, linker_3=None, out=sys.stdout, pct_identity=0.8, min_trim=4, min_len=25, verbose=False):
    '''
    fname - the fastq filename
    linker_5 - the 5' linker to remove
    linker_3 - the 3' linker to remove
    out - an output stream (eg: file, stdout)
    pct_identity - the percentage of matches that must be present in the alignment to strip away linkers
    min_trim - the distance away from the edges that the linkers much match w/in
    '''

    cs = is_colorspace_fastq(fname)
    sw = support.localalign.LocalAlignment(support.localalign.NucleotideScoringMatrix(2, -1), -1)
    removed = 0
    trimmed = 0
    for name, seq, qual in read_fastq(fname):
        if verbose:
            sys.stderr.write('Read: %s\n    : %s\n' % (name, seq))
        left = 0
        right = len(seq)

        if linker_5:
            aln = sw.align(seq, linker_5)
            if verbose:
                sys.stderr.write("5' alignment:\n")
                aln.dump(sys.stderr)
            if aln.r_pos < min_trim and aln.identity > pct_identity:
                left = aln.r_end

        if linker_3:
            aln = sw.align(seq, linker_3)
            if verbose:
                sys.stderr.write("3' alignment:\n")
                aln.dump(sys.stderr)
            if aln.r_end > len(seq) - min_trim and aln.identity > pct_identity:
                right = aln.r_pos

        s = seq[left:right]
        if len(s) >= min_len:
            if left > 0 or right < len(seq):
                trimmed += 1
            if cs and len(seq) != len(qual) and left == 0:
                out.write('%s\n%s\n+\n%s\n' % (name, s, qual[left:right - 1]))
            else:
                out.write('%s\n%s\n+\n%s\n' % (name, s, qual[left:right]))
        else:
            removed += 1

        # out.write('%s\n%s (%s-%s)\n' % (name,seq,left,right))
        # out.write('x'*left)
        # out.write(seq[left:right])
        # out.write('x' *(len(seq)-right))
        # out.write('\n')
    sys.stderr.write('Trimmed: %s\n' % trimmed)
    sys.stderr.write('Removed: %s (len)\n' % removed)


def usage():
    print __doc__
    print """Usage: fastqutils trim {opts} filename.fastq{.gz}

You must specify at least one of the following:
  -5 seq      5' linker sequence to remove
  -3 seq      3' linker sequence to remove

Options
  -len val    Minimum length of a read (discards shorter) [default: 25]
  -pct val    Required percent identity (0->1.0) [default: 0.8]
  -min val    Minumum number of bases to trim (or minumum dist. from the ends)
              [default: 4]
  -v          Verbose output for each alignment
"""
    sys.exit(1)

if __name__ == '__main__':
    linker_5 = None
    linker_3 = None
    fastq = None
    last = None
    min_len = 25
    min_trim = 4
    pct_identity = 0.8
    verbose = False

    for arg in sys.argv[1:]:
        if last == '-5':
            linker_5 = arg
            last = None
        elif last == '-3':
            linker_3 = arg
            last = None
        elif last == '-pct':
            pct_identity = float(arg)
            last = None
        elif last == '-len':
            min_len = int(arg)
            last = None
        elif last == '-min':
            min_trim = int(arg)
            last = None
        elif arg == '-v':
            verbose = True
        elif arg in ['-3', '-5', '-min', '-len', '-pct']:
            last = arg
        elif not fastq:
            fastq = arg
        else:
            usage('trim')

    if not fastq or not os.path.exists(fastq) or (not linker_5 and not linker_3):
        usage()
    else:
        if linker_5:
            sys.stderr.write("Removing %s from 5' end\n" % linker_5)
        if linker_3:
            sys.stderr.write("Removing %s from 3' end\n" % linker_3)

        fastq_trim(fastq, linker_5, linker_3, min_len=min_len, pct_identity=pct_identity, min_trim=min_trim, verbose=verbose)
