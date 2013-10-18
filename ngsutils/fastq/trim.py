#!/usr/bin/env python
## category General
## desc Remove 5' and 3' linker sequences (slow, S/W aligned)
'''
Removes linkers from 5' and 3' ends by performing local S/W alignments with
the linker sequences. This is best used when the sequencing is noisy with
a lot of low quality base calls. If there aren't expected to be many
mismatches or sequencing errors, use "fastqutils filter" as it uses a faster
sliding window method.
'''

import os
import sys

from ngsutils.fastq import FASTQ
import swalign


def fastq_trim(fastq, linker_5=None, linker_3=None, out=sys.stdout, pct_identity=0.8, min_trim=4, min_len=25, verbose=False, quiet=False, failed_out=None):
    '''
    fname - the fastq filename
    linker_5 - the 5' linker to remove
    linker_3 - the 3' linker to remove
    out - an output stream (eg: file, stdout)
    pct_identity - the percentage of matches that must be present in the alignment to strip away linkers
    min_trim - the distance away from the edges that the linkers much match w/in
    failed_out - an output for failed reads
    '''

    sw = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(2, -1), -1)
    removed = 0
    trimmed = 0
    is_colorspace = fastq.is_colorspace  # preload to keep reader happy.
    for read in fastq.fetch(quiet=quiet):
        retval = seq_trim(read.name, read.seq, read.qual, linker_5, linker_3, is_colorspace, sw, pct_identity, min_trim, min_len, verbose)
        if not retval:
            if failed_out:
                read.write(failed_out)
            removed += 1
        else:
            n_seq, n_qual = retval

            if len(read.qual) != n_qual:
                trimmed += 1

            read.clone(seq=n_seq, qual=n_qual).write(out)

    if not quiet:
        sys.stderr.write('Trimmed: %s\n' % trimmed)
        sys.stderr.write('Removed: %s (len)\n' % removed)


def seq_trim(name, seq, qual, linker_5, linker_3, cs, sw, pct_identity, min_trim, min_len, verbose):
    '''
    Returns (newseq, newqual) if there is a match, otherwise: None
    '''
    if verbose:
        sys.stderr.write('\nRead: %s\n    : %s\n' % (name, seq))
    left = 0
    right = len(seq)

    if linker_5:
        aln = sw.align(seq, linker_5)
        if verbose:
            sys.stderr.write("5' alignment:\n")
            aln.dump(out=sys.stderr)
        if aln.r_pos < min_trim and aln.identity >= pct_identity:
            left = aln.r_end

    if linker_3:
        aln = sw.align(seq, linker_3)
        if verbose:
            sys.stderr.write("3' alignment:\n")
            aln.dump(out=sys.stderr)
        if aln.r_end > len(seq) - min_trim and aln.identity >= pct_identity:
            right = aln.r_pos

    s = seq[left:right]
    if len(s) >= min_len:
        if cs and len(seq) != len(qual) and left == 0:
            return (s, qual[left:right - 1])
        else:
            return (s, qual[left:right])
    else:
        return None


def usage():
    print __doc__
    print """Usage: fastqutils trim {opts} filename.fastq{.gz}

You must specify at least one of the following:
  -5 seq      5' linker sequence to remove
  -3 seq      3' linker sequence to remove

Options
  -len val         Minimum length of a read (discards shorter) [default: 25]
  -pct val         Required percent identity (0->1.0) [default: 0.8]
  -min val         Minumum number of bases to trim (or minumum dist. from the
                   ends) [default: 4]
  -failed fname    Write failed reads to file
  -v               Verbose output for each alignment
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
    failed = None
    verbose = False

    if '-test' in sys.argv[1:]:
        import doctest
        doctest.testmod()
        sys.exit()

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
        elif last == '-failed':
            if not os.path.exists(arg):
                failed = arg
                last = None
            else:
                sys.stderr.write("Failed output file: %s exists!\nNot overwriting!\n" % arg)
                sys.exit(1)
        elif arg == '-v':
            verbose = True
        elif arg in ['-3', '-5', '-min', '-len', '-pct', '-failed']:
            last = arg
        elif not fastq:
            fastq = arg
        else:
            usage()

    if not fastq or not os.path.exists(fastq) or (not linker_5 and not linker_3):
        usage()
    else:
        if linker_5:
            sys.stderr.write("Removing %s from 5' end\n" % linker_5)
        if linker_3:
            sys.stderr.write("Removing %s from 3' end\n" % linker_3)

        failed_out = None
        if failed:
            failed_out = open(failed, 'w')

        fq = FASTQ(fastq)
        fastq_trim(fq, linker_5, linker_3, min_len=min_len, pct_identity=pct_identity, min_trim=min_trim, verbose=verbose, failed_out=failed_out)
        fq.close()

        if failed_out:
            failed_out.close()
