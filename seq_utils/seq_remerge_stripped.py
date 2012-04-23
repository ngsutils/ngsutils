#!/usr/bin/env python
'''
Takes a FASTA file stripped of homopolymers and an accompanying index file
and recreates the original FASTA file.
'''

import os
import sys

from support.eta import ETA
from seq_strip_homopolymers import HPSIndex, FASTAWriter


def seq_remerge_stripped(faname, idxname, quiet=False):
    fawriter = FASTAWriter()
    fa = open(faname)
    idx = HPSIndex(idxname)

    if not quiet:
        eta = ETA(os.stat(faname).st_size, fileobj=fa)

    ref = None
    lastbase = None

    idx_ref, idx_pos, idx_count, idx_offset = idx.next()
    for line in fa:
        if not quiet:
            eta.print_status(extra=ref)
        line = line.strip()
        if line[0] == '>':
            ref = line[1:].split(' ')[0]  # remove comments...
            strip_pos = 0
            fawriter.write_ref(ref)
        else:
            for base in line:
                if idx_ref == ref and idx_pos == strip_pos:
                    fawriter.write(lastbase * idx_count)
                    try:
                        idx_ref, idx_pos, idx_count, idx_offset = idx.next()
                    except StopIteration:
                        idx_ref = None

                fawriter.write(base)
                lastbase = base
                strip_pos += 1

    # if this ends in a repeat...
    if idx_ref == ref and idx_pos == strip_pos:
        fawriter.write(lastbase * idx_count)

    sys.stdout.write('\n')

    fa.close()
    idx.close()
    eta.done()


def usage():
    print __doc__
    print "Usage: sequtils remerge_stripped stripped.fa stripped.idx"
    sys.exit(1)


if __name__ == '__main__':
    faname = None
    idxname = None
    quiet = False

    for arg in sys.argv[1:]:
        if arg == '-q':
            quiet = True
        elif not faname:
            faname = arg
        elif not idxname:
            idxname = arg

    if not faname or not idxname:
        usage()

    seq_remerge_stripped(faname, idxname, quiet=quiet)
