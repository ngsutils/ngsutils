#!/usr/bin/env python
'''
Removes runs of homopolymers in a FASTA file

This will remove all runs of repeated nucleotides in a FASTA file. This can be
used to help map Ion Torrent or 454 reads, which are succeptible to errant
homopolymer calls.

This produces two files:
    1) a new FASTA file with the homopolymers removed
    2) an index file that allows for converting the new coordinates to the
       original genomic positions

The index is currently a gzipped text file in the format: (subject to change)
>reference
pos\trepeat_count\ttotal_offset
pos\trepeat_count\ttotal_offset
'''

import os
import sys
import struct
import gzip

from support.eta import ETA


class FASTAWriter(object):
    def __init__(self, fileobj=sys.stdout, wrap=50):
        self.fileobj = fileobj
        self.wrap = wrap

        self._line_count = 0
        self._first = True

    def write_ref(self, ref):
        if not self._first:
            self.fileobj.write('\n')

        self.fileobj.write('>%s\n' % ref)
        self._first = False
        self._line_count = 0

    def write(self, seq):
        for s in seq:
            self.fileobj.write(s)
            self._line_count += 1
            if self._line_count >= self.wrap:
                self.fileobj.write('\n')
                self._line_count = 0

    def close(self):
        self.write('\n')
        if self.fileobj != sys.stdout:
            self.fileobj.close()


class HPSIndex(object):
    'Class for reading the homopolymer stripped index file - format subject to change'
    def __init__(self, fname, mode='r'):
        self.fname = fname
        self.mode = mode

        self._cur_pos = 0

        if mode == 'r':
            self.fileobj = open(fname)
        elif mode == 'w':
            self.fileobj = open(fname, 'w')
            self.fileobj.write(struct.pack('<I', 0xCCBB601C))

        self.refs = []
        self._ref_offsets = {}
        self._ref_counts = {}

        self._cur_ref = None
        self._cur_count = 0

    def write_ref(self, ref):
        if self.mode != 'w':
            raise ValueError

        if self._cur_ref:
            self._ref_counts[self._cur_ref] = self._cur_count

        self.refs.append(ref)
        self._cur_ref = ref
        self._cur_count = 0
        self._ref_offsets[ref] = self._cur_pos

    def write(self, pos, count):
        if self.mode != 'w':
            raise ValueError

        self.fileobj.write(struct.pack('<II', pos, count))
        self._cur_count += 1

    def close(self):
        if self.mode == 'w':
            s = ''
            for ref in self.refs:
                count = 0
                offset = 0
                if ref in self._ref_counts:
                    count = self._ref_counts[ref]
                if ref in self._ref_offsets:
                    offset = self._ref_offsets[ref]
                s += struct.pack('<HsII', len(ref), ref, count, offset)
            self.fileobj.write(s)
            self.fileobj.write(struct.pack('<I', len(s)))
        self.fileobj.close()


def read_fasta_bases(fname):
    if fname[-3:] == '.gz':
        fobj = gzip.open(fname)
    else:
        fobj = open(fname)

    eta = ETA(os.stat(fname).st_size, fileobj=fobj)

    ref = None
    for line in fobj:
        eta.print_status(extra=ref)
        line = line.strip()
        if line[0] == '>':
            ref = line[1:]
        else:
            for base in line:
                yield (ref, base.upper())

    fobj.close()
    eta.done()


def seq_strip_homopolymer(fname, outfa_name=None, outidx_name=None, suffix=None):
    if outfa_name:
        outwriter = FASTAWriter(open(outfa_name, 'w'))

    if outidx_name:
        outidx = HPSIndex(outidx_name, 'w')

    lastref = None
    lastbase = None

    strip_pos = 0
    repeat_count = 0

    for ref, base in read_fasta_bases(fname):
        if ref != lastref:
            if outfa_name:
                if suffix:
                    outwriter.write_ref('%s%s' % (ref, suffix))
                else:
                    outwriter.write_ref(ref)

            if outidx_name:
                outidx.write_ref(ref)

            lastref = ref
            lastbase = None

            strip_pos = 0
            repeat_count = 0

        if base == lastbase:
            repeat_count += 1
        else:
            lastbase = base
            if repeat_count > 0 and outidx_name:
                outidx.write(strip_pos, repeat_count)
            strip_pos += 1
            repeat_count = 0
            if outfa_name:
                outwriter.write(base)

    # if this ends in a repeat...
    if repeat_count > 0 and outidx_name:
        outidx.write(strip_pos, repeat_count)

    if outfa_name:
        outwriter.close()
    if outidx_name:
        outidx.close()


def usage():
    print __doc__
    print """Usage: sequtils strip_homopolymer {opts} infile.fa

Options:
    ** You must select at least one of these **
    -fa  fname   Output a FASTA file with the homopolymers removed
    -idx fname   Output an index file with the location of the homopolymers

    -suf val     Suffix for reference names (include space for comment)
"""
    sys.exit(1)


if __name__ == '__main__':
    fname = None
    outfa = None
    outidx = None
    suffix = None
    last = None
    for arg in sys.argv[1:]:
        if last == '-suf':
            suffix = arg
            last = None
        elif last == '-fa':
            outfa = arg
            last = None
        elif last == '-idx':
            outidx = arg
            last = None
        elif arg in ['-suf', '-fa', '-idx']:
            last = arg
        elif not fname:
            if not os.path.exists(arg):
                print "Missing file: %s" % arg
                usage()
            fname = arg

    if not fname or (not outfa and not outidx):
        usage()

    seq_strip_homopolymer(fname, outfa, outidx, suffix)
