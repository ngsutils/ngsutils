#!/usr/bin/env python
'Support classes for dealing with homopolymer files'

import sys
import struct
import io


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
    'Class for reading the homopolymer stripped index file'
    _magic = 0xCCBB601C
    record_size = struct.calcsize('<III')

    def __init__(self, fname, mode='r', verbose=False):
        self.fname = fname
        self.mode = mode
        self.verbose = verbose

        self._cur_pos = 0
        self._forest = {}  # Hash of RBTrees - one for each chrom
        self.refs = []
        self._ref_offsets = {}
        self._ref_counts = {}
        self._ref_max = {}

        if mode == 'r':
            self.fileobj = io.open(fname, 'rb')
            isize = struct.calcsize('<I')
            hsize = struct.calcsize('<H')
            qsize = struct.calcsize('<H')

            filemagic, = self.__read_bytes('<I')
            assert filemagic == HPSIndex._magic

            self.fileobj.seek(-isize, 2)
            epilog_len, = self.__read_bytes('<I')

            self.fileobj.seek(-(epilog_len + isize), 2)
            epi_count = 0
            while epi_count < epilog_len:
                reflen, = self.__read_bytes('<H')
                refname, = self.__read_bytes('<%ss' % reflen)
                count, offset, refmax = self.__read_bytes('<IQI')

                self.refs.append(refname)
                self._ref_offsets[refname] = offset
                self._ref_counts[refname] = count
                self._ref_max[refname] = refmax

                epi_count += hsize + isize + isize + qsize + reflen
                if self.verbose:
                    print '%s\t%s\t%s\t%s' % (refname, offset, count, refmax)

        elif mode == 'w':
            self.fileobj = io.open(fname, 'wb', buffering=4 * 1024 * 1024)  # use 4MB buffer
            self.fileobj.write(struct.pack('<I', HPSIndex._magic))
            self._cur_pos += struct.calcsize('<I')

        self._cur_ref = None
        self._cur_count = 0
        self._cur_genome_offset = 0
        self._cur_genome_pos = 0

    def __read_bytes(self, fmt):
        return struct.unpack(fmt, self.fileobj.read(struct.calcsize(fmt)))

    def find(self, ref, start, end=None):
        if not end:
            end = start

        record_guess = self._ref_counts[ref] * start / self._ref_max[ref]
        records = []

        low = None
        high = None

        while True:
            records.append(record_guess)
            if self.verbose:
                print 'guess: %s' % record_guess,
            if record_guess < 0:
                return
            pos, count, genome_offset = self._read_record(ref, record_guess)
            if self.verbose:
                print 'pos: %s' % pos,

            if pos == start:
                if self.verbose:
                    print 'MATCH'
                break
            else:
                if self.verbose:
                    print 'low: %s high: %s' % (low, high),

                diff = (start - pos)
                if diff > 0:
                    low = record_guess
                else:
                    high = record_guess

                # naive method just bisect low/high until you hit the right one...
                # takes into account diff b/w guess and record only to establish high/low

                # print 'low: %s high: %s \t' % (low, high),
                # if low is None or high is None:
                #     record_guess = record_guess + (self._ref_counts[ref] * diff / self._ref_max[ref])
                #     if record_guess == low:
                #         record_guess = low + 1
                #     elif record_guess == high:
                #         record_guess = high - 1
                # else:
                #     record_guess = (high + low) / 2
                # print '\tguess: %s' % record_guess

                # more complicated method - takes into account the diff b/w guess and record

                record_guess = record_guess + (self._ref_counts[ref] * diff / self._ref_max[ref])
                if self.verbose:
                    print 'guess: %s' % record_guess,
                if record_guess == low:
                    record_guess += 1
                elif record_guess == high:
                    record_guess -= 1

                if self.verbose:
                    print record_guess,
                if record_guess in records:
                    record_guess = records[-1]
                    if self.verbose:
                        print 'checked!'
                    break

        if start <= pos <= end:
            yield pos, count, genome_offset

        while pos < end:
            record_guess += 1
            pos, count, genome_offset = self._read_record(ref, record_guess)
            if start <= pos <= end:
                yield pos, count, genome_offset

    def _read_record(self, ref, record):
        offset = self._ref_offsets[ref] + (record * HPSIndex.record_size)
        self.fileobj.seek(offset, 0)
        return self.__read_bytes('<III')

    def write_ref(self, ref):
        if self.mode != 'w':
            raise ValueError

        if self._cur_ref:
            self._ref_counts[self._cur_ref] = self._cur_count
            self._ref_max[self._cur_ref] = self._cur_genome_pos

        self.refs.append(ref)
        self._cur_ref = ref
        self._cur_count = 0
        self._ref_offsets[ref] = self._cur_pos

    def write(self, pos, count):
        if self.mode != 'w':
            raise ValueError

        if count > 0xFFFFFFFF:
            raise ValueError("Repeat-count is too high at position: %s (%s)" % (pos, count))
        # elif count > 0x7FFF:
        #     low = (count & 0x7FFF) | 0x8000  # low 15 bits, plus flag on bit 16
        #     high = count >> 15
        #     self._cur_pos += struct.calcsize('<IHH')
        #     self.fileobj.write(struct.pack('<IHH', pos, low, high))
        else:
            self.fileobj.write(struct.pack('<III', pos, count, self._cur_genome_offset))
            self._cur_pos += struct.calcsize('<III')

        self._cur_count += 1
        self._cur_genome_offset += count
        self._cur_genome_pos = pos

    def close(self):
        if self.mode == 'w':
            if self._cur_ref:
                self._ref_counts[self._cur_ref] = self._cur_count
                self._ref_max[self._cur_ref] = self._cur_genome_pos
            footer_count = 0
            for ref in self.refs:
                count = 0
                offset = 0
                refmax = 0
                if ref in self._ref_counts:
                    count = self._ref_counts[ref]
                if ref in self._ref_offsets:
                    offset = self._ref_offsets[ref]
                if ref in self._ref_max:
                    refmax = self._ref_max[ref]

                self.fileobj.write(struct.pack('<H%ssIQI' % len(ref), len(ref), ref, count, offset, refmax))
                footer_count += struct.calcsize('<H%ssIQI' % len(ref))

            self.fileobj.write(struct.pack('<I', footer_count))
        self.fileobj.close()

if __name__ == '__main__':
    idx = HPSIndex(sys.argv[1], verbose=True)
    if len(sys.argv) > 2:
        for arg in sys.argv[2:]:
            print arg
            chrom, startend = arg.split(':')
            if '-' in startend:
                start, end = [int(x) for x in startend.split('-')]
            else:
                start = int(startend)
                end = None

            for tup in idx.find(chrom, start, end):
                print '**** Found:',
                print tup
