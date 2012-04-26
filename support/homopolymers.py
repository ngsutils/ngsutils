#!/usr/bin/env python
'Support classes for dealing with homopolymer files'

import sys
import struct


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

    def __init__(self, fname, mode='r'):
        self.fname = fname
        self.mode = mode

        self._cur_pos = 0
        self._forest = {}  # Hash of RBTrees - one for each chrom
        self.refs = []
        self._ref_offsets = {}
        self._ref_counts = {}

        if mode == 'r':
            self.fileobj = open(fname)
            isize = struct.calcsize('<I')
            hsize = struct.calcsize('<H')

            filemagic, = self.__read_bytes('<I')
            assert filemagic == HPSIndex._magic

            self.fileobj.seek(isize, 2)
            epilog_len, = self.__read_bytes('<I')

            self.fileobj.seek(epilog_len + isize, 2)
            epi_count = 0
            while epi_count < epilog_len:
                reflen, = self.__read_bytes('<H')
                refname, = self.__read_bytes('<s', reflen)
                count, offset = self.__read_bytes('<II')

                self.refs.append(refname)
                self._ref_offsets[refname] = offset
                self._ref_counts[refname] = count

                print 'ref: %s count: %s offset: %s' % (refname, count, offset)

                epi_count += hsize + isize + isize + reflen

        elif mode == 'w':
            self.fileobj = open(fname, 'w')
            self.fileobj.write(struct.pack('<I', HPSIndex._magic))
            self._cur_pos += struct.calcsize('<I')

        self._cur_ref = None
        self._cur_count = 0

    def __read_bytes(self, fmt, size=0):
        if size:
            return struct.unpack(fmt, self.fileobj.read(size))
        return struct.unpack(fmt, self.fileobj.read(struct.calcsize(fmt)))

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

        if count < 32768:
            self.fileobj.write(struct.pack('<IH', pos, count | 0x8000))
            self._cur_pos += struct.calcsize('<IH')
        elif count < 2147483648:
            self.fileobj.write(struct.pack('<II', pos, count))
            self._cur_pos += struct.calcsize('<II')
        else:
            raise ValueError("Repeat-count is too high at position: %s (%s)"
                             % (pos, count))

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

if __name__ == '__main__':
    idx = HPSIndex(sys.argv[1])
