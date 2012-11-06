from mock import *
import ngsutils.bam


class MockBam(object):
    def __init__(self, refs):
        self._refs = {}
        for i, val in enumerate(refs):
            self._refs[i] = val
        self._reads = {}
        self._read_keys = []

    def fetch(self, ref=None, start=-1, end=-1, region=None):
        if region:
            print "Unsupported option!"
            sys.exit(1)

        tid = -1

        if ref and start > -1 and end > -1:
            for i, val in enumerate(self._refs):
                if val == ref:
                    tid = i
                    break
            if tid == -1:
                raise ValueError("Missing reference: %s" % ref)

        for k in self._read_keys:
            read = self._reads[k]
            if tid == -1:
                yield read
            else:
                if self._reads[k].tid == tid:
                    if start <= read.pos <= end or start <= read.aend <= end or (read.pos < start and reads.aend > end):
                        yield read

    def getrname(self, tid):
        return self._refs[tid]

    def close(self):
        pass


class MockRead(object):
    def __init__(self, qname, seq=None, qual=None, tid=-1, pos=-1, aend=None, cigar=None, tags=None, is_reverse=False, is_paired=False, is_mate_unmapped=False, is_secondary=False, is_qcfail=False, flags=0):
        self.qname = qname
        self.seq = seq
        self.qual = qual
        self.tid = tid
        self.pos = pos
        self.aend = aend

        self.is_paired = is_paired
        self.is_mate_unmapped = is_mate_unmapped
        self.is_reverse = is_reverse
        self.is_secondary = is_secondary
        self.is_qcfail = is_qcfail

        if flags:
            self._flags = flags
            self.is_paired = 0x1 & self.flags > 0
            self.is_unmapped = 0x4 & self.flags > 0
            self.is_mate_unmapped = 0x8 & self.flags > 0
            self.is_reverse = 0x10 & self.flags > 0
            self.is_secondary = 0x100 & self.flags > 0
            self.is_qcfail = 0x200 & self.flags > 0
        else:
            self._flags = 0
            self._flags |= 0x1 if self.is_paired else 0
            self._flags |= 0x4 if self.is_unmapped else 0
            self._flags |= 0x8 if self.is_mate_unmapped else 0
            self._flags |= 0x10 if self.is_reverse else 0
            self._flags |= 0x100 if self.is_secondary else 0
            self._flags |= 0x200 if self.is_qcfail else 0

        if self.tid == -1 or self.pos == -1:
            self.is_unmapped = True

        if cigar:
            self.cigar = ngsutils.bam.cigar_fromstr(cigar)
        else:
            self.cigar = None

        if tags:
            self.tags = tags
        else:
            self.tags = []

    def opt(self, k):
        for k1, val in self.tags:
            if k1 == k:
                return val
        return None

    @property
    def flags(self):
        return self._flags

    @property
    def is_unmapped(self):
        if self.tid and self.pos > -1:
            return False
        return True

    @is_unmapped.setter
    def is_unmapped(self, val):
        if val:
            self.tid = -1
            self.pos = -1
            self._flags |= 0x4
        else:
            self.tid = 1
            self.pos = 0
            self._flags = self._flags & 0xFFFB  # removes 0x4
