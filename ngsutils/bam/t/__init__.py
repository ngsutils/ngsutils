from mock import *
import ngsutils.bam


class MockBam(object):
    def __init__(self, refs):
        self._refs = {}
        for i, val in enumerate(refs):
            self._refs[i + 1] = val

    def getrname(self, tid):
        return self._refs[tid]


class MockRead(object):
    def __init__(self, qname, seq=None, qual=None, tid=-1, pos=-1, aend=None, is_reverse=False, cigar=None, tags=None):
        self.qname = qname
        self.seq = seq
        self.qual = qual
        self.tid = tid
        self.pos = pos
        self.aend = aend
        self.is_reverse = is_reverse

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
    def is_unmapped(self):
        if self.tid and self.pos > -1:
            return False
        return True

    @is_unmapped.setter
    def is_unmapped(self, val):
        if val:
            self.tid = -1
            self.pos = -1
        else:
            self.tid = 1
            self.pos = 0
