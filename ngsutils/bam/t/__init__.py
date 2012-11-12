from mock import *
import ngsutils.bam


def _matches(valid, queries):
    extra = False
    check = [False, ] * len(valid)
    for q in queries:
        found = False
        for i, v in enumerate(valid):
            if v == q:
                check[i] = True
                found = True
                break
        if not found:
            extra = True

    if not check or extra:
        return False

    return True


class MockBam(object):
    def __init__(self, refs):
        self._refs = refs[:]
        self._reads = {}
        self._read_keys = []
        self.filename = ''
        self.__iter = None
        self.__pos = 0

    def tell(self):
        return self.__pos

    def add_read(self, read, *args, **kwargs):
        if type(read) != MockRead:
            read = MockRead(read, *args, **kwargs)

        k = (read.tid, read.pos, len(self._reads))
        self._reads[k] = read
        self._read_keys.append(k)
        self._read_keys.sort()

    def __iter__(self):
        self.__iter = self.fetch()
        return self

    def next(self):
        try:
            return self.__iter.next()
        except:
            self.__iter = None
            raise StopIteration

    def fetch(self, ref=None, start=-1, end=-1, region=None):
        if region:
            print "Unsupported option!"
            sys.exit(1)
        self.__pos = 0
        tid = -1

        if ref and start > -1 and end > -1:
            for i, val in enumerate(self._refs):
                if val == ref:
                    tid = i
                    break
            if tid == -1:
                raise ValueError("Missing reference: %s" % ref)

        for k in self._read_keys:
            self.__pos += 1
            read = self._reads[k]
            if tid == -1:  # no limits, return them all
                yield read
            else:
                if self._reads[k].tid == tid:
                    if start <= read.pos <= end or start <= read.aend <= end or (read.pos < start and read.aend > end):
                        yield read

    @property
    def references(self):
        return self._refs[:]

    def getrname(self, tid):
        return self._refs[tid]

    def close(self):
        pass


class MockRead(object):
    def __init__(self, qname, seq=None, qual=None, tid=-1, pos=-1, aend=None, cigar=None, tags=None, mapq=0, rnext=-1, pnext=-1, isize=0, tlen=0, is_reverse=False, is_paired=False, is_secondary=False, is_qcfail=False, mate_is_unmapped=False, mate_is_reverse=False, is_read1=False, is_read2=False, flag=0):
        self.qname = qname
        self.seq = seq
        self.qual = qual
        self.tid = tid
        self.pos = pos
        self.aend = aend
        self.mapq = mapq
        self.rnext = rnext
        self.pnext = pnext
        self.isize = isize
        self.tlen = tlen

        self.is_paired = is_paired
        self.is_reverse = is_reverse
        self.is_secondary = is_secondary
        self.is_qcfail = is_qcfail
        self.mate_is_reverse = mate_is_reverse
        self.is_read1 = is_read1
        self.is_read2 = is_read2

        if flag:
            self._flag = flag
            self.is_paired = 0x1 & self.flag > 0
            self.is_unmapped = 0x4 & self.flag > 0
            self.mate_is_unmapped = 0x8 & self.flag > 0
            self.is_reverse = 0x10 & self.flag > 0
            self.mate_is_reverse = 0x20 & self.flag > 0
            self.is_read1 = 0x40 & self.flag > 0
            self.is_read2 = 0x80 & self.flag > 0
            self.is_secondary = 0x100 & self.flag > 0
            self.is_qcfail = 0x200 & self.flag > 0
        else:
            self._flag = 0
            self._flag |= 0x1 if self.is_paired else 0
            self._flag |= 0x4 if self.is_unmapped else 0
            self._flag |= 0x8 if self.mate_is_unmapped else 0
            self._flag |= 0x10 if self.is_reverse else 0
            self._flag |= 0x20 if self.mate_is_reverse else 0
            self._flag |= 0x40 if self.is_read1 else 0
            self._flag |= 0x80 if self.is_read2 else 0
            self._flag |= 0x100 if self.is_secondary else 0
            self._flag |= 0x200 if self.is_qcfail else 0

        if self.tid == -1 or self.pos == -1:
            self.is_unmapped = True

        if self.is_paired:
            if self.rnext == -1 or self.pnext == -1:
                self.mate_is_unmapped = True

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
    def flag(self):
        return self._flag

    @property
    def is_unmapped(self):
        if self.tid > -1 and self.pos > -1:
            return False
        return True

    @is_unmapped.setter
    def is_unmapped(self, val):
        if val:
            self.tid = -1
            self.pos = -1
            self._flag |= 0x4
        else:
            self.tid = 0
            self.pos = 0
            self._flag = self._flag & 0xFFFB  # removes 0x4

    @property
    def mate_is_unmapped(self):
        if not self.is_paired:
            return False

        if self.rnext > -1 and self.pnext > -1:
            return False
        return True

    @mate_is_unmapped.setter
    def mate_is_unmapped(self, val):
        if val:
            self.rnext = -1
            self.pnext = -1
            self._flag |= 0x8
        else:
            self.rnext = 0
            self.pnext = 0
            self._flag = self._flag & 0xFFF7  # removes 0x8
