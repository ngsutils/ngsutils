import gzip


class BedFile(object):
    def __init__(self, fname):
        self.fname = fname
        self._fobj = None
        self._bins = {}

    def fetch(self, chrom, start, end, strand=None):
        ''' For non-TABIX indexed BED files, find all regions w/in a range '''
        if not self._bins:
            for region in self:
                startbin = region.start / 100000
                endbin = region.end / 100000

                for bin in xrange(startbin, endbin + 1):
                    if not (region.chrom, bin) in self._bins:
                        self._bins[(chrom, bin)] = []
                self._bins[(chrom, bin)].append(region)

        startbin = start / 100000
        endbin = end / 100000

        for bin in xrange(startbin, endbin + 1):
            if (chrom, bin) in self._bins:
                for region in self._bins[(chrom, bin)]:
                    if strand and strand != region.strand:
                        continue
                    if start <= region.start <= end or start <= region.end <= end:
                        yield region
                    elif region.start < start and region.end > end:
                        yield region

    def __iter__(self):
        if self.fname[-4:] == '.bgz' or self.fname[-3:] == '.gz':
            self._fobj = gzip.open(self.fname)
        else:
            self._fobj = open(self.fname)
        return self

    def tell(self):
        if self._fobj:
            return self._fobj.tell()
        else:
            return -1

    def close(self):
        if self._fobj:
            self._fobj.close()
        self._bins = None

    def next(self):
        valid = False
        while not valid:
            line = self._fobj.next()
            if not line:
                self._fobj.close()
                self._fobj = None
                raise StopIteration
            line = line.strip()
            if line and line[0] != '#':
                cols = line.split('\t')
                while len(cols) < 6:
                    cols.append('')

                return BedRegion(*cols)


class BedRegion(object):
    def __init__(self, chrom, start, end, name, score, strand):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        if score == '':
            self.score = 0
        else:
            self.score = int(score)
        if strand == '':
            self.strand = None
        else:
            self.strand = strand

    def __repr__(self):
        return '%s|%s|%s|%s|%s|%s' % (self.chrom, self.start, self.end, self.name, self.score, self.strand)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
