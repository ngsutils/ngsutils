import ngsutils.support.ngs_utils


class BedFile(object):
    '''
    For non-Tabix indexed BED files, read in the entire file, allowing for
    iteration from region to region.

    This reads the entire file into memory, in a series of bins. Each bin
    is ~100K in size. Each bin can then be iterated over.

    This is less efficient than using a proper Tree, but in reality, this
    usually isn't an issue.
    '''

    _bin_const = 100000

    def __init__(self, fname=None, fileobj=None):
        self._bins = {}
        self._bin_list = []
        self._cur_bin_idx = 0
        self._cur_bin_pos = 0

        if not fname and not fileobj:
            raise ValueError("Must specify either filename or fileobj")

        if fname:
            with ngsutils.support.ngs_utils.gzip_opener(fname) as fobj:
                self._read(fobj)
        else:
            self._read(fileobj)

    def _read(self, fobj):
        for line in fobj:
            line = line.strip()
            if line and line[0] != '#':
                cols = line.split('\t')
                while len(cols) < 6:
                    cols.append('')

                region = BedRegion(*cols)
                startbin = region.start / BedFile._bin_const
                endbin = region.end / BedFile._bin_const

                for bin in xrange(startbin, endbin + 1):
                    if not (region.chrom, bin) in self._bins:
                        self._bin_list.append((region.chrom, bin))
                        self._bins[(region.chrom, bin)] = []
                self._bins[(region.chrom, bin)].append(region)
        self._bin_list.sort()

    def fetch(self, chrom, start, end, strand=None):
        ''' For non-TABIX indexed BED files, find all regions w/in a range '''
        startbin = start / BedFile._bin_const
        endbin = end / BedFile._bin_const

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
        self._cur_bin_idx = 0
        self._cur_bin_pos = 0
        return self

    def tell(self):
        return self._cur_bin_idx

    def close(self):
        pass

    def next(self):
        if self._cur_bin_idx >= len(self._bin_list):
            raise StopIteration

        binvals = self._bins[self._bin_list[self._cur_bin_idx]]
        if self._cur_bin_pos < len(binvals):
            val = binvals[self._cur_bin_pos]
            self._cur_bin_pos += 1
            return val
        else:
            self._cur_bin_idx += 1
            self._cur_bin_pos = 0
            return self.next()


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
        if self.score and self.name and self.strand:
            return '\t'. join([str(x) for x in [self.chrom, self.start, self.end, self.name, self.score, self.strand]])
        return '\t'. join([str(x) for x in [self.chrom, self.start, self.end]])
