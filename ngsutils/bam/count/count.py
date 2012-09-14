import ngsutils.support.stats
import sys
import pysam


class Model(object):
    def get_source(self):
        raise NotImplemented

    def get_name(self):
        raise NotImplemented

    def get_headers(self):
        raise NotImplemented

    def get_regions(self):
        '''
        Generator that yields a tuple:
        chrom,
        starts (list),
        ends (list),
        strand,
        columns to include (list-column values that will be output with counts),
        callback (function): args (bam, count, reads, outcols)
                             This should return (or yield) a list of cols to be written
                             to stdout.

        The start/end values are lists so that they can include multiple exons
        at once. For most usecases this will only have one value, but it makes
        checking exons in a gene easier.
        '''
        raise NotImplemented

    def get_postheaders(self):
        return None

    def count(self, bamfile, stranded, coverage, uniq_only, rpkm, norm, multiple, whitelist, blacklist):
        bam = pysam.Samfile(bamfile, 'rb')

        region_counts = []
        multireads = set()
        single_count = 0

        for chrom, starts, ends, strand, cols, callback in self.get_regions():
            outcols = cols[:]

            coding_len = 0
            for s, e in zip(starts, ends):
                coding_len += e - s
            outcols.append(coding_len)

            count, reads = _fetch_reads(bam, chrom, strand if stranded else None, starts, ends, multiple, False, whitelist, blacklist, uniq_only)
            outcols.append(None)

            for read in reads:
                if read.tags and 'IH' in read.tags:
                    ih = int(read.opt('IH'))
                else:
                    ih = 1

                if ih == 1:
                    single_count += 1
                else:
                    multireads.add(read.qname)

            if coverage:
                mean, stdev, median = calc_coverage(bam, chrom, strand, starts, ends, whitelist, blacklist)
                outcols.append(mean, stdev, median)

            if callback:
                for callback_cols in callback(bam, count, reads, outcols):
                    region_counts.append((count, coding_len, callback_cols))
            else:
                region_counts.append((count, coding_len, outcols))

        sys.stderr.write('Calculating normalization...')

        norm_val = None
        norm_val_orig = None

        if norm == 'all':
            norm_val_orig = _find_mapped_count(bam, whitelist, blacklist)
        elif norm == 'mapped':
            norm_val_orig = single_count + len(multireads)
        elif norm == 'quantile':
            norm_val_orig = _find_mapped_count_pcts([x[0] for x in region_counts])
        elif norm == 'median':
            norm_val_orig = _find_mapped_count_median([x[0] for x in region_counts])

        if norm_val_orig:
            norm_val = float(norm_val_orig) / 1000000

        sys.stderr.write('\n')

        sys.stdout.write('## input %s\n' % bamfile)
        sys.stdout.write('## model %s %s\n' % (self.get_name(), self.get_source()))
        sys.stdout.write('## stranded %s\n' % stranded)
        sys.stdout.write('## multiple %s\n' % multiple)
        if norm_val:
            sys.stdout.write('## norm %s %s\n' % (norm, norm_val_orig))
            sys.stdout.write('## CPM-factor %s\n' % norm_val)

        sys.stdout.write('#')
        sys.stdout.write('\t'.join(self.get_headers()))
        sys.stdout.write('\tlength\tcount')
        if norm_val:
            sys.stdout.write('\tcount (CPM)')
            if rpkm:
                sys.stdout.write('\tRPKM')

        if coverage:
            sys.stdout.write('\tcoverage mean\tcoverage stdev\tcoverage median')

        if self.get_postheaders():
            sys.stdout.write('\t')
            sys.stdout.write('\t'.join(self.get_postheaders()))

        sys.stdout.write('\n')

        for count, coding_len, outcols in region_counts:
            first = True
            for col in outcols:
                if not first:
                    sys.stdout.write('\t')
                first = False

                if col is None:  # this is the marker for the 'count' col
                    sys.stdout.write('%s\t' % count)

                    if norm_val:
                        sys.stdout.write(str(count / norm_val))
                        if rpkm:
                            sys.stdout.write('\t')
                            sys.stdout.write(str(count / (coding_len / 1000.0) / norm_val))

                else:
                    sys.stdout.write(str(col))

            sys.stdout.write('\n')


def _calc_read_regions(read):
    regions = []
    start = read.pos
    end = read.pos
    for op, length in read.cigar:
        if op == 0:
            end += length
        elif op == 1:
            pass
        elif op == 2:
            end += length
        elif op == 3:
            regions.append((start, end))
            end += length
            start = end

    regions.append((start, end))

    return regions


def _fetch_reads_excluding(bam, chrom, strand, start, end, multiple, whitelist=None, blacklist=None):
    '''
    Find reads that exclude this region.

    Example:  This read excludes region 'B'

    ------AAAAAAAAA---------BBBBBBBBB---------------------CCCCCCCC---------
              +++++.......................................+++

    '''

    reads = set()
    count = 0

    if not chrom in bam.references:
        return count, reads

    for read in bam.fetch(chrom, start, end):
        if not strand or (strand == '+' and not read.is_reverse) or (strand == '-' and read.is_reverse):
            excl = True
            for s, e in _calc_read_regions(read):
                if start <= s <= end or start <= e <= end:
                    excl = False
                    break
            if excl and not read in reads:
                reads.add(read)
                count += 1
    return count, reads


def _fetch_reads(bam, chrom, strand, starts, ends, multiple, exclusive, whitelist=None, blacklist=None, uniq=False):
    '''
    Find reads that match within the given regions...

    starts and ends are lists of start and end positions for the various regions
    we are interested in.

    Multiple dictates how reads that map to multiple locations are handled.
    Possible values:
        'complete' - a full 1.0 is added to the count, regardless of how many other places
                     this read maps.
        'partial'  - a fraction of a point is added to the count, based upon how many
                     other places this read maps.  For example, if a read mapped to 3 places
                     on the genome, it would add only 1/3 to the read count.
        'ignore'   - reads that map to multiple places on the genome are ignored.

    Exclusive is a boolean.  If true, then reads must start AND end w/in the given regions in order to
    count.  If false, the read must just be present in the region (start or end).  This is helpful
    in finding reads that map to constant regions as opposed to alternatively expressed regions.

    White/blacklist are used to restrict the counting to:
        include only reads in the whitelist
        or exclude those in the blacklist
    (mainly useful for only benchmarking)

    uniq will only include reads with uniq starting positions (strand specific)

    '''
    assert multiple in ['complete', 'partial', 'ignore']

    reads = set()
    start_pos = set()
    count = 0

    if not chrom in bam.references:
        return count, reads

    for s, e in zip(starts, ends):
        for read in bam.fetch(chrom, s, e):
            if blacklist and read.qname in blacklist:
                continue
            if not whitelist or read.qname in whitelist:
                if read.is_reverse:
                    k = (read.aend, '-')
                else:
                    k = (read.pos, '+')

                if uniq and k in start_pos:
                        continue

                if not strand or (strand == '+' and not read.is_reverse) or (strand == '-' and read.is_reverse):
                    if exclusive:
                        start_ok = False
                        end_ok = False

                        for s1, e1 in zip(starts, ends):
                            if s1 <= read.pos <= e1:
                                start_ok = True
                                break

                        if not start_ok:
                            continue

                        for s1, e1 in zip(starts, ends):
                            if s1 <= read.aend <= e1:
                                end_ok = True
                                break

                        if not end_ok:
                            continue

                    if not read in reads:
                        start_pos.add(k)
                        reads.add(read)

                        if read.tags and 'IH' in read.tags:
                            ih = int(read.opt('IH'))
                        else:
                            ih = 1

                        if ih == 1 or multiple == 'complete':
                            count += 1
                        elif multiple == 'partial':
                            count += (1.0 / ih)
                        else:  # multiple = ignore
                            pass
    return count, reads


def calc_coverage(bam, chrom, strand, starts, ends, whitelist, blacklist):
    coverage = []

    for start, end in zip(starts, ends):
        for pileup in bam.pileup(chrom, start, end):
            count = 0
            for pileupread in pileup.pileups:
                if blacklist and pileupread.qname in blacklist:
                    continue
                if not whitelist or pileupread.qname in whitelist:
                    if strand:
                        if strand == '+' and pileupread.alignment.is_reverse:
                            continue
                        elif strand == '-' and not pileupread.alignment.is_reverse:
                            continue

                    if not pileupread.is_del:
                        count += 1

            coverage.append(count)

    mean, stdev = ngsutils.support.stats.mean_stdev(coverage)
    coverage.sort()

    if len(coverage) % 2 == 1:
        median = coverage[len(coverage) / 2]
    else:
        a = coverage[(len(coverage) / 2) - 1]
        b = coverage[(len(coverage) / 2)]
        median = (a + b) / 2
    return mean, stdev, median


def _find_mapped_count_median(counts):
    counts.sort()
    filtered = [x for x in counts if x > 0]

    return filtered[len(filtered) / 2]


def _find_mapped_count_pcts(counts, min_pct=0.0, max_pct=0.75):
    counts.sort()
    filtered = [x for x in counts if x > 0]
    acc = 0
    min_i = int(min_pct * len(filtered))
    max_i = int(max_pct * len(filtered))

    for count in filtered[min_i:max_i]:
        acc += count

    return acc


def _find_mapped_count(bam, whitelist=None, blacklist=None):
    sys.stderr.write('Finding number of mapped reads\n')
    bam.seek(0)
    mapped_count = 0
    multinames = set()
    for read in bam.fetch():
        if blacklist and read.qname in blacklist:
            continue
        if not whitelist or read.qname in whitelist:
            if not read.is_unmapped and not read.qname in multinames:
                mapped_count += 1
                if read.tags and 'IH' in read.tags:
                    ih = int(read.opt('IH'))
                else:
                    ih = 1
                if ih > 1:
                    multinames.add(read.qname)
    bam.seek(0)
    sys.stderr.write("%s mapped reads\n" % mapped_count)
    return mapped_count
