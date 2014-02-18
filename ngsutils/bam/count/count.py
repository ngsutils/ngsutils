import ngsutils.support.stats
import sys
import tempfile
import ngsutils

from ngsutils.bam.t import MockBam
assert(MockBam)  # just for linting... it is used in a doctest


class TmpCountFile(object):
    def __init__(self):
        self.tmpfile = tempfile.TemporaryFile()

    def write(self, count, coding_len, cols):
        self.tmpfile.write('%s\t%s\t%s\n' % (count, coding_len, '\t'.join([str(x) for x in cols])))

    def fetch(self):
        self.tmpfile.flush()
        self.tmpfile.seek(0)
        for line in self.tmpfile:
            cols = line.strip('\n').split('\t')
            yield (int(cols[0]), int(cols[1]), cols[2:])

    def close(self):
        self.tmpfile.close()


class Model(object):
    def __init__(self):
        pass

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

    def count(self, bam, library_type='FR', coverage=False, uniq_only=False, fpkm=False, norm='', multiple='complete', whitelist=None, blacklist=None, out=sys.stdout, quiet=False, start_only=False):
        # bam = pysam.Samfile(bamfile, 'rb')

        # region_counts = []
        # multireads = set()
        # single_count = 0
        tmpcounts = TmpCountFile()

        counts_tally = {}
        total_count = 0.0

        if library_type in ['FR', 'RF']:
            stranded = True
        else:
            stranded = False

        for chrom, starts, ends, strand, cols, callback in self.get_regions():
            outcols = cols[:]

            coding_len = 0
            for s, e in zip(starts, ends):
                coding_len += e - s
            outcols.append(coding_len)

            count, reads = _fetch_reads(bam, chrom, strand if stranded else None, starts, ends, multiple, False, whitelist, blacklist, uniq_only, library_type, start_only)
            outcols.append('')
            total_count += count

            # for read in reads:
            #     if read.tags and 'IH' in read.tags:
            #         ih = int(read.opt('IH'))
            #     else:
            #         ih = 1

            #     if ih == 1:
            #         single_count += 1
            #     else:
            #         multireads.add(read.qname)

            if coverage:
                mean, stdev, median = calc_coverage(bam, chrom, strand if stranded else None, starts, ends, whitelist, blacklist, library_type=library_type)
                outcols.append(mean)
                outcols.append(stdev)
                outcols.append(median)

            if count > 0:
                if not count in counts_tally:
                    counts_tally[count] = 1
                else:
                    counts_tally[count] += 1

            if callback:
                for callback_cols in callback(bam, count, reads, outcols):
                    tmpcounts.write(count, coding_len, callback_cols)
                    # region_counts.append((count, coding_len, callback_cols))
            else:
                    tmpcounts.write(count, coding_len, outcols)
                # region_counts.append((count, coding_len, outcols))

        if not quiet:
            sys.stderr.write('Calculating normalization...')

        norm_val = None
        norm_val_orig = None

        if norm == 'all':
            norm_val_orig = _find_mapped_count(bam, whitelist, blacklist, quiet)
        elif norm == 'mapped':
            # norm_val_orig = single_count + len(multireads)
            norm_val_orig = total_count
        # elif norm == 'quantile':
        #     norm_val_orig = _find_mapped_count_pcts([x[0] for x in region_counts])
        elif norm == 'median':
            norm_val_orig = ngsutils.support.stats.count_median(counts_tally)
            # norm_val_orig = _find_mapped_count_median([x[0] for x in region_counts])

        if norm_val_orig:
            norm_val = float(norm_val_orig) / 1000000

        if not quiet:
            sys.stderr.write('\n')

        out.write('## %s\n' % (ngsutils.version()))
        out.write('## input%s%s\n' % (' ' if bam.filename else '', bam.filename))
        out.write('## model %s %s\n' % (self.get_name(), self.get_source()))
        out.write('## library_type %s\n' % library_type)
        out.write('## multiple %s\n' % multiple)
        if start_only:
            out.write('## start_only\n')
        if norm_val:
            out.write('## norm %s %s\n' % (norm, float(norm_val_orig)))
            out.write('## CPM-factor %s\n' % norm_val)

        # out.write('#')
        out.write('\t'.join(self.get_headers()))
        out.write('\tlength\tcount')
        if norm_val:
            out.write('\tcount (CPM)')
            if fpkm:
                out.write('\tRPKM')

        if coverage:
            out.write('\tcoverage mean\tcoverage stdev\tcoverage median')

        if self.get_postheaders():
            out.write('\t')
            out.write('\t'.join(self.get_postheaders()))

        out.write('\n')

        for count, coding_len, outcols in tmpcounts.fetch():
            first = True
            for col in outcols:
                if not first:
                    out.write('\t')
                first = False

                if col == '' or col is None:  # this is the marker for the 'count' col
                    out.write('%s' % count)

                    if norm_val:
                        out.write('\t')
                        out.write(str(count / norm_val))
                        if fpkm:
                            out.write('\t')
                            out.write(str(count / (coding_len / 1000.0) / norm_val))

                else:
                    out.write(str(col))

            out.write('\n')
        tmpcounts.close()


def _calc_read_regions(read):
    'Find regions of reference the read covers - breaking on long gaps (N)'
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


def _fetch_reads_excluding(bam, chrom, strand, start, end, multiple, whitelist=None, blacklist=None, library_type='FR'):
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
        frag_strand = None
        if library_type == 'FR':
            if read.is_read2:
                frag_strand = '+' if read.is_reverse else '-'
            else:
                frag_strand = '-' if read.is_reverse else '+'
        elif library_type == 'RF':
            if read.is_read2:
                frag_strand = '-' if read.is_reverse else '+'
            else:
                frag_strand = '+' if read.is_reverse else '-'

        if not strand or strand == frag_strand:
            excl = True
            for s, e in _calc_read_regions(read):
                if start <= s <= end or start <= e <= end:
                    excl = False
                    break
            if excl and not read in reads:
                reads.add(read.qname)
                count += 1
    return count, reads


def _fetch_reads(bam, chrom, strand, starts, ends, multiple, exclusive, whitelist=None, blacklist=None, uniq=False, library_type='FR', start_only=False):
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

    start_only means that only the start of the read is used to determine presence in a start-end region.

    paired-end options:
        if rev_read2 is True, then the reads flagged as read2 will
        have their strand information reversed.

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

                frag_strand = None
                # readpos = None

                if library_type == 'FR':
                    if read.is_read2:
                        frag_strand = '+' if read.is_reverse else '-'
                    else:
                        frag_strand = '-' if read.is_reverse else '+'
                elif library_type == 'RF':
                    if read.is_read2:
                        frag_strand = '-' if read.is_reverse else '+'
                    else:
                        frag_strand = '+' if read.is_reverse else '-'

                if not strand or strand == frag_strand:
                    if start_only:
                        start_ok = False
                        for s1, e1 in zip(starts, ends):
                            if not read.is_reverse:
                                if s1 <= read.pos <= e1:
                                    start_ok = True
                                    break
                            else:
                                if s1 <= read.aend <= e1:
                                    start_ok = True
                                    break

                        if not start_ok:
                            continue

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
                        reads.add(read.qname)

                        ih = 0
                        for tag, val in read.tags:
                            if tag == 'IH':
                                ih = int(val)
                                break
                            elif tag == 'NH':
                                ih = int(val)
                                break

                        if not ih:
                            ih = 1

                        if ih == 1 or multiple == 'complete':
                            count += 1
                        elif multiple == 'partial':
                            count += (1.0 / ih)
                        else:  # multiple = ignore
                            pass
    return count, reads


def calc_coverage(bam, chrom, strand, starts, ends, whitelist, blacklist, library_type='FR'):
    if not chrom in bam.references:
        return 0, 0, 0

    coverage = []
    for start, end in zip(starts, ends):
        for pileup in bam.pileup(chrom, start, end):
            count = 0
            for pileupread in pileup.pileups:
                if blacklist and pileupread.alignment.qname in blacklist:
                    continue
                if not whitelist or pileupread.alignment.qname in whitelist:
                    if strand:
                        frag_strand = None
                        if library_type == 'FR':
                            if pileupread.alignment.is_read2:
                                frag_strand = '-' if not pileupread.alignment.is_reverse else '+'
                            else:
                                frag_strand = '+' if not pileupread.alignment.is_reverse else '-'
                        elif library_type == 'RF':
                            if pileupread.alignment.is_read2:
                                frag_strand = '+' if not pileupread.alignment.is_reverse else '-'
                            else:
                                frag_strand = '-' if not pileupread.alignment.is_reverse else '+'

                        if strand != frag_strand:
                            continue

                    if not pileupread.is_del:
                        count += 1

            coverage.append(count)

    if coverage:
        mean, stdev = ngsutils.support.stats.mean_stdev(coverage)
        median = ngsutils.support.stats.median(coverage)

        return mean, stdev, median
    else:
        return 0, 0, 0


def _find_mapped_count_median(counts):
    '''
    >>> _find_mapped_count_median([10, 20, 30, 10, 30])
    20
    >>> _find_mapped_count_median([10, 20, 30, 10, 30, 0])
    20
    >>> _find_mapped_count_median([10, 20, 30, 10, 30, 22])
    21.0
    '''

    return ngsutils.support.stats.median([x for x in counts if x > 0])


def _find_mapped_count_pcts(counts, min_pct=0.0, max_pct=0.75):
    '''
    >>> _find_mapped_count_pcts([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    280

    >>> _find_mapped_count_pcts([0, 20, 30, 40, 50, 60, 70, 80, 90, 100, 0, 0, 10])
    280

    >>> _find_mapped_count_pcts([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], 0.1, 0.2)
    20

    >>> _find_mapped_count_pcts([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], 0.0, 1.0)
    550
    '''
    counts.sort()
    filtered = [x for x in counts if x > 0]
    acc = 0
    min_i = int(min_pct * len(filtered))
    max_i = int(max_pct * len(filtered))

    for count in filtered[min_i:max_i]:
        acc += count

    return acc


def _find_mapped_count(bam, whitelist=None, blacklist=None, quiet=False):
    '''
    >>> _find_mapped_count(MockBam(['chr1']).add_read('foo1', tid=0, pos=100, cigar='50M').add_read('foo2', tid=0, pos=100, cigar='50M').add_read('foo3', tid=0, pos=100, cigar='50M').add_read('foo4'), quiet=True)
    3
    >>> _find_mapped_count(MockBam(['chr1']).add_read('foo1', tid=0, pos=100, cigar='50M').add_read('foo2', tid=0, pos=100, cigar='50M').add_read('foo3', tid=0, pos=100, cigar='50M').add_read('foo4'), whitelist=['foo1', 'foo2'], quiet=True)
    2
    >>> _find_mapped_count(MockBam(['chr1']).add_read('foo1', tid=0, pos=100, cigar='50M').add_read('foo2', tid=0, pos=100, cigar='50M').add_read('foo3', tid=0, pos=100, cigar='50M').add_read('foo4'), blacklist=['foo1', 'foo2'], quiet=True)
    1
    >>> _find_mapped_count(MockBam(['chr1']).add_read('foo1', tid=0, pos=100, cigar='50M', tags=[('IH', 2)]).add_read('foo1', tid=0, pos=200, cigar='50M', tags=[('IH', 2)]).add_read('foo2', tid=0, pos=100, cigar='50M').add_read('foo3', tid=0, pos=100, cigar='50M').add_read('foo4'), quiet=True)
    3
    '''
    if not quiet:
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
                try:
                    ih = int(read.opt('IH'))
                except:
                    try:
                        ih = int(read.opt('NH'))
                    except:
                        ih = 1
                if ih > 1:
                    multinames.add(read.qname)
    bam.seek(0)
    if not quiet:
        sys.stderr.write("%s mapped reads\n" % mapped_count)
    return mapped_count
