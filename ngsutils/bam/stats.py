#!/usr/bin/env python
## category General
## desc Calculates simple stats for a BAM file
"""
Calculates simple stats for a BAM file
"""

import os
import sys
from ngsutils.bam import read_calc_mismatches, bam_iter, bam_open
from ngsutils.gtf import GTF
from ngsutils.support.regions import RegionTagger
from ngsutils.support.stats import counts_mean_stdev


class FeatureBin(object):
    '''track feature stats'''
    def __init__(self, tag):
        spl = tag.split(':')
        self.tag = spl[0]
        self.asc = True
        if len(spl) > 1:
            self.asc = True if spl[1] == '+' else False

        if tag in ["LENGTH", "LEN"]:
            self.asc = False

        self.bins = {}
        self._keys = []
        self._min = None
        self._max = None
        self.__cur_pos = -1
        self.missing = 0

    def __iter__(self):
        self._keys.sort()
        if not self.asc:
            self._keys.reverse()

        self.__cur_pos = -1
        return self

    def next(self):
        self.__cur_pos += 1
        if len(self._keys) <= self.__cur_pos:
            raise StopIteration
        else:
            return (self._keys[self.__cur_pos], self.bins[self._keys[self.__cur_pos]])

    @property
    def mean(self):
        acc = 0
        count = 0
        for k in self.bins:
            try:
                acc += (k * self.bins[k])
                count += self.bins[k]
            except:
                return 0

        if count > 0:
            return float(acc) / count
        else:
            return 'Not present'

    @property
    def max(self):
        if self._max is None:
            return 'Not present'

        return self._max

    def add(self, read):
        if self.tag in ['LENGTH', 'LEN']:
            val = len(read.seq)

        elif self.tag == 'MAPQ':
            val = read.mapq

        elif self.tag == 'MISMATCH':
            val = read_calc_mismatches(read)

        else:
            try:
                val = read.opt(self.tag)
            except KeyError:
                self.missing += 1
                return

        if not val in self.bins:
            self.bins[val] = 0
            self._keys.append(val)

        self.bins[val] += 1
        if not self._min or self._min > val:
            self._min = val
        if not self._max or self._max < val:
            self._max = val


def usage():
    print __doc__
    print """
Usage: bamutils stats {options} file1.bam {file2.bam...}

If a region is given, only reads that map to that region will be counted.
Regions should be be in the format: 'ref:start-end' or 'ref:start' using
1-based start coordinates.

Options:
    -all    Show the stats for all fragments (defaults to just the first fragment)

    -region chrom:start-end
            Only calculate statistics for this region

    -tags tag_name{:sort_order},tag_name{:sort_order},...

            For each tag that is given, the values for that tag will be
            tallied for all reads. Then a list of the counts will be presented
            along with the mean and maximum values. The optional sort order
            should be either '+' or '-' (defaults to +).

            There are also special case tags that can be used as well:
                MAPQ     - use the mapq score
                LENGTH   - use the length of the read
                MISMATCH - use the mismatch score (# mismatches) + (# indels)
                           where indels count for 1 regardless of length

                           Note: this requires the 'NM' tag (edit distance)
                           to be present

            Common tags:
                AS    Alignment score
                IH    Number of stored alignments in file for a read
                NH    Number of reported alignments for a read
                NM    Edit distance (each indel counts as many as its length)

            For example, to tally the "IH" tag (number of alignments) and the
            read length:
                -tags IH,LENGTH

    -delim char

            If delimiter is given, the reference names are split by this
            delimiter and only the first token is summarized.

    -gtf model.gtf

            If a GTF gene model is given, counts corresponding to exons,
            introns, promoters, junctions, intergenic, and mitochondrial
            regions will be calculated.

            Note: For paired-end reads, only the first fragment is counted
                  regardless of the {-all} option above
"""
    sys.exit(1)

flag_descriptions = {
    0x1: 'Multiple fragments',
    0x2: 'All fragments aligned',
    0x4: 'Unmapped',
    0x8: 'Next unmapped',
    0x10: 'Reverse complimented',
    0x20: 'Next reverse complimented',
    0x40: 'First fragment',
    0x80: 'Last fragment',
    0x100: 'Secondary alignment',
    0x200: 'QC Fail',
    0x400: 'PCR/Optical duplicate'
}


class FlagCounts(object):
    def __init__(self):
        self.counts = {}
        for fd in flag_descriptions:
            self.counts[fd] = 0

    def add(self, flag):
        for fd in flag_descriptions:
            if (fd & flag) > 0:
                self.counts[fd] += 1


class BamStats(object):
    def __init__(self, bamfile, gtf=None, region=None, delim=None, tags=[], show_all=False):
        regiontagger = None
        flag_counts = FlagCounts()

        ref = None
        start = None
        end = None

        if gtf:
            regiontagger = RegionTagger(gtf, bamfile.references, only_first_fragment=True)

        if region:
            ref, startend = region.rsplit(':', 1)
            if '-' in startend:
                start, end = [int(x) for x in startend.split('-')]
                start = start - 1
                sys.stderr.write('Region: %s:%s-%s\n' % (ref, start + 1, end))
            else:
                start = int(startend) - 1
                end = int(startend)
                sys.stderr.write('Region: %s:%s\n' % (ref, start + 1))

        total = 0
        mapped = 0
        unmapped = 0
        
        tlen_counts = {}

        names = set()
        refs = {}

        tagbins = {}
        for tag in tags:
            tagbins[tag] = FeatureBin(tag)

        for rname in bamfile.references:
            if delim:
                refs[rname.split(delim)[0]] = 0
            else:
                refs[rname] = 0

        # setup region or whole-file readers
        def _foo1():
            for read in bamfile.fetch(ref, start, end):
                yield read

        def _foo2():
            for read in bam_iter(bamfile):
                yield read

        if region:
            read_gen = _foo1
        else:
            read_gen = _foo2

        has_ih = True
        has_nh = True

        try:
            for read in read_gen():
                if not show_all and read.is_paired and not read.is_read1:
                    # only operate on the first fragment
                    continue

                try:
                    if has_ih and read.opt('IH') > 1:
                        if read.qname in names:
                            # reads only count once for this...
                            continue
                        names.add(read.qname)
                except KeyError:
                    if not read.is_unmapped:
                        has_ih = False
                    #missing IH tag - ignore
                    pass

                try:
                    if has_nh and read.opt('NH') > 1:
                        if read.qname in names:
                            # reads only count once for this...
                            continue
                        names.add(read.qname)
                except KeyError:
                    if not read.is_unmapped:
                        has_nh = False
                    #missing NH tag - ignore
                    pass

                flag_counts.add(read.flag)

                total += 1
                if read.is_unmapped:
                    unmapped += 1
                    continue

                mapped += 1

                if read.is_proper_pair and read.tid == read.mrnm:
                    # we don't care about reads that don't map to the same reference

                    # note: this doesn't work for RNA mapped to a reference genome...
                    # for RNA, you'd need to map to a transcript library (refseq) to get
                    # an accurate template length
                    if read.is_reverse:
                        k = -read.tlen
                    else:
                        k = read.tlen

                    if not k in tlen_counts:
                        tlen_counts[k] = 1
                    else:
                        tlen_counts[k] += 1

                if delim:
                    refs[bamfile.getrname(read.rname).split(delim)[0]] += 1
                else:
                    refs[bamfile.getrname(read.rname)] += 1

                if regiontagger:
                    regiontagger.add_read(read, bamfile.getrname(read.rname))

                for tag in tagbins:
                    tagbins[tag].add(read)

        except KeyboardInterrupt:
            sys.stderr.write('*** Interrupted - displaying stats up to this point! ***\n\n')

        self.total = total
        self.mapped = mapped
        self.unmapped = unmapped
        self.flag_counts = flag_counts
        self.tagbins = tagbins
        self.refs = refs
        self.regiontagger = regiontagger
        self.tlen_counts = tlen_counts

    def distribution_gen(self, tag):
        acc = 0.0
        for val, count in self.tagbins[tag]:
            acc += count
            pct = acc / self.mapped
            yield (val, count, pct)


def bam_stats(infiles, gtf_file=None, region=None, delim=None, tags=[], show_all=False):
    if gtf_file:
        gtf = GTF(gtf_file)
    else:
        gtf = None

    sys.stderr.write('Calculating Read stats...\n')

    stats = [BamStats(bam_open(x), gtf, region, delim, tags, show_all=show_all) for x in infiles]

    sys.stdout.write('\t')
    for fname, stat in zip(infiles, stats):
        sys.stdout.write('%s\t\t' % fname)
    sys.stdout.write('\n')

    sys.stdout.write('Reads:\t')
    for stat in stats:
        sys.stdout.write('%s\t\t' % stat.total)
    sys.stdout.write('\n')

    sys.stdout.write('Mapped:\t')
    for stat in stats:
        sys.stdout.write('%s\t\t' % stat.mapped)
    sys.stdout.write('\n')

    sys.stdout.write('Unmapped:\t')
    for stat in stats:
        sys.stdout.write('%s\t\t' % stat.unmapped)
    sys.stdout.write('\n')

    sys.stdout.write('\nFlag distribution\n')
    validflags = set()
    maxsize = 0
    for flag in flag_descriptions:
        for stat in stats:
            if stat.flag_counts.counts[flag] > 0:
                validflags.add(flag)
                maxsize = max(maxsize, len(flag_descriptions[flag]))

    for flag in sorted(validflags):
        sys.stdout.write("[0x%03x] %-*s" % (flag, maxsize, flag_descriptions[flag]))
        for stat in stats:
            sys.stdout.write('\t%s\t%0.2f%%' % (stat.flag_counts.counts[flag], (float(stat.flag_counts.counts[flag]) * 100 / stat.total)))
        sys.stdout.write('\n')
    sys.stdout.write('\n')

    if stats[0].tlen_counts:
        sys.stdout.write('Template length:')
        for stat in stats:
            mean, stdev = counts_mean_stdev(stat.tlen_counts)
            sys.stdout.write('\t%0.2f\t+/- %0.2f' % (mean, stdev))
        sys.stdout.write('\n')
    sys.stdout.write('\n')

    stat_tags = {}
    for tag in stats[0].tagbins:
        stat_tags[tag] = []
        for stat in stats:
            stat_tags[tag].append(stat.tagbins[tag])

    for tag in stat_tags:
        asc = stats[0].tagbins[tag].asc
        sys.stdout.write("Ave %s:" % tag)
        for i, tagbin in enumerate(stat_tags[tag]):
            sys.stdout.write('\t%s' % tagbin.mean)
            if i != len(stats):
                sys.stdout.write('\t')
        sys.stdout.write('\n')

        sys.stdout.write("Max %s:" % tag)
        for i, tagbin in enumerate(stat_tags[tag]):
            sys.stdout.write('\t%s' % tagbin.max)
            if i != len(stats):
                sys.stdout.write('\t')
        sys.stdout.write('\n')

        sys.stdout.write('%s distribution:\n' % tag)

        gens = []
        gen_vals = []
        last_pcts = []

        for stat in stats:
            gens.append(stat.distribution_gen(tag))
            gen_vals.append(None)
            last_pcts.append(0.0)

        good = True

        last = None

        while good:
            good = False
            for i, stat in enumerate(stats):
                if not gen_vals[i]:
                    try:
                        gen_vals[i] = gens[i].next()
                    except StopIteration:
                        pass
            vals = [tup[0] for tup in gen_vals if tup]
            if not vals:
                continue
            if asc:
                minval = min(vals)
            else:
                minval = max(vals)

            if last:
                if asc:
                    last += 1
                    # fill in missing values
                    while last < minval:
                        sys.stdout.write('%s' % last)
                        for i, stat in enumerate(stats):
                            sys.stdout.write('\t0\t%s' % last_pcts[i])
                        sys.stdout.write('\n')
                        last += 1
                else:
                    last -= 1
                    # fill in missing values
                    while last > minval:
                        sys.stdout.write('%s' % last)
                        for i, stat in enumerate(stats):
                            sys.stdout.write('\t0\t%s' % last_pcts[i])
                        sys.stdout.write('\n')
                        last -= 1

            last = minval
            sys.stdout.write(str(minval))

            for i, tup in enumerate(gen_vals):
                if tup and tup[0] == minval:
                    sys.stdout.write('\t%s\t%s' % (tup[1], tup[2]))
                    last_pcts[i] = tup[2]
                    gen_vals[i] = None
                    good = True
                else:
                    sys.stdout.write('\t0\t%s' % (last_pcts[i]))
            sys.stdout.write('\n')
        sys.stdout.write('\n')

    sys.stdout.write('Reference counts')
    for stat in stats:
        sys.stdout.write('\tcount\t')
    sys.stdout.write('\n')
    for k in sorted([x for x in stats[0].refs]):
        sys.stdout.write('%s' % k)
        for stat in stats:
            sys.stdout.write('\t%s\t' % stat.refs[k])
        sys.stdout.write('\n')

    if gtf_file:
        sys.stdout.write('Mapping regions')
        for stat in stats:
            sys.stdout.write('\tcount\tCPM')
        sys.stdout.write('\n')
        sorted_keys = [x for x in stats[0].regiontagger.counts]
        sorted_keys.sort()
        for k in sorted_keys:
            sys.stdout.write('%s' % k)
            for stat in stats:
                sys.stdout.write('\t%s\t%s' % (stat.regiontagger.counts[k], float(stat.regiontagger.counts[k]) / stat.mapped / 1000000))
            sys.stdout.write('\n')


if __name__ == '__main__':
    infiles = []
    gtf = None
    region = None
    delim = None
    show_all = False
    tags = []

    last = None
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-region':
            region = arg
            last = None
        elif last == '-gtf':
            gtf = arg
            last = None
        elif last == '-delim':
            delim = arg
            last = None
        elif last == '-tags':
            tags = arg.split(',')
            last = None
        elif arg == '-all':
            show_all = True
        elif arg in ['-gtf', '-delim', '-tags', '-region']:
            last = arg
        elif os.path.exists(arg):
            infiles.append(arg)
        else:
            sys.stderr.write('Unknown option: %s\n' % arg)
            usage()

    if not infiles:
        usage()
    else:
        bam_stats(infiles, gtf, region, delim, tags, show_all=show_all)
