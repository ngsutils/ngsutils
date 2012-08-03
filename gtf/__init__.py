'''
GTF gene models

See: http://mblab.wustl.edu/GTF22.html

GTF/GFF is a 1-based format, meaning that start sites start counting at 1. We
read in those start sites and subtract one so that the internal representation
is 0-based.

All positions returned are 0-based
'''


import sys
import os
from support.ngs_utils import  gzip_opener
from support import symbols
from support.eta import ETA


class GTF(object):
    def __init__(self, filename):
        self._genes = {}
        self._pos = 0
        self._sources = []
        self._gene_order = []
        warned = False

        sys.stderr.write('Reading GTF file...\n')
        with gzip_opener(filename) as f:
            eta = ETA(os.stat(filename).st_size, fileobj=f)
            for line in f:
                try:
                    idx = line.find('#')
                    if idx > -1:
                        if idx == 0:
                            continue
                        line = line[:-idx]
                    chrom, source, feature, start, end, score, strand, frame, attrs = line.rstrip().split('\t')
                    source = symbols[source]
                    start = int(start) - 1  # Note: 1-based
                    end = int(end)
                    attributes = {}
                    for key, val in [x.split(' ', 1) for x in [x.strip() for x in attrs.split(';')] if x]:
                        if val[0] == '"' and val[-1] == '"':
                            val = val[1:-1]
                        attributes[key] = val

                    gid = None
                    if 'isoform_id' in attributes:
                        gid = attributes['isoform_id']
                    else:
                        gid = attributes['gene_id']
                        if not warned:
                            sys.stderr.write('\nGTF file missing isoform annotation! Each transcript will be treated separately. (%s)\n' % gid)
                            warned = True

                    eta.print_status(extra=gid)
                except:
                    import traceback
                    sys.stderr.write('Error parsing line:\n%s\n' % line)
                    traceback.print_exc()
                    sys.exit(1)

                if not gid in self._genes or chrom != self._genes[gid].chrom:
                    self._genes[gid] = _GTFGene(gid, chrom, source, **attributes)
                    self._gene_order.append((chrom, start, gid))

                self._genes[gid].add_feature(attributes['transcript_id'], feature, start, end, strand)

        eta.done()
        self._gene_order.sort()

    def fsize(self):
        return len(self._genes)

    def tell(self):
        return self._pos

    @property
    def genes(self):
        self._pos = 0
        for chrom, start, gene_id in self._gene_order:
            yield self._genes[gene_id]
            self._pos += 1


class _GTFGene(object):
    """
    Stores info for a single gene_id

    A gene consists of one or more transcripts. It *must* have a gene_id, and chrom
    It may also contain a gene_name and an isoform_id.

    The gid is the 'gene_id' from the GTF file *unless* there is an isoform_id
    attribute. If 'isoform_id' is present, this is used for the gid.


    """

    def __init__(self, gid, chrom, source, gene_id, gene_name=None, isoform_id=None, **extras):
        self.gid = gid
        self.gene_id = gene_id
        self.gene_name = gene_name if gene_name else gene_id
        self.isoform_id = isoform_id if isoform_id else gene_id

        self.chrom = chrom
        self.source = source

        self._transcripts = {}
        self._regions = []

        self.start = None
        self.end = None
        self.strand = None

    @property
    def transcripts(self):
        for t in self._transcripts:
            yield self._transcripts[t]

    def add_feature(self, transcript_id, feature, start, end, strand):
        if not transcript_id in self._transcripts:
            self._transcripts[transcript_id] = _GTFTranscript(transcript_id, strand)

        t = self._transcripts[transcript_id]

        if feature == 'exon':
            t.exons.append((start, end))

            if self.start is None or start < self.start:
                self.start = start
            if self.end is None or end > self.end:
                self.end = end
            if self.strand is None:
                self.strand = strand

            if t.start is None or start < t.start:
                t.start = start
            if t.end is None or end > t.end:
                t.end = end

        elif feature == 'CDS':
            t.cds.append((start, end))
        elif feature == 'start_codon':
            t.start_codon = (start, end)
        elif feature == 'stop_codon':
            t.stop_codon = (start, end)
        else:
            # this is an unsupported feature - possibly add a debug message
            pass

    @property
    def regions(self):
        # these are potentially memory-intensive, so they are calculated on the
        # fly when needed.
        if not self._regions:
            all_starts = []
            all_ends = []
            tids = []

            for tid in self._transcripts:
                tids.append(tid)
                starts = []
                ends = []
                for start, end in self._transcripts[tid].exons:
                    starts.append(start)
                    ends.append(end)
                all_starts.append(starts)
                all_ends.append(ends)

            self._regions = calc_regions(self.start, self.end, tids, all_starts, all_ends)

        i = 0
        for start, end, const, names in self._regions:
            i += 1
            yield (i, start, end, const, names)


def calc_regions(txStart, txEnd, kg_names, kg_starts, kg_ends):
    '''
        This takes a list of start/end positions (one set per isoform)

        It splits these into regions by assigning each base in the
        txStart-txEnd range a number.  The number is a bit-mask representing
        each isoform that includes that base in a coding region.  Each
        isoform has a number 2^N, so the bit-mask is just a number.

        Next, it scans the map for regions with the same bit-mask.
        When a different bitmask is found, the previous region (if not intron)
        is added to the list of regions.

        Returns a list of tuples:
        (start,end,is_const,names) where names is a comma-separated string
                                   of accns that make up the region

    '''

    map = [0, ] * (txEnd - txStart + 1)
    mask = 1

    mask_start_end = {}
    mask_names = {}

    for name, starts, ends in zip(kg_names, kg_starts, kg_ends):
        mask_start = None
        mask_end = None
        mask_names[mask] = name

        for start, end in zip(starts, ends):
            if not mask_start:
                mask_start = int(start)
            mask_end = int(end)

            for i in xrange(int(start) - txStart, int(end) - txStart + 1):
                map[i] = map[i] | mask

        mask_start_end[mask] = (mask_start, mask_end)
        mask = mask * 2

    last_val = 0
    regions = []
    region_start = 0

    def _add_region(start, end, value):
        rstart = start + txStart
        rend = end + txStart
        const = True
        names = []

        '''
        This code only calls a region alt, if the transcript actually spans this region.

        example - two transcripts:

        1)        xxxxxx------xxxx------xxx---xxxx--xxxxxxxxx
        2)           xxx------xxxx--xxx-------xxxx--xxxxxxxxxxxxx
        const/alt cccccc------cccc--aaa-aaa---cccc--ccccccccccccc

        I'm not sure if this is a good idea or not... What this *should* be is:

        1)        xxxxxx------xxxx------xxx---xxxx----xxxxxxxxx
        2)           xxx------xxxx--xxx-------xxxx----xxxxxxxxxxxxx
        3)           xxx------xxxx--xxx-------xxxxxa
        const/alt aaaccc------cccc--aaa-aaa---ccccaa--cccccccccaaaa

        Where alt-5 and alt-3 are only dictated by transcripts that contain them...
        There should be extra code to handle this situation...

        ## TODO: Add alt-5/alt-3 code
        '''

        for mask in mask_start_end:
            mstart, mend = mask_start_end[mask]
            if rstart >= mstart and rend <= mend:
                if value & mask == 0:
                    const = False
                else:
                    names.append(mask_names[mask])

        regions.append((rstart, rend, const, ','.join(names)))

    for i in xrange(0, len(map)):
        if map[i] == last_val:
            continue

        if last_val:
            _add_region(region_start, i - 1, last_val)

            region_start = i - 1
        else:
            region_start = i

        last_val = map[i]

    if last_val:
        _add_region(region_start, i, last_val)  # extend by one...

    return regions


class _GTFTranscript(object):
    def __init__(self, transcript_id, strand):
        self.transcript_id = transcript_id
        self.strand = strand
        self.exons = []
        self.cds = []
        self.start_codon = None
        self.stop_codon = None

        self.start = None
        self.end = None
