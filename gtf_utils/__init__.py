'''
GTF Utility classes

See: http://mblab.wustl.edu/GTF22.html

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
                    start = int(start)  # Note: 1-based
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
    """

    def __init__(self, gid, chrom, source, gene_id, gene_name=None, isoform_id=None, **extras):
        self.gid = gid
        self.gene_id = gene_id
        self.gene_name = gene_name if gene_name else gene_id
        self.isoform_id = isoform_id if isoform_id else gene_id

        self.chrom = chrom
        self.source = source

        self._transcripts = {}
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
