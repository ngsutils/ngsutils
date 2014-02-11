#!/usr/bin/env python
## category General
## desc Build a junction library from FASTA and GTF model
'''
Takes a GTF gene model and a genome FASTA file and produces a splice-
junction library in FASTA format.

This will build junctions of a minimum size with a minimum overlap. Meaning,
if the size of an exon is smaller than the minumum size, the next splice-site
will be appended. This can result in more than one potential junction being
included. This way potential RNA reads can be properly mapped across as many
splice sites as is needed. The longer the reads, the more junctions that may
be potentially included. (Maximum of 5)

This will also take into account the sides of a junction to include so that
sufficient overlap with the junction is established when reads are mapped.
'''

import os
import sys
import pysam
from ngsutils.gtf import GTF
from eta import ETA


def gtf_junctions(gtf, refname, fragment_size, min_size, max_exons=5, known=False, out=sys.stdout, quiet=False, scramble=False, retain_introns=False):
    ref = pysam.Fastafile(refname)

    references = []
    with open('%s.fai' % refname) as f:
        for line in f:
            cols = line.split('\t')
            references.append(cols[0])

    if not quiet:
        eta = ETA(gtf.fsize(), fileobj=gtf)
    else:
        eta = None

    exporter = JunctionExporter(ref, fragment_size, min_size, max_exons, out, scramble)

    for gene in gtf.genes:
        if not gene.chrom in references:
            continue

        if eta:
            eta.print_status(extra='%s:%s %s' % (gene.chrom, gene.start, gene.gene_name))

        if known:
            for txpt in gene.transcripts:
                last = None
                for exon in txpt.exons:
                    if last:
                        exporter.export(gene.chrom, [last, exon])
                    last = exon
        else:
            exons = set()
            for txpt in gene.transcripts:
                for exon in txpt.exons:
                    exons.add(exon)

            exons = list(exons)
            exons.sort()

            if retain_introns:
                exporter.export_retained_introns(gene.chrom, exons, gene.strand)

            if scramble:
                # We can just pretend the transcript is repeated
                # and then let the set take care of removing the duplicates
                exons = exons * 2

            exporter.export(gene.chrom, exons)



    if eta:
        eta.done()
    ref.close()


class JunctionExporter(object):
    def __init__(self, ref, fragment_size, min_size, max_exons, out, scramble=False):
        self.ref = ref
        self.fragment_size = fragment_size
        self.min_size = min_size
        self.max_exons = max_exons
        self.out = out
        self._junctions = set()
        self._cur_chrom = None
        self.scramble = scramble

    def export_retained_introns(self, chrom, exons, strand):
        # Retain introns from both the alt-3 and alt-5 side.
        if chrom != self._cur_chrom:
            self._junctions = set()
            self._cur_chrom = chrom

        for i, (start, end) in enumerate(exons):
            # alt-3' extension (rel to + strand)
            if i < len(exons):
                frag_start = max(start, end - self.fragment_size)
                frag_end = end + self.fragment_size

                name = '%s:%s-%s,%s-%s' % (chrom, frag_start, end, end, frag_end)
                if not name in self._junctions:
                    self._junctions.add(name)
                    seq = self.ref.fetch(chrom, frag_start, frag_end)
                    self.out.write('>%s\n%s\n' % (name, seq))

            # alt-5 extension
            if i > 0:
                frag_start = start - self.fragment_size
                frag_end = min(end, start + self.fragment_size)

                name = '%s:%s-%s,%s-%s' % (chrom, frag_start, start, start, frag_end)
                if not name in self._junctions:
                    self._junctions.add(name)
                    seq = self.ref.fetch(chrom, frag_start, frag_end)
                    self.out.write('>%s\n%s\n' % (name, seq))

    def export(self, chrom, exons):
        if chrom != self._cur_chrom:
            self._junctions = set()
            self._cur_chrom = chrom

        for i, (start, end) in enumerate(exons):
            if i == len(exons) - 1:
                # can't splice the last exon
                continue

            frag_start = start

            if end - start > self.fragment_size:
                frag_start = end - self.fragment_size

            # print '[%s] %s:%s-%s' % (i,gene.chrom,frag_start,end)
            seq3 = self.ref.fetch(chrom, frag_start, end)
            for j in xrange(len(exons) - i - 1):
                # print '   [%s]' % (j+i+1),
                # print '%s-%s' % exons[j+i+1]
                for name, seq in self._extend_junction(seq3, '%s:%s-%s' % (chrom, frag_start, end), chrom, exons[j + i + 1:], end):
                    if not name in self._junctions:
                        self._junctions.add(name)
                        self.out.write('>%s\n%s\n' % (name, seq))

    def _extend_junction(self, seq, name, chrom, exons, anchor_frag_end, counter=1):
        if counter >= self.max_exons:
            return

        start, end = exons[0]
        if not self.scramble and start <= anchor_frag_end:
            # this shouldn't happen in non-scrambled 
            # conditions, but just in case.
             return

        frag_end = end
        if end - start > self.fragment_size:
            frag_end = start + self.fragment_size

        seq5 = self.ref.fetch(chrom, start, frag_end)
        newname = '%s,%s-%s' % (name, start, frag_end)
        newseq = seq + seq5
        if len(newseq) >= self.min_size:
            yield newname, newseq
            return
        elif len(exons) > 1 and (counter + 1) < self.max_exons:
            for i in xrange(1, len(exons)):
                for nn_name, nn_seq in self._extend_junction(newseq, newname, chrom, exons[i:], frag_end, counter + 1):
                    yield nn_name, nn_seq


def usage(msg=None):
    if msg:
        print msg
    print __doc__
    print '''\
Usage: gtfutils junctions {opts} genes.gtf{.gz} genome.fasta

Arguments
  genes.txt       Gene model in GTF format
  genome.fasta    Reference genome in FASTA or RAGZ format
                  (samtools indexed ref.fasta.fai req'd)

Options
  -frag size        Number of bases on either side of the junction to include
                    [default 46]
  -min size         Minimum size of a junction
                    [default 50]
  -known            Only export known junctions

  -scramble         Include potential circular junctions
  -retain-introns   Include retained introns (retains introns from both the 
                    5' and 3' splice side)

                        ________           _______
                    ---|        |---------|       |---
                        --------           -------
                              ******   ++++++

'''
    sys.exit(1)

if __name__ == '__main__':
    gtf = None
    fasta = None
    frag_size = 46
    min_size = 50
    known = False
    scramble = False
    retain_introns = False
    last = None

    for arg in sys.argv[1:]:
        if last == '-frag':
            frag_size = int(arg)
            last = None
        elif last == '-min':
            min_size = int(arg)
            last = None
        elif arg in ['-frag', '-min']:
            last = arg
        elif arg == '-scramble':
            scramble = True
        elif arg == '-retain-introns':
            retain_introns = True
        elif arg == '-known':
            known = True
        elif arg == '-h':
            usage()
        elif gtf is None and os.path.exists(arg):
            gtf = arg
        elif fasta is None and os.path.exists(arg):
            if not os.path.exists('%s.fai' % arg):
                usage('Missing .fai FASTA index for file: %s' % arg)
            fasta = arg

    if not gtf or not fasta:
        usage()

    if known and scramble:
        usage("You can not use both -known and -scramble at the same time!")

    gtf_junctions(GTF(gtf), fasta, frag_size, min_size, known=known, scramble=scramble, retain_introns=retain_introns)
