#!/usr/bin/env python
## category Conversion
## desc Extract BED regions from a reference FASTA file
'''
Extract BED regions from a reference FASTA file
'''

import sys
import os
from ngsutils.bed import BedFile
from ngsutils.support import revcomp
import pysam


def bed_tofasta(bed, ref_fasta, min_size=50, stranded=True, out=sys.stdout):
    if not os.path.exists('%s.fai' % ref_fasta):
        pysam.faidx(ref_fasta)

    fasta = pysam.Fastafile(ref_fasta)

    refs = set()
    with open('%s.fai' % ref_fasta) as f:
        for line in f:
            refs.add(line.split('\t')[0].strip())

    for region in bed:
        if region.end - region.start >= min_size and region.chrom in refs:
            seq = fasta.fetch(region.chrom, region.start, region.end)
            if stranded and region.strand:
                if region.strand == '-':
                    seq = revcomp(seq)
                out.write('>%s:%d-%d[%s]\n%s\n' % (region.chrom, region.start, region.end, region.strand, seq))
            else:
                out.write('>%s:%d-%d\n%s\n' % (region.chrom, region.start, region.end, seq))

    fasta.close()


def usage():
    print __doc__
    print """\
Usage: bedutils tofasta {-min size} {-ns} bedfile ref.fasta

Outputs the sequences of each BED region to FASTA format.

Option:
-min  The minumum size of a region
-ns   Ignore the strand of a region (always return seq from the + strand)
"""

if __name__ == "__main__":

    min_size = 50
    bed = None
    ref = None
    stranded = True

    last = None
    for arg in sys.argv[1:]:
        if last == '-min':
            min_size = int(arg)
            last = None
        elif arg in ['-min']:
            last = arg
        elif arg == '-ns':
            stranded = False
        elif not bed and os.path.exists(arg):
            bed = arg
        elif not ref and os.path.exists(arg):
            ref = arg

    if not bed or not ref:
        usage()
        sys.exit(1)

    bed_tofasta(BedFile(bed), ref, min_size=min_size, stranded=stranded)
