#!/usr/bin/env python
## category Conversion
## desc Extract BED regions from a reference FASTA file
'''
Extract BED regions from a reference FASTA file.

Note: Sequences that are extracted will be in the same orientation as the
BED region, unless the {-ns} option is given.
'''

import sys
import os
from ngsutils.bed import BedFile
from ngsutils.support import revcomp
import pysam


def bed_tofasta(bed, ref_fasta, min_size=50, stranded=True, include_name=False, out=sys.stdout):
    if not os.path.exists('%s.fai' % ref_fasta):
        pysam.faidx(ref_fasta)

    fasta = pysam.Fastafile(ref_fasta)

    refs = set()
    with open('%s.fai' % ref_fasta) as f:
        for line in f:
            refs.add(line.split('\t')[0].strip())

    name = ''
    for region in bed:
        if include_name:
            name = '%s|' % (region.name.strip())

        if region.end - region.start >= min_size and region.chrom in refs:
            seq = fasta.fetch(region.chrom, region.start, region.end)
            if stranded and region.strand:
                if region.strand == '-':
                    seq = revcomp(seq)
                out.write('>%s%s:%d-%d[%s]\n%s\n' % (name, region.chrom, region.start, region.end, region.strand, seq))
            else:
                out.write('>%s%s:%d-%d%s\n%s\n' % (name, region.chrom, region.start, region.end, seq))

    fasta.close()


def usage():
    print __doc__
    print """\
Usage: bedutils tofasta {-min size} {-name} {-ns} bedfile ref.fasta

Outputs the sequences of each BED region to FASTA format.

Option:
-min    The minumum size of a region

-name   Include the name field of the BED region in the FASTA sequence name
        If used, the final name will be in the form:
            name|chrX:start-end[strand]

        The default is to not include the BED region name (only the genomic
        coordinates will be exported).

-ns     Ignore the strand of a region (always return seq from the + strand)
"""

if __name__ == "__main__":

    min_size = 50
    bed = None
    ref = None
    stranded = True
    include_name = False

    last = None
    for arg in sys.argv[1:]:
        if last == '-min':
            min_size = int(arg)
            last = None
        elif arg in ['-min']:
            last = arg
        elif arg == '-name':
            include_name = True
        elif arg == '-ns':
            stranded = False
        elif not bed and os.path.exists(arg):
            bed = arg
        elif not ref and os.path.exists(arg):
            ref = arg

    if not bed or not ref:
        usage()
        sys.exit(1)

    bed_tofasta(BedFile(bed), ref, min_size=min_size, stranded=stranded, include_name=include_name)
