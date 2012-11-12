#!/usr/bin/env python
## category General
## desc Extracts reads based on regions in a BED file
"""
Extracts reads based on regions in a BED file

Given a BAM file and a BED file, this script will extract only reads that
map to the given BED regions.

Specifically, if a read starts or ends within a BED region, it is extracted.
However, if a read completely spans a region (starting before and ending
after), it is ignored.

This tends to be faster than running 'bamutils filter' because this extracts
reads from a region, without iterating over all of the reads in the BAM file.
But, if you are already filtering the BAM file for other criteria, it is
probably worth it to just to the filtering with your other criteria.
"""

import os
import sys
from eta import ETA
import pysam

from ngsutils.bam import read_alignment_fragments_gen


def usage():
    print __doc__
    print """
Usage: bamutils extract {opts} in.bam out.bam regions.bed

Options:
  -ns    Ignore strandedness of reads and regions
"""
    sys.exit(1)


def bam_extract(infile, outfile, bedfile, nostrand=False):
    with open(bedfile) as f:
        bamfile = pysam.Samfile(infile, "rb")
        outfile = pysam.Samfile(outfile, "wb", template=bamfile)
        eta = ETA(os.stat(bedfile).st_size, fileobj=f)

        passed = 0

        for line in f:
            eta.print_status(extra="extracted:%s" % (passed))
            if line[0] == '#':
                continue
            cols = line.strip().split('\t')
            if not cols:
                continue

            chrom = cols[0]

            if not chrom in bamfile.references:
                continue

            start = int(cols[1])
            end = int(cols[2])

            if not nostrand and len(cols) >= 6:
                strand = cols[5]
            else:
                strand = None

            for read in bam_extract_reads(bamfile, chrom, start, end, strand):
                outfile.write(read)
                passed += 1

        eta.done()
        bamfile.close()
        outfile.close()

    sys.stderr.write("%s extracted\n" % (passed,))


def bam_extract_reads(bamfile, chrom, start, end, strand=None):
    for read in bamfile.fetch(chrom, start, end):
        if strand is None or (read.is_reverse and strand == '-') or (not read.is_reverse and strand == '+'):
            good = False
            for frag_start, frag_end in read_alignment_fragments_gen(read):
                if start <= frag_start <= end or start <= frag_end <= end:
                    good = True
                    break
            if good:
                yield read


if __name__ == '__main__':
    infile = None
    outfile = None
    bedfile = None
    nostrand = False

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif arg == '-ns':
            nostrand = True
        elif not infile and os.path.exists(arg):
            infile = arg
        elif not outfile:
            outfile = arg
        elif not bedfile and os.path.exists(arg):
            bedfile = arg
        else:
            print "Unknown argument: %s" % arg
            usage()

    if not infile or not outfile or not bedfile:
        usage()
    else:
        bam_extract(infile, outfile, bedfile, nostrand)
