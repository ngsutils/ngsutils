#!/usr/bin/env python
## category General
## desc Given a number of BED files, calculate the number of samples that overlap regions in a reference BED file
'''
Given a BED file containing genomic regions of interest (reference) and a set
of other BED files (tabix-indexed), calculate how many samples also
contain the region of interest.

Example:
ref           (A)========        (B)============================
sample 1       ---    ----              -----  -- ------   ----
sample 2           ----           ------------
sample 3                                                     -----
sample 4          ----

Would return:
        sample 1        sample 2        sample 3        sample 4        total
ref A      2               1               0               1              3
ref B      4               1               1               0              3

'''

import os
import sys
import ngsutils.support.ngs_utils
from ngsutils.bed import BedFile


def usage():
    print __doc__
    print """\
Usage: bedutils refcount {-ns} ref.bed sample1.bed sample2.bed ...

Where the sample BED files are Tabix indexed (sample1.bed.tgz.tbi exists).
The sample names are calculated using the names of the files.

The reference BED file doesn't need to be Tabix indexed (it is read directly).

Options:
  -ns   Ignore strandedness in counts
"""


def bed_refcount(refbed, bedfiles, stranded=True, out=sys.stdout):
    samples = ngsutils.support.ngs_utils.filenames_to_uniq([os.path.basename(x.filename) for x in bedfiles])

    out.write('#chrom\tstart\tend\tname\tscore\tstrand')
    for sample in samples:
        out.write('\t')
        out.write(sample)
    out.write('\tpresent\n')

    for region in refbed:
        out.write(region)
        present = 0
        for bed in bedfiles:
            count = 0
            for sample_region in bed.fetch(region.chrom, region.start, region.end, region.strand if stranded else None):
                count += 1

            out.write('\t')
            out.write(str(count))
            if count > 0:
                present += 1

        out.write('\t%s\n' % present)


if __name__ == '__main__':
    refname = None
    bedfiles = []
    stranded = True
    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif arg == '-ns':
            stranded = False
        elif not refname and os.path.exists(arg):
            refname = arg
        elif os.path.exists(arg):
            bedfiles.append(BedFile(arg))
        else:
            print "Bad argument: %s" % arg
            usage()
            sys.exit(1)

    if not refname or not bedfiles:
        usage()
        sys.exit(1)

    bed_refcount(BedFile(refname), bedfiles, stranded=stranded)
