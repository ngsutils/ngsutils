#!/usr/bin/env python
## category General
## desc Combine multiple BAM files together (taking best-matches)
"""
Combine multiple BAM files together (taking best-matches)

Given a number of BAM files, this script will merge them together, taking
only the best matches.  There can be any number of files, but the BAM header
will be taken from the first one.  The input files should be sorted by read
name, or at least have reads in the same order.

The first input file should have a record for every read in the other files.
However, the secondary files *may* have missing lines, so long as they are in
the same order as the first file.

The value of the attribute/tag given will be used to determine which reads
should be kept and which should be discarded. The tag should be a numberic
(int/float) type. This defaults to 'AS'.

Additionally, each file can have more than one record for each read, but they
should all have the same value for the tag used in determining which reads to
keep. For example, if the AS tag is used (default), then each read in a file
should have the same AS value. Reads in different files will have different
values.

"""

import os
import sys
import pysam


def usage():
    print __doc__
    print """
Usage: bamutils merge {opts} out.bam in1.bam in2.bam ...

Options
  -tag VAL    Tag to use to determine from which file reads will be taken
              (must be type :i or :f) [default: AS]
  -discard    Discard reads that aren't mapped in any file.
"""
    sys.exit(1)


def bam_reads_batch(bam):
    reads = []
    last = None
    for read in bam:
        if last and read.qname != last:
            yield reads
            reads = []
        last = read.qname
        reads.append(read)

    if reads:
        yield reads


def bam_merge(fname, infiles, tag='AS', discard=False, quiet=False):
    bams = []
    last_reads = []
    bamgens = []
    counts = []
    unmapped = 0

    for infile in infiles:
        bam = pysam.Samfile(infile, "rb")
        last_reads.append(None)
        bams.append(bam)
        counts.append(0)
        bamgens.append(bam_reads_batch(bam))

    outfile = pysam.Samfile('%s.tmp' % fname, "wb", template=bams[0])

    while True:
        found = False
        for i, bamgen in enumerate(bamgens):
            if last_reads[i] == None:
                try:
                    last_reads[i] = bamgen.next()
                    if last_reads[i]:
                        found = True
                except:
                    pass
            else:
                found = True
        if not found:
            break

        best_val = None
        best_reads = None
        best_source = 0

        first_group = last_reads[0]
        for i in xrange(len(last_reads)):
            if not last_reads[i]:
                continue

            match = False
            for read in last_reads[i]:
                if read.qname == first_group[0].qname:
                    match = True
                    if not read.is_unmapped:
                        tag_val = int(read.opt(tag))
                        if not best_val or tag_val > best_val:
                            best_val = tag_val
                            best_reads = last_reads[i]
                            best_source = i
                            break
            if match:
                last_reads[i] = None

        if best_reads:
            counts[best_source] += 1
            for read in best_reads:
                outfile.write(read)
        else:
            unmapped += 1

            if not discard:
                outfile.write(first_group[0])

    if not quiet:
        for fn, cnt in zip(infiles, counts):
            print "%s\t%s" % (fn, cnt)
        print "unmapped\t%s" % unmapped

    outfile.close()
    for bam in bams:
        bam.close()

    os.rename('%s.tmp' % fname, fname)

if __name__ == '__main__':
    infiles = []
    outfile = None
    last = None
    discard = False
    tag = 'AS'

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-tag':
            tag = arg
            last = None
        elif arg in ['-tag']:
            last = arg
        elif arg == '-discard':
            discard = True
        elif not outfile:
            outfile = arg
        elif os.path.exists(arg):
            infiles.append(arg)

    if not infiles or not outfile:
        usage()
    else:
        bam_merge(outfile, infiles, tag, discard)
