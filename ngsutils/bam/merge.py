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
should be kept and which should be discarded. The tag should be a numeric
(int/float) type. More than one tag may be used. This defaults to ['AS+', 'NM-'].

Additionally, each file can have more than one record for each read, that may
all have the same value for the tag used in determining which reads to keep.
For example, if the AS tag is used (default), then each read in a file
may have the same AS value. In this case, all reads with the best AS score
will be kept.

"""

import os
import sys
import pysam
import ngsutils.bam


def usage():
    print __doc__
    print """
Usage: bamutils merge {opts} out.bam in1.bam in2.bam ...

Options
  -tag VAL    Tag to use to determine from which file reads will be taken.
              (must be type :i or :f) You may have more than one of these,
              in which case they will be sorted in order. You can add a +/-
              at the end of the name to signify sort order (asc/desc). 
              [default: AS+, NM-]

  -discard    Discard reads that aren't mapped in any file.

  -keepall    Keep all mappings for each read, not just the best one.
              (Note: only one mapping to each ref/pos will be kept)
"""
    sys.exit(1)


def bam_merge(fname, infiles, tags=['AS+', 'NM-'], discard=False, keepall=False, quiet=False):
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
        bamgens.append(ngsutils.bam.bam_batch_reads(bam))

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
        best_reads = []
        best_source = 0

        mappings = {}
        unmapped = None

        first_group = last_reads[0]
        for i in xrange(len(last_reads)):
            if not last_reads[i]:
                continue

            match = False
            for read in last_reads[i]:                        
                if read.qname == first_group[0].qname:
                    match = True

                    if not read.is_unmapped:
                        tag_val = []
                        for tag in tags:
                            val = float(read.opt(tag[:2]))
                            if tag[-1] == '-':
                                val = -val
                            tag_val.append(val)

                        if keepall:
                            if not (read.tid, read.pos) in mappings:
                                mappings[(read.tid, read.pos)] = (tag_val, i, read)
                            elif tag_val > mappings[(read.tid, read.pos)][0]:
                                mappings[(read.tid, read.pos)] = (tag_val, i, read)
                        else:        
                            if not best_val or tag_val > best_val:
                                best_val = tag_val
                                best_reads = [read]
                                best_source = i
                            elif tag_val == best_val:
                                best_reads.append(read)
                    elif not discard:
                        unmapped = read
            if match:
                last_reads[i] = None

        if keepall and mappings:
            outs = []
            for k in mappings:
                outs.append(mappings[k])

            for tagval, i, read in sorted(outs):
                counts[i] += 1
                outfile.write(read)

        elif best_reads:
            counts[best_source] += 1
            for read in best_reads:
                outfile.write(read)
        elif unmapped:
            unmapped += 1

            if not discard:
                outfile.write(unmapped)

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
    keepall = False
    tags = []

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-tag':
            tags.append(arg)
            last = None
        elif arg in ['-tag']:
            last = arg
        elif arg == '-keepall':
            keepall = True
        elif arg == '-discard':
            discard = True
        elif not outfile:
            outfile = arg
        elif os.path.exists(arg):
            infiles.append(arg)

    if not tags:
        tags = ['AS+', 'NM-']

    if not infiles or not outfile:
        usage()
    else:
        bam_merge(outfile, infiles, tags, discard, keepall)
