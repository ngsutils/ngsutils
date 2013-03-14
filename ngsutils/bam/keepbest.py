#!/usr/bin/env python
## category General
## desc Parses BAM file and keeps the best mapping for reads that have multiple mappings
'''
Parses BAM file and keeps the best mapping for reads that have multiple
mappings. If an aligner outputs multiple mappings for each read, this will
remove mappings, keeping only the best one.
'''
import sys
import os

import ngsutils.bam


def get_read_tag_value(read, tag):
    if tag == 'MAPQ':
        return read.mapq

    for name, value in read.tags:
        if name == tag:
            return value

    return -1


def _output_best(lastreads, outfile):
    if not lastreads:
        return

    lastreads.sort()

    best = lastreads[-1][0]
    for score, read in lastreads[::-1]:
        if score == best:
            outfile.write(read)
        else:
            break


def bam_keepbest(fname, outname, tag="AS"):
    bamfile = ngsutils.bam.bam_open(fname)
    outfile = ngsutils.bam.bam_open(outname, "w", template=bamfile)

    lastreads = []
    for read in ngsutils.bam.bam_iter(bamfile):
        if lastreads and read.qname == lastreads[0][1].qname:
            lastreads.append((get_read_tag_value(read, tag), read))
        else:
            if lastreads:
                _output_best(lastreads, outfile)

            lastreads = [(get_read_tag_value(read, tag), read), ]

    _output_best(lastreads, outfile)

    bamfile.close()
    outfile.close()


def usage():
    print __doc__
    print """Usage: bamutils keepbest {-tag tag} infile.bam outfile.bam

Options:
   -tag tag    Use {tag} to determine which mappings to keep. This can be any
               tag present in the BAM file or "MAPQ" (default: AS).

"""
    sys.exit(-1)

if __name__ == "__main__":
    fname = None
    outname = None
    tag = "AS"
    last = None

    for arg in sys.argv[1:]:
        if last == '-tag':
            tag = arg
            last = None
        elif arg in ['-tag']:
            last = arg
        elif arg == "-h":
            usage()
        elif not fname and os.path.exists(arg):
            fname = arg
        elif not outname:
            if os.path.exists(arg):
                sys.stderr.write("Output file: %s exists!" % arg)
                usage()

            outname = arg
        else:
            sys.stderr.write("Unknown option: %s!\n" % arg)
            usage()

    if not fname:
        usage()

    bam_keepbest(fname, outname, tag)
