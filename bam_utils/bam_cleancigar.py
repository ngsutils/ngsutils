#!/usr/bin/env python
'''
Fixes BAM files where the CIGAR alignment has a zero length element
'''
import sys
import os
import pysam


def bam_cleancigar(infile, outfile):
    bam = pysam.Samfile(infile, "rb")
    out = pysam.Samfile(outfile, "wb", template=bam)
    total = 0
    count = 0
    for read in bam:
        newcigar = []
        changed = False
        if not read.is_unmapped:
            for op, length in read.cigar:
                if length > 0:
                    newcigar.append((op, length))
                else:
                    changed = True

            if changed:
                read.cigar = newcigar
                count += 1
        total += 1
        out.write(read)

    bam.close()
    out.close()
    sys.stderr.write('Wrote: %s reads\nAltered: %s\n' % (total, count))


def usage():
    print __doc__
    print """Usage: bamutils cleancigar {-f} inbamfile outbamfile

Options:
  -f     Force overwriting an existing outfile
"""
    sys.exit(-1)

if __name__ == "__main__":
    infile = None
    outfile = None
    force = False

    for arg in sys.argv[1:]:
        if arg == "-h":
            usage()
        elif arg == "-f":
            force = True
        elif not infile:
            if os.path.exists(os.path.expanduser(arg)):
                infile = os.path.expanduser(arg)
            else:
                sys.stderr.write("File: %s not found!" % arg)
                usage()
        elif not outfile:
            if force or not os.path.exists(os.path.expanduser(arg)):
                outfile = arg
            else:
                sys.stderr.write("File: %s exists! Not overwriting without -f force." % arg)
                usage()
        else:
            usage()

    if not infile or not outfile:
        usage()

    bam_cleancigar(infile, outfile)
