#!/usr/bin/env python
## category General
## desc Postprocesses a BAM file to rename pairs that have an extra /N value
'''
Postprocesses a BAM file to rename pairs that have an extra /N value

Some aligners output paired end reads with names ending in /1 or /2 to signify
where the read came from in a paired end experiment. This can cause problems
with downstream analysis packages that expect paired end reads to have the
exact same name.
'''
import sys
import os
import pysam


def bam_renamepair(infile, outfile, delim='/'):
    bam = pysam.Samfile(infile, "rb")
    out = pysam.Samfile(outfile, "wb", template=bam)
    for read in bam:
        name, num = read.qname.rsplit(delim, 1)
        newtags = list(read.tags)
        newtags.append('ZN', num)
        read.tags = newtags
        read.qname = name
        out.write(read)
    bam.close()
    out.close()


def usage():
    print __doc__
    print """Usage: bamutils renamepair {opts} inbamfile outbamfile

Options:
  -f           Force overwriting an existing outfile
  -delim val   The trailing delimiter to use (default '/')
"""
    sys.exit(-1)

if __name__ == "__main__":
    infile = None
    outfile = None
    delim = '/'
    last = None
    force = False

    for arg in sys.argv[1:]:
        if last == '-delim':
            delim = arg
            last = None
        elif arg == "-h":
            usage()
        elif arg == "-delim":
            last = arg
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
                sys.stderr.write(
                    "File: %s exists! Not overwriting without -f force." % arg
                    )
                usage()
        else:
            usage()

    if not infile or not outfile:
        usage()

    bam_renamepair(infile, outfile, delim)
