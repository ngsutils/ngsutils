#!/usr/bin/env python
## category General
## desc Export reads, mapped positions, and other tags
'''
Outputs information about reads in a BAM file.

This can be used to optionally extract various fields from a BAM file on a
read-by-read basis. By default this will export the read-name, the mapped
reference, the mapped position, and the CIGAR alignment string for a read.
'''

import sys
import os
from ngsutils.support.eta import ETA
import pysam

bam_cigar = ['M', 'I', 'D', 'N', 'S', 'H', 'P']


def bam_export(fname, mapped=False, unmapped=False, whitelist=None, blacklist=None, fields=None):
    bamfile = pysam.Samfile(fname, "rb")
    eta = ETA(0, bamfile=bamfile)

    for read in bamfile:
        eta.print_status(extra=read.qname, bam_pos=(read.rname, read.pos))
        if whitelist and not read.qname in whitelist:
            continue
        if blacklist and read.qname in blacklist:
            continue

        try:
            if mapped and not read.is_unmapped:
                export_read(bamfile, read, fields)
            elif unmapped and read.is_unmapped:
                export_read(bamfile, read, fields)
        except IOError:
            break

    eta.done()
    bamfile.close()


def export_read(bamfile, read, fields):
    cols = []
    for field in fields:
        if field == '-name':
            cols.append(read.qname)
        elif field == '-ref':
            cols.append(bamfile.getrname(read.rname))
        elif field == '-pos':
            cols.append(read.pos + 1)  # output 1-based
        elif field == '-cigar':
            cols.append(''.join(['%s%s' % (length, bam_cigar[op]) for op, length in read.cigar]))
        elif field == '-flags':
            cols.append(read.flag)
        elif field == '-seq':
            cols.append(read.seq)
        elif field == '-qual':
            cols.append(read.qual)
        elif field == '-mapq':
            cols.append(read.mapq)
        elif field == '-nextref':
            if read.rnext == -1:
                cols.append('*')
            else:
                cols.append(bamfile.getrname(read.rnext))
        elif field == '-nextpos':
            if read.rnext == -1:
                cols.append(-1)
            else:
                cols.append(read.pnext + 1)  # output 1-based
        elif field == '-tlen':
            cols.append(read.tlen)
        elif field == '-isize':
            cols.append(read.isize)
        elif field[:5] == '-tag:':
            if field[5:] == '*':
                for tag, val in read.tags:
                    if type(val) == int:
                        cols.append('%s:i:%s' % (tag, val))
                    elif type(val) == float:
                        cols.append('%s:f:%s' % (tag, val))
                    else:
                        cols.append('%s:Z:%s' % (tag, val))
            else:
                cols.append(read.opt(field[5:]))

    print '\t'.join([str(x) for x in cols])


def usage():
    print __doc__
    print """\
Usage: bamutils export {opts} {fields} bamfile

Options:
  -mapped              Output only mapped reads
  -unmapped            Output only unmapped reads

  -whitelist file.txt  Output only reads that are listed in a text file
  -blacklist file.txt  Output only reads that are not listed in a text file

Fields:
  -name          Read name
  -ref           Mapped reference (chrom)
  -pos           Mapped position (1-based)
  -cigar         CIGAR alignment string
  -flags         Mapping flags (base 10 number)
  -seq           Sequence
  -qual          Quality sequence
  -mapq          MAPQ score
  -nextref       Next mapped reference (paired-end)
  -nextpos       Next mapped position (paired-end)
  -tlen          Template length (paired, observed)
  -tag:tag_name  Any tag
                 For example:
                   -tag:AS -tag:NH
                   Outputs the alignment score and the edit distance
  -tag:*         Outputs all tags in SAM format


  Default fields: -name -ref -pos -cigar
"""
    sys.exit(1)

if __name__ == "__main__":
    mapped = False
    unmapped = False
    fname = None
    wl = None
    bl = None
    last = None
    fields = []

    for arg in sys.argv[1:]:
        if last == '-whitelist':
            if not os.path.exists(arg):
                print "Error: %s missing!" % arg
                usage()
            with open(arg) as f:
                wl = [x.strip() for x in f]
            last = None
        elif last == '-blacklist':
            if not os.path.exists(arg):
                print "Error: %s missing!" % arg
                usage()
            with open(arg) as f:
                bl = [x.strip() for x in f]
            last = None
        elif arg in ['-blacklist', '-whitelist']:
            last = arg
        elif arg == '-h':
            usage()
        elif arg == '-unmapped':
            unmapped = True
        elif arg == '-mapped':
            mapped = True
        elif os.path.exists(arg):
            fname = arg
        elif arg[0] == '-':
            fields.append(arg)

    if not fname:
        usage()
        sys.exit(1)

    if not fields:
        fields = ['-name', '-ref', '-pos', '-cigar']

    if not unmapped and not mapped:
        bam_export(fname, True, True, wl, bl, fields)
    else:
        bam_export(fname, mapped, unmapped, wl, bl, fields)
