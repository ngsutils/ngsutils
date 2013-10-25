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
from ngsutils.bam import bam_iter, cigar_tostr, bam_open


def bam_export(bam, mapped=True, unmapped=True, whitelist=None, blacklist=None, fields=None, out=sys.stdout, quiet=False):
    for read in bam_iter(bam, quiet=quiet):
        if whitelist and not read.qname in whitelist:
            continue
        if blacklist and read.qname in blacklist:
            continue

        try:
            if mapped and not read.is_unmapped:
                export_read(bam, read, fields, out)
            elif unmapped and read.is_unmapped:
                export_read(bam, read, fields, out)
        except IOError:
            break


def export_read(bamfile, read, fields, out=sys.stdout):
    cols = []
    for field in fields:
        if field == '-name':
            cols.append(read.qname)
        elif field == '-ref':
            if read.tid == -1:
                cols.append('*')
            elif bamfile:
                cols.append(bamfile.getrname(read.tid))
            else:
                cols.append('?')
        elif field == '-pos':
            cols.append(read.pos + 1)  # output 1-based
        elif field == '-strand':
            cols.append('-' if read.is_reverse else '+')
        elif field == '-cigar':
            if read.cigar:
                cols.append(cigar_tostr(read.cigar))
            else:
                cols.append('*')
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
            elif bamfile:
                cols.append(bamfile.getrname(read.rnext))
            else:
                cols.append('?')
        elif field == '-nextpos':
            if read.rnext == -1:
                cols.append(0)
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
                    elif len(val) == 1:
                        cols.append('%s:A:%s' % (tag, val))
                    else:
                        cols.append('%s:Z:%s' % (tag, val))
            else:
                try:
                    cols.append(read.opt(field[5:]))
                except KeyError:
                    cols.append('')

    out.write('%s\n' % '\t'.join([str(x) for x in cols]))


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
  -strand        Mapped strand (+/-)
  -cigar         CIGAR alignment string
  -flags         Mapping flags (base 10 number)
  -seq           Sequence
  -qual          Quality sequence
  -mapq          MAPQ score
  -nextref       Next mapped reference (paired-end)
  -nextpos       Next mapped position (paired-end)
  -tlen          Template length (paired-end)
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

    if not fields:
        fields = ['-name', '-ref', '-pos', '-cigar']

    if not unmapped and not mapped:
        unmapped = True
        mapped = True

    bamfile = bam_open(fname)
    bam_export(bamfile, mapped, unmapped, wl, bl, fields)
    bamfile.close()
