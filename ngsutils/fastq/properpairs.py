#!/usr/bin/env python
## category General
## desc Find properly paired reads (when fragments are filtered separately)
'''

When you filter paired FASTQ files separately, fragments my be filtered or
removed from either file. This program will then find all of the properly
paired reads. Valid pairs will be written to new output files that may be
optionally gzip compressed.

'''

import os
import sys
import gzip

import Queue

from ngsutils.fastq import FASTQ


def find_fastq_pairs(fq1, fq2, out1, out2, quiet=False):
    buffer1 = Queue.Queue()
    buffer2 = Queue.Queue()

    readnames1 = set()
    readnames2 = set()

    total1 = 0
    total2 = 0
    matched = 0
    removed = 0

    def callback():
        return 'ma/rm:%s/%s' % (matched, removed)

    gen1 = fq1.fetch(quiet=quiet, callback=callback)
    gen2 = fq2.fetch(quiet=True)

    while True:

        try:
            read1 = gen1.next()
            total1 += 1

        except:
            read1 = None

        try:
            read2 = gen2.next()
            total2 += 1
        except:
            read2 = None

        if read1 is None and read2 is None:
            break

        if read1 and read2 and read1.name == read2.name:
            read1.write(out1)
            read2.write(out2)
            matched += 1

        elif read1 and read1.name in readnames2:
            read1.write(out1)
            buffer1 = Queue.Queue()
            readnames1 = set()
            matched += 1

            cur = None
            while cur is None or cur.name != read1.name:
                if cur:
                    removed += 1
                cur = buffer2.get()
                readnames2.remove(cur.name)

            if cur:
                cur.write(out2)

            if read2:
                buffer2.put(read2)
                readnames2.add(read2.name)

        elif read2 and read2.name in readnames1:
            read2.write(out2)
            buffer2 = Queue.Queue()
            readnames2 = set()
            matched += 1

            cur = None
            while cur is None or cur.name != read2.name:
                if cur:
                    removed += 1
                cur = buffer1.get()
                readnames1.remove(cur.name)

            if cur:
                cur.write(out1)

            if read1:
                buffer1.put(read1)
                readnames1.add(read1.name)

        elif read1 and read2:
            buffer1.put(read1)
            buffer2.put(read2)

            readnames1.add(read1.name)
            readnames2.add(read2.name)

    return total1, total2, matched


def usage(msg=""):
    if msg:
        print '%s\n' % msg
    print """Usage: fastqutils properpairs filename1.fastq{.gz} filename2.fastq{.gz} output1 output2

Options:
  -f    Force overwriting output file (if it exists)
  -z    Output files should be gzip compressed
"""
    sys.exit(1)

if __name__ == '__main__':
    fqname1 = None
    fqname2 = None
    outname1 = None
    outname2 = None
    force = False

    gz = False

    for arg in sys.argv[1:]:
        if arg == '-z':
            gz = True
        elif arg == '-f':
            force = True
        elif not fqname1:
            if not os.path.exists(arg):
                usage("File %s doesn't exist!" % arg)
            fqname1 = arg
        elif not fqname2:
            if not os.path.exists(arg):
                usage("File %s doesn't exist!" % arg)
            fqname2 = arg
        elif not outname1:
            outname1 = arg
        elif not outname2:
            outname2 = arg

    if not fqname1 or not fqname2 or not outname1 or not outname2:
        usage()

    if not force:
        for fname in [outname1, outname2]:
            if os.path.exists(fname):
                usage("File %s exists!" % fname)

    fq1 = FASTQ(fqname1)
    fq2 = FASTQ(fqname2)

    if gz:
        out1 = gzip.open(outname1, 'w')
        out2 = gzip.open(outname2, 'w')
    else:
        out1 = open(outname1, 'w')
        out2 = open(outname2, 'w')

    total1, total2, matched = find_fastq_pairs(fq1, fq2, out1, out2)

    print "Totals: %s, %s" % (total1, total2)
    print "Proper pairs: %s" % matched

    fq1.close()
    fq2.close()
    out1.close()
    out2.close()
