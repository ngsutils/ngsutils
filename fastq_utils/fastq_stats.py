#!/usr/bin/env python
'''
Calculates summary statistics for a FASTQ file.

It can calculate distribution of read lengths, average quality by base 
position, if a file is in colorspace, if it contains paired end data, and 
what encoding is used for the quality values (Sanger or Illumina).
'''

import os,sys

from fastq_utils import read_fastq,is_colorspace_fastq,is_paired_fastq,fastq_qualtype

def fastq_stats(fname):
    cs = is_colorspace_fastq(fname)
    if cs:
        print "Space:\tcolorspace"
    else:
        print "Space:\tnucleotide"

    pairs = is_paired_fastq(fname)
    if pairs > 0:
        print "Pairing:\tPaired-end (%s)" % pairs
    else:
        print "Pairing:\tFragmented"
    
    qual_totals = fastq_qualtype(fname)
    
    print "Quality scale:\t%s" % qual_totals[-1][1]
    print ' '.join(['(%s,%s)' % (q[1],q[0]) for q in qual_totals])
    
    lengths = []
    qualities = []
    total = 0
    line = 0
    try:
        for name,seq,qual in read_fastq(fname):
            if not name[0] == '@':
                print 'Invalid formatted record [line %s]' % line
                break

            if cs:
                if len(seq) != len(qual)+1:
                    print 'Seq / qual mismatch [line %s]' % line
                    break
            else:
                if len(seq) != len(qual):
                    print 'Seq / qual mismatch [line %s]' % line
                    break
            
            line += 4
            total += 1
            while len(qual) >= len(lengths):
                lengths.append(0)
                qualities.append(0)
        
            lengths[len(qual)] += 1
        
            for idx,q in enumerate([ord(x) for x in qual]):
                qualities[idx]+=q
    except KeyboardInterrupt:
        pass
        
    print "Number of reads:\t%s" % total
    print ""
    print "Lengths"
    for idx,count in enumerate(lengths[::-1]):
        if count:
            print "%s\t%s" % (len(lengths)-idx-1,count)

    print ""
    print "Quality by position"
    print ''.join([chr(x/total) for x in qualities])
    

def usage():
    print __doc__
    print "Usage: fastqutils stats filename.fastq{.gz}"
    sys.exit(1)

if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if os.path.exists(arg):
            fname = arg

    if not fname:
        usage()
        
    fastq_stats(fname)
        
