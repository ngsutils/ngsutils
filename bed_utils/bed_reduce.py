#!/usr/bin/env python
'''
Reduces BED regions to overlap them. The BED file *must* be sorted in order 
to merge them.
'''

import os,sys

def usage():
    print __doc__
    print """\
Usage: %s {-nostrand} {-c} {-extend num[,num]} bedfile

-extend 
Can either a single number (extend the same in both direction) or a comma 
delimited pair of numbers where the first number extends the region in the 5' 
direction and the second number extends the region in the 3' direction.

-c
Output the number of regions merged as the score (count).  Otherwise, it adds
the scores for all of the combined regions together.

-nostrand
Ignore strand information when merging regions
""" % os.path.basename(sys.argv[0])

def bed_reduce(fname,extend=[0,0], stranded = True, count = False):
    if fname == '-':
        f = sys.stdin
    else:
        f = open(fname)
    last = None
    
    last_chrom = None
    last_start = None
    last_end = None
    last_name = None
    last_score = 0
    last_strand = None
    
    # these are just for checking that the file is sorted
    lchrom = None
    lstart = None
    lend = None
    
    for line in f:
        chrom,start,end,name,score,strand = line.strip().split('\t')
        start = int(start)
        end = int(end)
        if lchrom == chrom:
            if start < lstart or (start == lstart and end < lend):
                sys.stderr.write('BED file is not sorted!\n')
                sys.stderr.write('chrom: %s\t%s (= %s)\n' % (lchrom,chrom,(chrom==lchrom)))
                sys.stderr.write('start: %s\t%s (< %s)\n' % (lstart,start, (lstart < start)))
                sys.stderr.write('end: %s\t%s\n' % (lend,end))
                sys.exit(1)
        
        lchrom = chrom
        lstart = start
        lend = end

        if strand == '+':
            start = start - extend[0]
            end = end + extend[1]
            if start < 0 :
                start = 0
        else:
            start = start - extend[1]
            end = end + extend[0]
            if start < 0 :
                start = 0
        
        score = int(score)


        if last_chrom == chrom and start < last_end and (not stranded or last_strand == strand):
            last_end = end
            if count:
                last_score += 1
            else:
                last_score += score
                
            last_name.append(name)
        else:
            if last_chrom:
                sys.stdout.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (last_chrom,last_start,last_end,','.join(last_name),last_score,last_strand))
            last_chrom = chrom
            last_start = start
            last_end = end
            last_name = [name,]

            if count:
                last_score = 1
            else:
                last_score = score
                
            last_strand = strand
    
    sys.stdout.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (last_chrom,last_start,last_end,','.join(last_name),last_score,last_strand))

    if f != sys.stdin:
        f.close()

if __name__ == '__main__':
    fname = None
    extend = [0,0]
    stranded = True
    count=False
    last = None
    
    for arg in sys.argv[1:]:
        if last == '-extend':
            if ',' in arg:
                extend = [int(x) for x in arg.split(',')]
            else:
                extend = [int(arg),] * 2
            last = None
        elif arg in ['-extend']:
            last = arg
        elif arg == '-nostrand':
            stranded = False
        elif arg == '-c':
            count = True
        elif not fname and (arg == '-' or os.path.exists(arg)):
            fname = arg
        else:
            print "Unknown option: %s" % arg
            usage()
            sys.exit(1)
    
    if not fname:
        usage()
        sys.exit(1)
    
    bed_reduce(fname,extend,stranded,count)