#!/usr/bin/env python
'''
Takes two paired FASTQ files and merges them together.  The read names need to 
be the same and in the same order.  If the reads are from an Illumina machine, 
the pairs will be named the same, except for their last character.  This is 
allowed.

If either reads contains the phrase "QCFAIL" in the description both will be 
removed from the final FASTQ file.

Note: Illumina paired end reads are: Read1 (forward) Read2 (reverse).
'''

import os,sys,gzip

try:
    sys.path.append(os.path.normpath(os.path.join(os.path.dirname(__file__),'..')))
    from support.eta import ETA
except:
    ETA = None

_comps = { 'A':'T','T':'A','G':'C','C':'G' }

def _comp(n):
    if n in _comps:
        return _comps[n]
    return 'N'

def read_fastq(fileobj,rev=False,colorspace=False):
    try:
        name = fileobj.next().strip()
        seq = fileobj.next().strip()
        plus = fileobj.next()
        qual = fileobj.next().strip()
        if not rev:
            return (name,seq,qual)
        else:
            if colorspace:
                # colorspace just needs to be reversed, not complimented
                return (name,seq[::-1],qual[::-1])
            return (name,''.join([ _comp(x) for x in seq[::-1].upper()]),qual[::-1])
            
    except Exception,e:
        print e
        return None

def merge_fastq(filename1, filename2, rev_one=False, rev_two=False, colorspace=False):
    teller = None
    if filename1[-3:] == '.gz':
        f1 = gzip.open(filename1)
    else:
        f1 = open(filename1)
    
    if filename2[-3:] == '.gz':
        f2 = gzip.open(filename2)
    else:
        f2 = open(filename2)
    
    read1 = read_fastq(f1,rev_one,colorspace)
    read2 = read_fastq(f2,rev_two,colorspace)

    if ETA:
        eta = ETA(os.stat(filename1).st_size, fileobj=f1)
    else:
        eta = None

    while read1 and read2:
        if eta:
            eta.print_status()
        n1,s1,q1 = read1
        n2,s2,q2 = read2

        if ' ' in n1:
            k1,desc1 = n1.split(' ',1)
        else:
            k1 = n1
            desc1=''

        if ' ' in n2:
            k2,desc2 = n2.split(' ',1)
        else:
            k2 = n2
            desc2=''
        

        if 'QCFAIIL' not in desc1 and 'QCFAIL' not in desc2:
                if k1 == k2:
                    print '%s\n%s\n+\n%s' % (n1,s1,q1)
                    print '%s\n%s\n+\n%s' % (n2,s2,q2)
                elif k1[:-1] == k2[:-1]:
                    print '%s %s\n%s\n+\n%s' % (k1[:-1],k1[-1],s1,q1)
                    print '%s %s\n%s\n+\n%s' % (k2[:-1],k2[-1],s2,q2)
                 
        # get next reads...
        read1 = read_fastq(f1,rev_one,colorspace)
        read2 = read_fastq(f2,rev_two,colorspace)
        
    if eta:
        eta.done()

    f1.close()
    f2.close()

def usage():
    print __doc__
    print "Usage: %s {-cs} {-rev} file1 {-rev} file2" % os.path.basename(sys.argv[0])
    print ""
    sys.exit(-1)

if __name__ == '__main__':
    rev_one = False
    first = None
    rev_two = False
    second = None
    colorspace = False
    
    for arg in sys.argv[1:]:
        if arg == '-cs':
            colorspace = True
        elif arg == '-rev' and not first:
            rev_one = True
        elif arg == '-rev' and not second:
            rev_two = True
        elif not first and os.path.exists(arg):
            first = arg
        elif not second and os.path.exists(arg):
            second = arg
        else:
            usage()

    if not first or not second:
        usage()

    merge_fastq(first,second,rev_one,rev_two,colorspace)
