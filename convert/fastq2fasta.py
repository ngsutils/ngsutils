#!/usr/bin/env python

import sys,gzip

def fastq_to_fasta(fname):
    if fname == '-':
        f = sys.stdin
    elif fname[-3:] == '.gz':
        f = gzip.open(fname)
    else:
        f = open(fname)

    try:
        while f:
            n=f.next()
            s=f.next()
            f.next() # plus
            f.next() # qual
        
            sys.stdout.write('>%s%s' % (n[1:],s))
    except:
        pass
        
    if f != sys.stdin:
        f.close()

if __name__ == '__main__':
    fastq_to_fasta(sys.argv[1])
