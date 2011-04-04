#!/usr/bin/env python
'''
Convert FASTQ to FASTA.  Optionally outputs just the quality values.

Same format as SOLiD csfasta / qual files.
'''

import os,sys,gzip
from support.eta import ETA

def read_fastq(fname):
    if fname[-3:].lower() == '.gz':
        f = gzip.open(fname)
    else:
        f = open(fname)
    eta = ETA(os.stat(fname).st_size, fileobj = f)
    
    while f:
        eta.print_status()

        try:
            name = f.next().strip()
            seq = f.next().strip()
            f.next() # plus
            qual = f.next().strip()
        except:
            break

        quals = ' '.join([str(ord(x)-33) for x in qual])
        yield (name[1:],seq,quals)
        
    f.close()
    eta.done()

def export_seq(fname):
    for name,seq,quals in read_fastq(fname):
        sys.stdout.write('>%s\n%s\n' % (name,seq))

def export_qual(fname):
    for name,seq,quals in read_fastq(fname):
        sys.stdout.write('>%s\n%s\n' % (name,quals))

def usage():
    print """Usage: %s {-qual} filename.fastq{.gz}
Options:
  -qual    Export the quality values (space separated numbers)
""" % os.path.basename(sys.argv[0])
    sys.exit(1)

if __name__ == '__main__':
    qual = False
    fname = None
    
    for arg in sys.argv[1:]:
        if arg == '-qual':
            qual = True
        elif os.path.exists(arg):
            fname = arg
    
    if not fname:
        usage()
        
    if qual:
        export_qual(fname)
    else:
        export_seq(fname)