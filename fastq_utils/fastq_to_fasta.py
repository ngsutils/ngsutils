#!/usr/bin/env python
'''
Convert FASTQ to FASTA.  Optionally outputs just the quality values.

Same format as SOLiD csfasta / qual files.
'''

import os,sys,gzip

def read_fastq(fname):
    
    if fname[-3:].lower() == '.gz':
        f = gzip.open()
    else:
        f = open(fname)
        
    try:
        while f:
            name = f.next()
            seq = f.next()
            f.next() # plus
            qual = f.next()
            quals = ' '.join([ord(x)-33 for x in qual])
            
            yield (name[1:],seq,quals)
    except:
        pass
    f.close()

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