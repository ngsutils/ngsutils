#!/usr/bin/env python
'''
Takes a RefIso gene model and a genomic FASTA file and produces a splice-
junction library in FASTA format.
'''

import os,sys
import support.refiso

def usage():
    print __doc__
    print """\
Usage: %s {opts} refiso.txt{.gz} ref.fasta

Arguments
  refiso.txt      Gene model in RefIso format
  ref.fasta       Reference genome in FASTA or RAGZ format 
                  (samtools indexed ref.fasta.fai req'd)

Options
  -frag size      Number of bases on either side of the junction to include
                  [default 46]
  -min size       Minimum size of a junction            
                  [default 50]

""" % (os.path.basename(sys.argv[0]),)
    sys.exit(1)

if __name__ == '__main__':
    model = None
    fasta = None
    frag_size = 46
    min_size = 50
    last = None
    
    for arg in sys.argv[1:]:
        if last == '-frag':
            frag_size = int(arg)
            last = None
        elif last == '-min':
            min_size = int(arg)
            last = None
        elif arg in ['-frag','-min']:
            last = arg
        elif arg == '-h':
            usage()
        elif model is None and os.path.exists(arg):
            model = arg
        elif fasta is None and os.path.exists(arg) and os.path.exists('%s.fai' % arg):
            fasta = arg
    
    if not model or not fasta:
        usage()

    support.refiso.refiso_junctions(model,fasta,frag_size,min_size)
