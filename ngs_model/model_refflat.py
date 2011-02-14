#!/usr/bin/env python
'''
Builds a RefIso format gene model from RefFlat formatted file (RefSeq)
'''

import os,sys
import support.refiso

def usage():
    print __doc__
    print """\
Usage: %s {opts} refflat.txt{.gz} 

Options
  -extend num     merge isoforms if they are {num} bases away from each other
                  [default 0]
""" % (os.path.basename(sys.argv[0]),)
    sys.exit(1)

if __name__ == '__main__':
    refflat = None
    extend = 0
    last = None
    
    for arg in sys.argv[1:]:
        if last == '-extend':
            extend = int(arg)
            last = None
        elif arg in ['-frag','-extend']:
            last = arg
        elif arg == '-h':
            usage()
        elif refflat is None and os.path.exists(arg):
            refflat = arg
    
    if not refflat:
        usage()

    support.refiso.refflat_to_refiso(refflat,extend)
