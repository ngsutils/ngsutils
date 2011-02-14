#!/usr/bin/env python
'''
Builds a RefIso format gene model from UCSC Known Gene annotations
'''

import os,sys
import support.refiso

def usage():
    print __doc__
    print """\
Usage: %s {opts} knownGene.txt{.gz} kgXref.txt{.gz} knownIsoforms.txt{.gz} 

Arguments
  knownGene       Known Gene model from UCSC
  kgXref          Known Gene cross-references from UCSC
  knownIsoforms   Known Gene isoformst from UCSC

Options
  -extend num     merge isoforms if they are {num} bases away from each other
                  [default 0]
""" % (os.path.basename(sys.argv[0]),)
    sys.exit(1)

if __name__ == '__main__':
    kg = None
    xref = None
    isof = None
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
        elif kg is None and os.path.exists(arg):
            kg = arg
        elif xref is None and os.path.exists(arg):
            xref = arg
        elif isof is None and os.path.exists(arg):
            isof = arg
    
    if not kg or not xref or not isof:
        usage()

    support.refiso.kg_to_refiso(kg,xref,isof,extend)
