#!/usr/bin/env python
'''
Builds a RefIso format gene model from UCSC Known Gene annotations
'''

import os,sys,urllib
import support.refiso
from seq_refflat import ETAHook

def usage():
    print __doc__
    print """\
Usage: %s {opts} knownGene.txt{.gz} kgXref.txt{.gz} knownIsoforms.txt{.gz} 

Arguments
  knownGene       Known Gene model from UCSC
  kgXref          Known Gene cross-references from UCSC
  knownIsoforms   Known Gene isoformst from UCSC

Options
  -org    org     Download Known Gene files from UCSC for org (hg18, mm9, etc)
  -extend num     merge isoforms if they are {num} bases away from each other
                  [default 0]
""" % (os.path.basename(sys.argv[0]),)
    sys.exit(1)

if __name__ == '__main__':
    kg = None
    xref = None
    isof = None
    org = None
    extend = 0
    last = None
    
    for arg in sys.argv[1:]:
        if last == '-extend':
            extend = int(arg)
            last = None
        elif last == '-org':
            org = arg
        elif arg in ['-extend','-org']:
            last = arg
        elif arg == '-h':
            usage()
        elif kg is None and os.path.exists(arg):
            kg = arg
        elif xref is None and os.path.exists(arg):
            xref = arg
        elif isof is None and os.path.exists(arg):
            isof = arg

    if org:
        path = os.path.join(os.path.expanduser('~'),'.annotations',org)
        try:
            os.makedirs(path)
        except:
            pass
            
        kg = os.path.join(path,'knownGene.txt.gz')
        xref = os.path.join(path,'kgXref.txt.gz')
        isof = os.path.join(path,'knownIsoforms.txt.gz')
        
        if not os.path.exists(kg):
            hook = ETAHook()
            url = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/knownGene.txt.gz' % org
            sys.stderr.write('%s\n' % url)
            urllib.urlretrieve(url, kg, hook)
            hook.done()

        if not os.path.exists(xref):
            hook = ETAHook()
            url = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/kgXref.txt.gz' % org
            sys.stderr.write('%s\n' % url)
            urllib.urlretrieve(url, xref, hook)
            hook.done()

        if not os.path.exists(isof):
            hook = ETAHook()
            url = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/knownIsoforms.txt.gz' % org
            sys.stderr.write('%s\n' % url)
            urllib.urlretrieve(url, isof, hook)
            hook.done()
    
    if not kg or not xref or not isof:
        usage()

    support.refiso.kg_to_refiso(kg,xref,isof,extend)
