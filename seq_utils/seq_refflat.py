#!/usr/bin/env python
'''
Builds a RefIso format gene model from RefFlat formatted file (RefSeq)
'''

import os,sys,urllib
import support.refiso
from support.eta import ETA

class ETAHook(object):
    def __init__(self):
        self.eta = None
        pass
    def __call__(self,block_count,block_size,total_size):
        if self.eta is None:
            self.eta = ETA(total_size)
        self.eta.print_status(block_count * block_size)
    
    def done(self):
        self.eta.done()
        


def usage():
    print __doc__
    print """\
Usage: %s {opts} refflat.txt{.gz} 

Options
  -org    org     Download refFlat file from UCSC for org (hg18, mm9, etc)
  -extend num     merge isoforms if they are {num} bases away from each other
                  [default 0]
""" % (os.path.basename(sys.argv[0]),)
    sys.exit(1)

if __name__ == '__main__':
    refflat = None
    org = None
    extend = 0
    last = None
    
    for arg in sys.argv[1:]:
        if last == '-extend':
            extend = int(arg)
            last = None
        elif last == '-org':
            org = arg
        elif arg in ['-org','-extend']:
            last = arg
        elif arg == '-h':
            usage()
        elif refflat is None and os.path.exists(arg):
            refflat = arg
    
    if org:
        path = os.path.join(os.path.expanduser('~'),'.annotations',org)
        try:
            os.makedirs(path)
        except:
            pass

        refflat = os.path.join(path,'refFlat.txt.gz')
        url = 'http://hgdownload.cse.ucsc.edu/goldenPath/%s/database/refFlat.txt.gz' % org
        
        if not os.path.exists(refflat):
            hook = ETAHook()
            sys.stderr.write('%s\n' % url)
            urllib.urlretrieve(url, '%s.tmp' % refflat, hook)
            hook.done()
            os.rename('%s.tmp' % refflat, refflat)

    if not refflat:
        usage()

    support.refiso.refflat_to_refiso(refflat,extend)
