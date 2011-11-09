#!/usr/bin/env python
'''
Tags each read in a BAM file

Currently supported tags:
  -suffix
  -xs
'''
import sys
import os
from support.eta import ETA
import pysam

def bam_tag(infname, outfname, suffix, xs, tag):
    # write to a tmp file, then move afterwards
    tmp = os.path.join(os.path.dirname(outfname),'.tmp.%s' % os.path.basename(outfname))

    bamreader = BamReader(infname)
    outbam = pysam.Samfile(tmp,'wb', template=bamreader.bamfile)
    
    chain = bamreader
    
    if suffix:
        chain = Suffix(chain, suffix)
    
    if xs:
        chain = CufflinksXS(chain)

    if tag:
        chain = Tag(chain,tag)

    if chain == bamreader:
        sys.stderr.write('Missing tag argument\n')
        usage()        

    for read in chain.filter():
        outbam.write(read)
    
    outbam.close()

    if os.path.exists(outfname):
        os.unlink(outfname)  # Not really needed on *nix
    os.rename(tmp,outfname)

class BamReader(object):
    def __init__(self, fname):
        self.bamfile = pysam.Samfile(fname,"rb")
        
    def filter(self):
        eta = ETA(0,bamfile=self.bamfile)
        for read in self.bamfile:
            eta.print_status(extra=read.qname,bam_pos=(read.rname,read.pos))
            yield read

        eta.done()
        self.bamfile.close()

class Suffix(object):
    def __init__(self, parent, suffix):
        self.parent = parent
        self.suffix = suffix
    
    def filter(self):
        for read in self.parent.filter():
            read.qname = "%s%s" % (read.qname, suffix)
            yield read

class Tag(object):
    def __init__(self, parent, tag):
        self.parent = parent
        self.tag = tag.rsplit(':',1)
        if self.tag[0][-2:].lower() == ':i':
            self.tag = (self.tag[0],int(self.tag[1]))
        if self.tag[0][-2:].lower() == ':f':
            self.tag = (self.tag[0],float(self.tag[1]))

    def filter(self):
        for read in self.parent.filter():
            newtags = []
            for key, val in read.tags:
                newtags.append((key, val))
            
            newtags.append(self.tag)
            read.tags = newtags
            yield read


class CufflinksXS(object):
    def __init__(self, parent):
        self.parent = parent
    
    def filter(self):
        for read in self.parent.filter():
            newtags = []
            for key, val in read.tags:
                newtags.append((key, val))

            if not read.is_unmapped:
                if read.is_reverse:
                    newtags.append(('XS:A','-'))
                else:
                    newtags.append(('XS:A','+'))

            read.tags = newtags
            yield read
            
            

def usage():
    print __doc__
    print """\
Usage: bamutils tag {opts} in.bamfile out.bamfile

Arguments:
  in.bamfile    The input BAM file
  out.bamfile   The name of the new output BAM file
  
Options:
  -suffix suff  A suffix to add to each read name
  -xs           Add the XS:A tag for +/- strandedness (req'd by Cufflinks)
  -tag tag      Add an arbitrary tag (ex: -tag XX:z:test)
  -f            Force overwriting the output BAM file if it exists
"""
    sys.exit(1)

if __name__ == "__main__":
    infname = None
    outfname = None
    suffix = None
    tag = None
    force = False
    xs = False
    last = None
    
    for arg in sys.argv[1:]:
        if arg == '-f':
            force = True
        elif last == '-suffix': 
            suffix = arg
            last = None
        elif last == '-tag': 
            tag = arg
            last = None
        elif arg in ['-suffix', '-tag']:
            last = arg
        elif arg == '-xs':
            xs = True
        elif not infname and os.path.exists(arg):
            infname = arg
        elif not outfname:
            if not force and os.path.exists(arg):
                sys.stderr.write('%s already exists! Not overwriting without force (-f)!' % arg)
                sys.exit(1)
            outfname = arg
    
    if not infname or not outfname:
        usage()

    bam_tag(infname, outfname, suffix, xs, tag)

