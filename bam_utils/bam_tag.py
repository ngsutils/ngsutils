#!/usr/bin/env python
'''
Tags each read in a BAM file with a suffix.  This enables merging multiple
BAM files and avoiding any name collision.
'''
import sys
import os
from support.eta import ETA
import pysam

def bam_tag(suffix, infname, outfname):
    bamfile = pysam.Samfile(infname,"rb")
    eta = ETA(0,bamfile=bamfile)
    
    # write to a tmp file, then move afterwards
    tmp = os.path.join(os.path.dirname(outfname),'.tmp.%s' % os.path.basename(outfname))
    outbam = pysam.Samfile(tmp,'wb')

    for read in bamfile:
        eta.print_status(extra=read.qname,bam_pos=(read.rname,read.pos))
        
        read.qname = "%s%s" % (read.qname, suffix)
        outbam.write(read)
    
    eta.done()
    bamfile.close()
    outbam.close()

    if os.path.exists(outfname):
        os.unlink(outfname)  # Not really needed on *nix
    os.rename(tmp,outfname)

def usage():
    print __doc__
    print """\
Usage: bamutils tag {opts} suffix in.bamfile out.bamfile

Arguments:
  suffix        The suffix to add to each read name
  in.bamfile    The input BAM file
  out.bamfile   The name of the new output BAM file
  
Options:
  -f            Force overwriting the output BAM file if it exists
"""
    sys.exit(1)

if __name__ == "__main__":
    infname = None
    outfname = None
    suffix = ""
    force = False
    
    for arg in sys.argv[1:]:
        if arg == '-f':
            force = True
        elif not suffix: 
            suffix = arg
        elif not infname and os.path.exists(arg):
            infname = arg
        elif not outfname:
            if not force and os.path.exists(arg):
                sys.stderr.write('%s already exists! Not overwriting without force (-f)!')
                sys.exit(1)
            outfname = arg
            
    if not infname or not suffix or not outfname:
        usage()
        sys.exit(1)

    bam_tag(suffix, infname, outfname)

