#!/usr/bin/env python
"""
Given a BAM file, this script will split it into smaller BAM files with a 
limit on the number of reads included.
"""

import os,sys,gzip
from support.eta import ETA
import pysam

def usage():
    base = os.path.basename(sys.argv[0])
    print __doc__
    print """
Usage: %s {-n num} in.bam out_template_name

out_template_name will be the template for the smaller BAM files.  They will 
be named "out_template_name.N.bam" where out_template_name is the given 
argument and N is the file number.

Options:
    -n      The number of reads to include in sub-files
            (default: 1000000)
""" % (base,)
    sys.exit(1)

def bam_split(infile,out_template,read_count = 1000000):
    bamfile = pysam.Samfile(infile,"rb")
    outfile = None
    
    file_count = 0
    
    eta = ETA(0,bamfile=bamfile)
    count = 0
    fname = ""
    for read in bamfile.fetch():
        if not outfile or count >= read_count:
            if outfile:
                outfile.close()
            file_count += 1
            count = 0
            fname = '%s.%s.bam' % (out_template,file_count)
            outfile = pysam.Samfile(fname,"wb",template=bamfile)
            
        eta.print_status(extra=os.path.basename(fname),bam_pos=(read.rname,read.pos))
        outfile.write(read)
        count += 1

    eta.done()
    bamfile.close()
    outfile.close()
    sys.stderr.write("Split into %s files" % (file_count))


if __name__ == '__main__':
    infile = None
    outfile = None
    num=1000000
    last = None

    for arg in sys.argv[1:]:
        if last == '-n':
            num = int(arg)
            last = None
        elif arg == '-h':
                usage()
        elif arg in ['-n',]:
            last = arg
        elif not infile:
            infile = arg
        elif not outfile:
            outfile = arg

    if not infile or not outfile:
        usage()
    else:
        bam_split(infile,outfile,num)
        
