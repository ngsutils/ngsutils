#!/usr/bin/env python
'''
Splits a FASTQ file into multiple smaller files

Output is a set of gzip compressed FASTQ files
'''

import os,sys,gzip

from fastq_utils import read_fastq,is_paired_fastq

def fastq_split(fname,outbase,chunks,ignore_pairs = False):
    i=0
    chunk=1
    if ignore_pairs:
        is_paired = False
    else:
        is_paired = is_paired_fastq(fname)

    outs = []
    for i in xrange(chunks):
        outs.append(gzip.open('%s.%s.fastq.gz' % (outbase,i+1),'w'))

    i = 0
    last_name = None

    for name,seq,qual in read_fastq(fname):
        sn = name.split()[0]
        if not is_paired:
            i += 1
        elif sn != last_name:
            i += 1
        
        if i >= len(outs):
            i = 0
            
        last_name = sn

        outs[i].write('%s\n%s\n+\n%s\n' % (name,seq,qual))
                
    for out in outs:
        out.close()

def usage():
    print __doc__
    print """\
Usage: fastqutils split {opts} filename.fastq{.gz} out_template num_chunks

Options:
  -ignorepaired    Normally for paired-end samples, each read of the pair is 
                   written to the same file. With this option, paired reads 
                   can be written to any file. This is useful for splitting a
                   paired FASTQ file back into separate files for each 
                   fragment.

"""
    sys.exit(1)

if __name__ == '__main__':
    fname = None
    outtemplate = None
    chunks = 0
    ignore_pairs = False

    for arg in sys.argv[1:]:
        if arg == '-ignorepaired':
            ignore_pairs = True
        elif not fname and os.path.exists(arg):
            fname = arg
        elif not outtemplate:
            outtemplate = arg
        else:
            chunks = int(arg)

    if not fname or not chunks or not outtemplate:
        usage()
        
    fastq_split(fname,outtemplate,chunks,ignore_pairs)
