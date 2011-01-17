#!/usr/bin/env python
"""
This takes a BAM file that has been mapped to a genomic region and converts 
the mapping to genomic coordinates.

For example: 
chr1:1000-2000    20    50M

is converted to:
chr1    1020    50M
"""

import os,sys,gzip
sys.path.append(os.path.join(os.path.dirname(__file__),"..","utils")) #eta
sys.path.append(os.path.join(os.path.dirname(__file__),"..","ext")) #pysam

from eta import ETA
import pysam

def usage():
    base = os.path.basename(sys.argv[0])
    print __doc__
    print """
Usage: %s in.bam out.bam chrom.sizes

""" % (base,)
    sys.exit(1)

bam_cigar = ['M','I','D','N','S','H','P']

_region_cache = {}

def region_pos_to_genomic_pos(name,start,cigar):
    '''
        converts a junction position to a genomic location given a junction
        ref name, the junction position, and the cigar alignment.

        returns: (genomic ref, genomic pos, genomic cigar)
    '''
    # example name: chr3R:17630851-17630897,17634338-17634384       17       39M
    if name in _region_cache:
        chrom,fragments = _region_cache[name]
    else:
        c1 = name.split(':')
        chrom = c1[0]

        fragments = []
        for fragment in c1[1].split(','):
            s,e = fragment.split('-')
            fragments.append((int(s),int(e)))
        
        _region_cache[name] = (chrom,fragments)
    
    chr_cigar = []
    chr_start = fragments[0][0]
    
    read_start = int(start)
    
    frag_idx = 0
    frag_start = 0
    frag_end = 0
    
    for i,(s,e) in enumerate(fragments):
        if chr_start+read_start < e:
            chr_start += read_start
            frag_idx = i
            frag_start = s
            frag_end = e
            break
            
        else:
            chr_start += (e-s)
            read_start -= (e-s)
    
    cur_pos = chr_start
    
    for op,length in cigar:
        if op == 1:
            chr_cigar.append((op,length))
            
        elif op in [0,2]:
            if cur_pos + length <= frag_end:
                cur_pos += length
                chr_cigar.append((op,length))
                
            else:
                while cur_pos + length > frag_end:
                    chr_cigar.append((op,frag_end-cur_pos))
                    length -= (frag_end-cur_pos)
                    cur_pos = frag_end
                    frag_idx += 1
                    frag_start, frag_end = fragments[frag_idx]
                    chr_cigar.append((3,frag_start-cur_pos))
                    cur_pos = frag_start

                cur_pos = cur_pos + length
                chr_cigar.append((op,length))
        else:
            print "Unsupported CIGAR operation (%s)" % bam_cigar[op]
            sys.exit(1)
    
    return (chrom,chr_start,chr_cigar)

def bam_region_to_ref(infile,outfile,chrom_sizes):
    bamfile = pysam.Samfile(infile,"rb")
    header = bamfile.header
    header['SQ'] = []

    with open(chrom_sizes) as f:
        for line in f:
            if line[0] != '#':
                cols = line.strip().split('\t')
                header['SQ'].append({'LN':int(cols[1]),'SN':cols[0]})
    
    outfile = pysam.Samfile(outfile,"wb",header=header)
    
    eta = ETA(0,bamfile=bamfile)
    count = 0
    
    for read in bamfile.fetch():
        eta.print_status(bam_pos=(read.rname,read.pos))
        chrom,pos,cigar = region_pos_to_genomic_pos(bamfile.getrname(read.rname),read.pos,read.cigar)
        read.pos = pos
        read.cigar = cigar
        
        found=False
        for i,name in enumerate(outfile.references):
            if name == chrom:
                read.rname = i
                found=True
                break
        if not found:
            print "Can't find chrom: %s" % chrom
            sys.exit(1)
        outfile.write(read)

    eta.done()
    bamfile.close()
    outfile.close()

def test():
    
    chrom,start,cigar = region_pos_to_genomic_pos('chr1:1000-1050,2000-2050,3000-4000',25,[(0,100)])
    
    assert chrom == 'chr1'
    assert start == 1025
    assert cigar == [(0, 25), (3, 950), (0, 50), (3, 950), (0, 25)]
    
if __name__ == '__main__':
    infile = None
    outfile = None
    chrom_sizes = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif not infile:
            infile = arg
        elif not outfile:
            outfile = arg
        elif not chrom_sizes:
            chrom_sizes = arg

    if not infile or not outfile or not chrom_sizes:
        usage()
    else:
        bam_region_to_ref(infile,outfile,chrom_sizes)
        
