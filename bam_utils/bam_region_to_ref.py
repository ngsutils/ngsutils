#!/usr/bin/env python
"""
This takes a BAM file that has been mapped to a genomic region and converts 
the mapping to genomic coordinates.  This can be used to convert reads mapped
against a junction library or targetted resequencing back to genomic 
coordinates.  The names of the reference sequences should be named:
chrom:start-end (0-based start).  If there are gaps (junctions), they should
be named: chrom:start-end,start-end,etc...

In this is the case, this script will ensure the proper conversion of 
reference, start position, and CIGAR alignment.

Example 1: 
chr1:1000-2000    20    50M

converted to:
chr1    1020    50M

Example 2:
chr1:1000-1050,2000-2050,3000-4000    25    100M

converted to:
chr1    1025    25M950N50M950M25M

"""

import os,sys,gzip
from support.eta import ETA
import pysam

def usage():
    base = os.path.basename(sys.argv[0])
    print __doc__
    print """
Usage: %s {-overlap} in.bam out.bam chrom.sizes

Options:
-overlap    Require that all reads must overlap a splice junction
            by 4 bases. (Also removes unmapped reads)

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

def is_junction_valid(cigar,min_overlap=4):
    '''
        Does the genomic cigar alignment represent a 'good' alignment.
        Used for checking junction->genome alignments

        1) the alignment must not start at a splice junction
        2) the alignment must not start or end with an overhang
        3) the alignment must overhang the splice junction by min_overlap (4)



        |-----------------|oooooooooooooooo|------------------|
                                            XXXXXXXXXXXXXXXXXXXXXXXX (bad 1)
      XXXXXXXXXXXXX (bad 2)                           XXXXXXXXXXXXXXXX (bad 2)                         
                        XX-----------------XXXXXXXXXXXXXXXXX (bad 3)

    '''
    first = True
    pre_gap = True

    pre_gap_count = 0
    post_gap_count = 0
    
    has_gap = False
    
    for op,length in cigar:
        # mapping can't start at a gap
        if first and op == 3:
            return (False,'Starts at gap (%s)'% cigar)
        first = False

        if op == 3:
            pre_gap = False
            post_gap_count = 0
            has_gap = True

        elif pre_gap:
            pre_gap_count += length
        else:
            post_gap_count += length

        # mapping must start with more than min_overlap base match

    if not has_gap:
        return (False,"Doesn't cover junction")
    elif pre_gap_count < min_overlap:
        return (False, "Too short overlap at 5' (%s)" % cigar)
    elif post_gap_count < min_overlap:
        return (False, "Too short overlap at 3' (%s)" % cigar)

    return True,''

def bam_batch_reads(bam):
    reads = []
    last = None
    for read in bam:
        if last and read.qname != last:
            yield reads
            reads = []
        last = read.qname
        reads.append(read)

    if reads:
        yield reads

def bam_region_to_ref(infile,outfile,chrom_sizes, enforce_overlap=False):
    bamfile = pysam.Samfile(infile,"rb")
    header = bamfile.header
    header['SQ'] = []

    with open(chrom_sizes) as f:
        for line in f:
            if line[0] != '#':
                cols = line.strip().split('\t')
                header['SQ'].append({'LN':int(cols[1]),'SN':cols[0]})
    
    outfile = pysam.Samfile(outfile,"wb",header=header)
    
    # eta = ETA(0,bamfile=bamfile)
    count = 0
    converted_count = 0
    invalid_count = 0
    unmapped_count = 0

    for batch in bam_batch_reads(bamfile):
        count += 1
        outreads = []

        if count > 100:
            count = 0
            # if enforce_overlap:
            #     eta.print_status(extra="conv:%d inv:%d un:%d" % (converted_count,invalid_count,unmapped_count),bam_pos=(batch[0].rname,batch[0].pos))
            # else:
            #     eta.print_status(bam_pos=(batch[0].rname,batch[0].pos))
            
        for read in batch:
            if read.is_unmapped and not read.is_secondary:
                unmapped_count+=1
                if not enforce_overlap:
                    outfile.write(read)
                continue
                
            chrom,pos,cigar = region_pos_to_genomic_pos(bamfile.getrname(read.rname),read.pos,read.cigar)
        
            read.pos = pos
            read.cigar = cigar
            
            chrom_found=False
            for i,name in enumerate(outfile.references):
                if name == chrom:
                    read.rname = i
                    chrom_found=True
                    break
            if not chrom_found:
                print "Can't find chrom: %s" % chrom
                sys.exit(1)

            if not enforce_overlap:
                outfile.write(read)
                continue

            valid,reason = is_junction_valid(cigar)
            if valid:
                converted_count += 1
                outreads.append(read)
            else:
                invalid_count += 1
        
        if enforce_overlap and outreads:
            for i,read in enumerate(outreads):
                newtags = []
                for key,val in read.tags:
                    if key == 'HI':
                        newtags.append(('HI',i+1))
                    elif key == 'IH':
                        newtags.append(('IH',len(outreads)))
                    else:
                        newtags.append((key,val))
                read.tags = newtags
                outfile.write(read)
        # else:
        #     read = pysam.AlignedRead()
        #     read.is_unmapped = True
        #     read.qname = batch[0].qname
        #     read.rname = -1
        #     read.mrnm = -1
        #     read.mpos = -1
        #     read.pos = -1
        #     read.mapq = 0
        #     read.isize = 0
        #     read.seq = batch[0].seq
        #     read.qual = batch[0].qual
        #     if batch[0].opt('CS'):
        #         # for some reason, pysam doesn't like it if you don't set them all at once
        #         read.tags = [('PG',batch[0].opt('PG')),('AS',2147483649),('NH',0),("IH",1),("HI",1),('CS',batch[0].opt('CS')),('CQ',batch[0].opt('CQ'))]
        #     else:
        #         read.tags = [('PG',batch[0].opt('PG')),('AS',2147483649),('NH',0),("IH",1),("HI",1),]
        # 
        #     outfile.write(read)
            
    # eta.done()
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
    overlap = False
    
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif arg == '-overlap':
            overlap = True
        elif not infile:
            infile = arg
        elif not outfile:
            outfile = arg
        elif not chrom_sizes:
            chrom_sizes = arg

    if not infile or not outfile or not chrom_sizes:
        usage()
    else:
        # import cProfile
        # pf = 'prof'
        # i=0
        # while os.path.exists(pf):
        #     pf = 'prof_%s' % i
        #     i+=1
        # cProfile.run('bam_region_to_ref(infile,outfile,chrom_sizes,overlap)',pf)
        bam_region_to_ref(infile,outfile,chrom_sizes,overlap)
        
