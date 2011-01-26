#!/usr/bin/env python
"""
Given a BAM file, this script will only allow reads that meet filtering 
criteria to be written to output. The output is another BAM file with the 
reads not matching the criteria removed.

Currently, the available filters are:
    mapped
    
    mismatch num            # mismatches or indels
                            indel always counts as 1 regardless of length
                            (requires NM tag in reads)
    mismatch_ref num ref.fa # mismatches or indel - looks up mismatches 
                            directly in a ref FASTA file (if NM not available)
    
    eq  tag_name value
    lt  tag_name value
    lte tag_name value
    gt  tag_name value
    gte tag_name value

    Where tag_name should be the full name, plus the type eg: AS:i
    
Common tags to filter by:
    AS:i    Alignment score
    NM:i    Edit distance (each indel counts as many as it's length)
    IH:i    Number of alignments
    
"""

import os,sys,gzip
from support.eta import ETA
import pysam

def usage():
    base = os.path.basename(sys.argv[0])
    print __doc__
    print """
Usage: %s in.bam out.bam {-failed out.txt} criteria...

If given, -failed, will be a text file containing the read names of all reads 
that were removed with filtering.

Example: %s filename.bam output.bam -mapped -gte AS:i 1000
This will remove all unmapped reads, as well as any reads that have an AS:i 
value less than 1000.
""" % (base,base)
    sys.exit(1)

bam_cigar = ['M','I','D','N','S','H','P']

class Mismatch(object):
    def __init__(self,num):
        self.num = int(num)
    def filter(self,bam,read):
        if read.is_unmapped:
            return False
        
        inserts = 0
        deletions = 0
        indels = 0
        edits = int(read.opt('NM'))
        #
        # NM counts the length of indels
        # We really just care about *if* there is an indel, not the size
        #
        
        for op,length in read.cigar:
            if op == 1:
                inserts += length
                indels += 1
            elif op == 2:
                deletions += length
                indels += 1
        
        mismatches = edits - inserts - deletions + indels
        
        if mismatches > self.num:
            return False
                
        return True
    def __repr__(self):
        return '%s mismatch%s in %s' % (self.num,'' if self.num == 1 else 'es',os.path.basename(self.refname))

class MismatchRef(object):
    def __init__(self,num,refname):
        self.num = int(num)
        self.refname = refname
        
        if not os.path.exists('%s.fai' % refname):
            pysam.faidx(refname)
        
        self.ref = pysam.Fastafile(refname)

    def filter(self,bam,read):
        if read.is_unmapped:
            return False
        chrom = bam.getrname(read.rname)
        start = read.pos
        end = read.aend
        
        edits = 0
        ref_pos = 0
        read_pos = 0

        for op,length in read.cigar:
            if op == 1:
                edits += 1
                read_pos += length
            elif op == 2:
                edits += 1
                ref_pos += length
            elif op == 3:
                ref_pos += length
            elif op == 0:
                refseq = self.ref.fetch(chrom,start+ref_pos,start+ref_pos+length)
                for refbase,readbase in zip(refseq,read.seq[read_pos:read_pos+length]):
                    if refbase.upper() != readbase.upper():
                        edits += 1
                ref_pos += length
                read_pos += length
                
            if edits > self.num:
                return False
        return True
    def __repr__(self):
        return '%s mismatch%s in %s' % (self.num,'' if self.num == 1 else 'es',os.path.basename(self.refname))
        

class Mapped(object):
    def __init__(self):
        pass
    def filter(self,bam,read):
        if read.is_paired and (read.is_unmapped or read.mate_is_unmapped):
            return False
        elif read.is_unmapped:
            return False
        return True
    def __repr__(self):
        return 'is mapped'

class _TagCompare(object):
    def __init__(self,tag,value):
        self.args = '%s %s' % (tag,value)
        if ':' in tag:
            self.tag = tag[0:tag.find(':')]
        else:
            self.tag = tag
            
        if tag[-1] == 'i':
            self.value = int(value)
        elif tag[-1] == 'f':
            self.value = float(value)
        else:
            self.value = value
    def __repr__(self):
        return "%s %s %s" % (self.tag,self.__class__.op,self.value)
            
class TagLessThan(_TagCompare):
    op = '<'
    def filter(self,bam,read):
        for name,value in read.tags:
            if name == self.tag and value < self.value:
                return True
        return False

class TagLessThanEquals(_TagCompare):
    op = '<='
    def filter(self,bam,read):
        for name,value in read.tags:
            if name == self.tag and value <= self.value:
                return True
        return False

class TagGreaterThan(_TagCompare):
    op = '>'
    def filter(self,bam,read):
        for name,value in read.tags:
            if name == self.tag and value > self.value:
                return True
        return False

class TagGreaterThanEquals(_TagCompare):
    op = '>='
    def filter(self,bam,read):
        for name,value in read.tags:
            if name == self.tag and value > self.value:
                return True
        return False
        
class TagEquals(_TagCompare):
    op = '='
    def filter(self,bam,read):
        for name,value in read.tags:
            if name == self.tag and value == self.value:
                return True
        return False

_criteria = {
    'mapped': Mapped,
    'lt': TagLessThan,
    'gt': TagGreaterThan,
    'lte': TagLessThanEquals,
    'gte': TagGreaterThanEquals,
    'eq': TagEquals,
    'mismatch': Mismatch,
    'mismatch_ref': MismatchRef
}

def bam_filter(infile,outfile,criteria,failedfile = None, verbose = False):
    if verbose:
        sys.stderr.write('Input file  : %s\n' % infile)
        sys.stderr.write('Output file : %s\n' % outfile)
        if failedfile:
            sys.stderr.write('Failed reads: %s\n' % failedfile)
        sys.stderr.write('Criteria:\n')
        for criterion in criteria:
            sys.stderr.write('    %s\n' % criterion)
    
        sys.stderr.write('\n')
    
    bamfile = pysam.Samfile(infile,"rb")
    outfile = pysam.Samfile(outfile,"wb",template=bamfile)
    if failedfile:
        failed_out = open(failedfile,'w')
    else:
        failed_out = None
    eta = ETA(0,bamfile=bamfile)
    
    passed = 0
    failed = 0
    
    for read in bamfile:
        eta.print_status(extra="kept:%s, failed:%s" % (passed,failed),bam_pos=(read.rname,read.pos))
        p=True

        for criterion in criteria:
            if not criterion.filter(bamfile,read):
                p = False
                failed+=1
                if failed_out:
                    failed_out.write('%s\n'%read.qname)
                #outfile.write(read_to_unmapped(read))
                break
        if p:
            passed += 1
            outfile.write(read)

    eta.done()
    bamfile.close()
    outfile.close()
    if failed_out:
        failed_out.close()
    sys.stderr.write("%s kept\n%s failed\n" % (passed,failed))

def read_to_unmapped(read):
    '''
    Example unmapped read
    2358_2045_1839    |  4 | *                             |  0    |  0   | *                  | *  |  0 |  0 | GGACGAGAAGGAGTATTTCTCCGAGAACACATTCACGGAGAGTCTAACTC           | 0QT\VNO^IIRJKXTIHIKTY\STZZ[XOJKPWLHJQQQ^XPQYIIRQII           | PG:Z:bfast   | AS:i:-                     | NH:i:0   | IH:i:1   | HI:i:1                                     | CS:Z:T10213222020221330022203222011113021130222212230122             | CQ:Z:019=2+><%@.'>52%)'15?90<7;:5+(-49('-7-<>5-@5%<.2%=            | CM:i:0   | XA:i:4
    
    Example mapped read:
    2216_1487_198     |  16 | chr11:19915291-19925392      |  5531 |   12 | 24M2D26M           | *  |  0 |  0 | TCTGTCTGGGTGCAGTAGCTATACGTAATCCCAGCACTTTGGGAGGCCAA           | 1C8FZ`M""""""WZ""%"#\I;"`R@D]``ZZ``\QKX\Y]````IK^`           | PG:Z:bfast   | AS:i:1150   | NM:i:2   | NH:i:10   | IH:i:1   | HI:i:1   | MD:Z:24^CT26                   | CS:Z:T00103022001002113210023031213331122121223321221122             | CQ:Z:A8=%;AB<;97<3/9:>>3>@?5&18@-%;866:94=14:=?,%?:7&)1            | CM:i:9   | XA:i:4   | XE:Z:-------23322---2-11----2----------------------------               
    
    '''
    # TODO: rewrite read properties to unmapped state
    #       if colorspace: convert CS to NT directly
    #       remove excess tags
    
    read.is_unmapped = True
    read.rname = None
    read.pos = 0
    read.mapq = 0
    return read

if __name__ == '__main__':
    infile = None
    outfile = None
    failed = None
    criteria = []

    crit_args=[]
    last = None
    verbose = False
    for arg in sys.argv[1:]:
            
        if last == '-failed':
            failed = arg
            last = None
        elif arg == '-h':
            usage()
        elif arg == '-failed':
            last = arg
        elif arg == '-v':
            verbose = True
        elif not infile:
            infile = arg
        elif not outfile:
            outfile = arg
        elif arg[0] == '-':
            if not arg[1:] in _criteria:
                print "Unknown criterion: %s" % arg
                usage()
            if crit_args:
                criteria.append(_criteria[crit_args[0][1:]](*crit_args[1:]))
            crit_args = [arg,]
        elif crit_args:
            crit_args.append(arg)

    if crit_args:
        criteria.append(_criteria[crit_args[0][1:]](*crit_args[1:]))
    
    if not infile or not outfile or not criteria:
        usage()
    else:
        bam_filter(infile,outfile,criteria,failed,verbose)
        
