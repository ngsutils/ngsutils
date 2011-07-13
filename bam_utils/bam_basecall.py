#!/usr/bin/env python
"""
Given a BAM file and a genomic reference, for each position covered in the 
BAM file, show the reference base, and the number of A/T/C/G's and InDels. 

You can also optionally filter out all bases whose quality score is below a 
threshold, bases that aren't covered by enough reads, bases that have no 
variation compared to reference, or bases whose variation is too low.

The output is a tab-delimited file that contains the following for each base:

chromosome
position (1-based)
reference base
# reads that contain this base
Entropy
# A calls
# C calls
# G calls
# T calls
# deletions
# gaps
# inserts

Entropy is sum(a..t) { p log2 p } where p = freq(+pseudocount) / genomic freq.
pseudo count = genomic freq * sqrt(N)

We use the following genomic frequencies: A 0.3, C 0.2, G 0.2, T 0.3
"""

import os,sys,math,collections
from support.eta import ETA
import pysam

def usage():
    base = os.path.basename(sys.argv[0])
    print __doc__
    print """
Usage: %s {opts} in.bam {chrom:start-end}

Options:
-ref   val    Include reference basecalls from this file
-qual  val    Minimum quality level to use in calculations
              (numeric, Sanger scale) (default 0)
            
-count val    Report only bases with this minimum number of reads covering it
              (matches, inserts, deletions counted) (default 0)

-mask  val    The bitmask to use for filtering reads (default 1540)

-showgaps     Report gaps/splice-junctions in RNA-seq data

""" % (base)
    sys.exit(1)

__genomic_freq = {'A': 0.3, 'C': 0.2, 'G': 0.2,'T': 0.3}

def calc_entropy(a,c,t,g):
    counts = {'A':a,'C':c,'G':g,'T':t,}

    N = counts['A'] + counts['C'] + counts['G'] + counts['T']
    if N == 0:
        return 0
        
    N_sqrt = math.sqrt(N)
    
    count_pseudo = {}
    N_pseudo = 0
    
    for base in 'ATCG':
        count_pseudo[base] = counts[base] + (__genomic_freq[base] * N_sqrt)
        N_pseudo += count_pseudo[base]

    acc = 0
    for base in 'ATCG':
        p = float(count_pseudo[base]) / N_pseudo / __genomic_freq[base]
        acc += (p * math.log(p,2))
    
    return acc

MappingRecord = collections.namedtuple('MappingRecord','qpos cigar_op base qual read')
MappingPos = collections.namedtuple('MappingPos','tid pos records')
BasePosition = collections.namedtuple('BasePosition','tid pos total a c g t n deletions gaps insertions reads')

class BamBaseCaller(object):
    def __init__(self, bam_fname, min_qual=0, min_count=0, chrom=None, start=0, end=0, mask=1540,quiet=False):
        self.bam = pysam.Samfile(bam_fname,'rb')
        self.min_qual = min_qual
        self.min_count = 0

        self.chrom = chrom
        self.start = start
        self.end = end

        self.mask = mask
        self.quiet = quiet
        

        if chrom and start and end:
            def _gen():
                for p in self.bam.fetch(chrom,start,end):
                    yield p
        else:
            def _gen():
                for p in self.bam:
                    yield p
        self._gen = _gen

        self.buffer = None
        self.current_tid = None
    
    def close(self):
        self.bam.close()
        
    def _calc_pos(self,tid,pos,records):
        counts = {'A':0,'C':0,'G':0,'T':0,'N':0}
        insertions = {}
        deletions = 0
        gaps = 0
        total = 0
        reads = []

        for record in records:
            qpos,cigar_op,base,qual,read = record
            if cigar_op == 0: # M
                if qual >= self.min_qual and (read.flag & self.mask) == 0 :
                    total += 1
                    reads.append(record)

                    counts[base] += 1
            elif cigar_op == 1: # I
                if qual >= self.min_qual and (read.flag & self.mask) == 0 :
                    reads.append(record)

                    if not base in insertions:
                        insertions[base] = 1
                    else:
                        insertions[base] += 1
            elif cigar_op == 2: # D
                deletions += 1
                #total += 1 # not sure these should be included, 
                           # samtools mpileup includes them
                           # IGV doesn't
                           
                reads.append(record)
            elif cigar_op == 3: # N
                gaps += 1
                reads.append(record)
        
        if total >= self.min_count:
            return BasePosition(tid,pos,total,counts['A'],counts['C'],counts['G'],counts['T'],counts['N'],deletions,gaps,insertions,reads)
            

    def fetch(self):
        self.current_tid = None
        self.buffer = collections.deque()
        if not self.quiet:
            eta = ETA(0,bamfile=self.bam)
        else:
            eta = None

        for read in self._gen():
            if eta:
                eta.print_status(extra='%s:%s (%s)' % (self.bam.references[read.tid],read.pos,len(self.buffer)),bam_pos=(read.tid,read.pos))
            
            if self.current_tid != read.tid: # new chromosome
                while self.buffer:
                    tid,pos,records = self.buffer.popleft()
                    yield self._calc_pos(tid,pos,records)
                
                self.current_tid = read.tid
            
            # handle all positions that are 5' of the current one
            while self.buffer and read.pos > self.buffer[0].pos:
                tid,pos,records = self.buffer.popleft()
                yield self._calc_pos(tid,pos,records)
                
            self._push_read(read)
        
        # flush buffer for the end
        while self.buffer:
            tid,pos,records = self.buffer.popleft()
            yield self._calc_pos(tid,pos,records)
        
        if eta:
            eta.done()
        
    def _push_read(self,read):
        
        if not self.buffer:
            self.buffer.append(MappingPos(read.tid, read.pos, []))
        
        while self.buffer[-1].pos < read.aend:
            self.buffer.append(MappingPos(read.tid, self.buffer[-1].pos+1, []))
        
        buf_idx = 0
        while self.buffer[0].pos < read.pos:
            buf_idx += 1
            
        read_idx = 0
        for op,length in read.cigar:
            if op == 0:
                for i in xrange(length):
                    self.buffer[buf_idx].records.append(MappingRecord(read_idx,op,read.seq[read_idx],read.qual[read_idx],read))
                    buf_idx += 1
                    read_idx += 1
            elif op == 1:
                inseq = ''
                inqual = 0
                for i in xrange(length):
                    inseq += read.seq[read_idx]
                    inqual += ord(read.qual[read_idx])-33
                    read_idx += 1
                
                inqual = inqual / len(inseq)
                
                self.buffer[buf_idx].records.append(MappingRecord(read_idx,op,inseq,inqual,read))
                
            elif op == 2:
                for i in xrange(length):
                    self.buffer[buf_idx].records.append(MappingRecord(read_idx,op,None,None,read))
                    buf_idx += 1
            
            elif op == 3:
                for i in xrange(length):
                    self.buffer[buf_idx].records.append(MappingRecord(read_idx,op,None,None,read))
                    buf_idx += 1
                
        
def bam_basecall(bam_fname,ref_fname,min_qual=0, min_count=0, chrom=None,start=None,end=None,mask=1540,quiet = False, showgaps=False):
    if ref_fname:
        ref = pysam.Fastafile(ref_fname)
    else:
        ref = None

    sys.stdout.write('chrom\tpos\tref\tcount\tave mappings\tentropy\tA\tC\tG\tT\tN\tDeletions\tGaps\tInsertions\tInserts\n')
    bbc = BamBaseCaller(bam_fname,min_qual,min_count,chrom,start,end,mask,quiet)
    for basepos in bbc.fetch():
        if start and end:
            if basepos.pos < start or basepos.pos >= end:
                continue
        
        big_total = basepos.total + basepos.deletions + len(basepos.insertions)
        
        if big_total < min_count:
            continue
        
        if big_total == 0 and not (showgaps and basepos.gaps > 0):
            continue
            
        if ref:
            refbase = ref.fetch(bbc.bam.references[basepos.tid],basepos.pos,basepos.pos+1).upper()
        else:
            refbase = 'N'

        entropy = calc_entropy(basepos.a,basepos.c,basepos.g,basepos.t)
    
        read_ih_acc = 0
        for qpos,cigar_op,base,qual,read in basepos.reads:
            if cigar_op in [0,1,2]:
                read_ih_acc += int(read.opt('IH'))
        
        inserts = []
        for insert in basepos.insertions:
            inserts.append((basepos.insertions[insert],insert))
        inserts.sort()
        inserts.reverse()
        
        insert_str_ar = []
        incount = 0
        for count,insert in inserts:
            insert_str_ar.append('%s:%s' % (insert,count))
            incount += count
        
        if big_total > 0:
            ave_mapping = (float(read_ih_acc) / big_total)
        else:
            ave_mapping = 0
        
        sys.stdout.write('%s\n' % '\t'.join([str(x) for x in (bbc.bam.references[basepos.tid],basepos.pos+1,refbase,basepos.total,ave_mapping,entropy,basepos.a,basepos.c,basepos.g,basepos.t,basepos.n,basepos.deletions,basepos.gaps,incount,','.join(insert_str_ar))]))

    bbc.close()
    if ref:
        ref.close()

if __name__ == '__main__':
    bam = None
    ref = None
    
    min_qual = 0
    min_count = 0
    mask = 1540
    chrom = None
    start = None
    end = None
    quiet = False
    showgaps = False
    
    last = None
    for arg in sys.argv[1:]:
        if last == '-qual':
            min_qual = int(arg)
            last = None
        elif last == '-ref':
            if os.path.exists(arg) and os.path.exists('%s.fai' % arg):
                ref = arg
            else:
                print "Missing FASTA file or index: %s" % arg
                usage()
            last = None
        elif last == '-count':
            min_count = int(arg)
            last = None
        elif last == '-mask':
            mask = int(arg)
            last = None
        elif arg == '-h':
            usage()
        elif arg == '-showgaps':
            showgaps = True
        elif arg == '-q':
            quiet = True
        elif arg in ['-qual','-count','-mask','-ref']:
            last = arg
        elif not bam and os.path.exists(arg) and os.path.exists('%s.bai' % arg):
            bam = arg
        elif not ref and os.path.exists(arg) and os.path.exists('%s.fai' % arg):
            ref = arg
        elif not chrom:
            chrom,startend = arg.split(':')
            if '-' in startend:
                start,end = [int(x) for x in startend.split('-')]
            else:
                start = int(startend)
                end = start
            start = start - 1
        else:
            print "Unknown option or missing index: %s" % arg
            usage()

    if not bam:
        usage()
    else:
        bam_basecall(bam,ref,min_qual,min_count,chrom,start,end,mask,quiet,showgaps)
        
