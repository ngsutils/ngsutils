#!/usr/bin/env python
"""
Given a BAM file and a genomic reference, for each position covered in the 
BAM file, show the reference base, and the number of A/T/C/G's and InDels. 

You can also optionally filter out all bases whose quality score is below a 
threshold, bases that aren't covered by enough reads, bases that have no 
variation compared to reference, or bases whose variation is too low.

The output is a tab-delimited file that contains the following for each base:

chromosome
position (0-based)
reference base (lowercase = no variation measured, uppercase = some variation)
# reads that contain this base
Variation [0..1]
Entropy
# A calls
# C calls
# G calls
# T calls
# inserts
# deletions

Variation is calculated as: (# calls not matching ref + # inserts + # dels)
                            -----------------------------------------------
                                 total calls + # inserts + # deletions
                                 
Entropy is sum(a..t) { p log2 p } where p = freq(+pseudocount) / genomic freq.
pseudo count = genomic freq * sqrt(N)

We use the following genomic frequencies: A 0.3, C 0.2, G 0.2, T 0.3

"""

import os,sys,math
from eta import ETA
import pysam

def usage():
    base = os.path.basename(sys.argv[0])
    print __doc__
    print """
Usage: %s {opts} in.bam ref.fa

Options:
-qual val     Minimum quality level to use in calculations
              (numeric, Sanger scale) (default 0)
            
-count val    Report only bases with this minimum number of reads covering it
              (default 0)

-noperfect    Don't output bases that perfectly match reference

-var val      Report only bases with this minimum variation
              Valid values: 0.0->1.0 (default 0)

""" % (base)
    sys.exit(1)

__genomic_freq = {'A': 0.3, 'C': 0.2, 'G': 0.2,'T': 0.3}

def calc_entropy(counts):
    '''
    counts = dict('A','T','C','G')
    '''
    for base in 'ATCG':
        if not base in counts:
            counts[base] = 0
    
    N = counts['A'] + counts['C'] + counts['G'] + counts['T']
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
    

def bam_basecall(bam_fname,ref_fname,min_qual=0, min_count=0, noperfect=False, min_var=0):
    bam = pysam.Samfile(bam_fname,"rb")
    ref = pysam.Fastafile(ref_fname)
    eta = ETA(0,bamfile=bam)
    sys.stdout.write('#chrom\tpos\tref\tcount\tave mappings\tvar\tentropy\tA\tC\tG\tT\tInserts\tDeletions\n')
    printed = False
    for pileup in bam.pileup():
        chrom = bam.getrname(pileup.tid)
        eta.print_status(extra='%s:%s' % (chrom,pileup.pos),bam_pos=(pileup.tid,pileup.pos))
        
        counts = {'A':0,'C':0,'G':0,'T':0}
        inserts = 0
        deletions = 0
        total = 0
        
        read_ih_acc = 0
        read_count = 0
        for pileupread in pileup.pileups:
            read_ih_acc += int(pileupread.alignment.opt('IH'))
            read_count += 1
            if pileupread.is_del:
                deletions += 1
            else:
                if min_qual:
                    if pileupread.alignment.qual[pileupread.qpos] < min_qual:
                        continue
                if pileupread.indel == 0:
                    base = pileupread.alignment.seq[pileupread.qpos].upper()
                    if base != 'N':
                        counts[base]+=1
                        total += 1
                elif pileupread.indel > 0:
                    inserts += 1
        
        if total > 0 or inserts > 0 or deletions > 0:
            refbase = ref.fetch(chrom,pileup.pos,pileup.pos+1).upper()
            if not refbase in counts:
                continue
                
            if noperfect and total == counts[refbase]:
                continue

            if total+inserts+deletions == 0:
                refbase = refbase.lower()
                var = 0
            elif total == counts[refbase]:
                refbase = refbase.lower()
                var = 0
            else:
                var = float(total - counts[refbase] + inserts + deletions)/(total + inserts + deletions)
            
            if total > 0:
                entropy = calc_entropy(counts)
            else:
                entropy = 0
            
            ave_mappings = float(read_ih_acc) / read_count
            
            if total >= min_count and var >= min_var:
                sys.stdout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom,pileup.pos,refbase,total,ave_mappings,var,entropy,counts['A'],counts['C'],counts['G'],counts['T'],inserts,deletions))
    
    eta.done()
    bam.close()
    ref.close()

if __name__ == '__main__':
    bam = None
    ref = None
    
    min_qual = 0
    min_count = 0
    min_var = 0
    noperfect = False
    
    last = None
    for arg in sys.argv[1:]:
        if last == '-qual':
            min_qual = int(arg)
            last = None
        elif last == '-count':
            min_count = int(arg)
            last = None
        elif last == '-var':
            min_var = float(arg)
            last = None
        elif arg == '-h':
            usage()
        elif arg in ['-qual','-count','-var']:
            last = arg
        elif arg == '-noperfect':
            noperfect = True
        elif not bam and os.path.exists(arg) and os.path.exists('%s.bai' % arg):
            bam = arg
        elif not ref and os.path.exists(arg) and os.path.exists('%s.fai' % arg):
            ref = arg
        else:
            print "Unknown option or missing index: %s" % arg
            usage()

    if not bam or not ref:
        usage()
    else:
        bam_basecall(bam,ref,min_qual,min_count,noperfect,min_var)
        
