#!/usr/bin/env python
"""
Calculates a minor allele frequency in pooled genomic sequencing.

Given a BAM file and a genomic reference, for each position covered in the 
BAM file, show the reference base, the potential minor allele, and probable
background.

This assumes that if a SNP exists, there is likely only one possible variation.
So, this calculation will fail if more than one minor allele is present.  This
also ignores indels.

Outputs a tab-delimited file with the following columns:
chrom
position (1-based)
reference base
alternative base
total counts
ref counts (1st or second highest)
alternative counts (1st or second highest)
background (3rd highest count)
ref-background
alt-background
"""

import os,sys,math
from support.eta import ETA
import pysam
try:
    __rsrc = os.path.join(os.path.dirname(__file__),'minorallele_cpci.R')
    import rpy2.robjects as robjects
    with open(__rsrc) as f:
        robjects.r(f.read())
except Exception:
    robjects = None

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

-alleles val  The number of alleles included in this sample
              If given, a Clopper-Pearson style confidence interval will be 
              calculated. (requires rpy2)
""" % (base)
    if robjects:
        print "rpy2 detected!"
    else:
        print "rpy2 not detected!"
    
    sys.exit(1)


def bam_minorallele(bam_fname,ref_fname,min_qual=0, min_count=0, num_alleles = 0):
    bam = pysam.Samfile(bam_fname,"rb")
    ref = pysam.Fastafile(ref_fname)
    eta = ETA(0,bamfile=bam)
    sys.stdout.write('chrom\tpos\tref\talt\ttotal\tref count\talt count\tbackground count\tref-background\talt-background')
    if robjects:
        sys.stdout.write('\tMean level\t95% CI low\t95% CI high\tCI Range\tlow count\thigh count')
    sys.stdout.write('\n')

    printed = False
    for pileup in bam.pileup():
        chrom = bam.getrname(pileup.tid)
        eta.print_status(extra='%s:%s' % (chrom,pileup.pos),bam_pos=(pileup.tid,pileup.pos))
        
        counts = {'A':0,'C':0,'G':0,'T':0}
        inserts = 0
        deletions = 0
        total = 0
        
        for pileupread in pileup.pileups:
            if not pileupread.is_del:
                if min_qual:
                    if pileupread.alignment.qual[pileupread.qpos] < min_qual:
                        continue
                if pileupread.indel == 0:
                    base = pileupread.alignment.seq[pileupread.qpos].upper()
                    if base != 'N':
                        counts[base]+=1
                        total += 1
        
        if total > min_count:
            refbase = ref.fetch(chrom,pileup.pos,pileup.pos+1).upper()
            if not refbase in counts:
                continue

            scounts = []
            for c in counts:
                scounts.append((counts[c],c))
            scounts.sort()
            scounts.reverse()
            
            if scounts[0][1] == refbase:
                refcount = scounts[0][0]
                altbase = scounts[1][1]
                altcount = scounts[1][0]
            else:
                refcount = scounts[1][0]
                altbase = scounts[0][1]
                altcount = scounts[0][0]
            
            background = scounts[2][0]
            refback = refcount-background
            altback = altcount-background

            cols = [chrom,(pileup.pos+1),refbase,altbase,total,refcount,altcount,background,refback,altback]
            if robjects and num_alleles:
                ci_low,ci_high = calc_cp_ci(refback+altback,altback,num_alleles)
                cols.append(float(altback) / (refback+altback))
                cols.append(ci_low)
                cols.append(ci_high)
                cols.append(ci_high-ci_low)
                cols.append(ci_low * (num_alleles))
                cols.append(ci_high * (num_alleles))
            
            sys.stdout.write('%s\n' % '\t'.join([str(x) for x in cols]))
    
    eta.done()
    bam.close()
    ref.close()

def calc_cp_ci(N,count,num_alleles):
    if robjects:
        return robjects.r['CP.CI'](N,count,num_alleles)
    return None,None

if __name__ == '__main__':
    bam = None
    ref = None
    
    min_qual = 0
    min_count = 0
    num_alleles = 0
    
    last = None
    for arg in sys.argv[1:]:
        if last == '-qual':
            min_qual = int(arg)
            last = None
        elif last == '-count':
            min_count = int(arg)
            last = None
        elif last == '-alleles':
            num_alleles = int(arg)
            last = None
        elif arg == '-h':
            usage()
        elif arg in ['-qual','-count','-alleles']:
            last = arg
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
        bam_minorallele(bam,ref,min_qual,min_count,num_alleles)
        
