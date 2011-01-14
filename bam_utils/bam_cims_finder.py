#!/usr/bin/env python
"""
Given a set of BAM files, we search for areas where there are an unusual 
amount of deletions.  For CLIP-Seq, this can be an indicator of the location 
of protein-RNA interaction.

Output is either a BED file or a FASTA format file containing these hotspots.
Only unique regions are returned across all files.
"""

import os,sys,math
sys.path.append(os.path.join(os.path.dirname(__file__),"..","utils")) #eta
sys.path.append(os.path.join(os.path.dirname(__file__),"..","ext")) #pysam

from eta import ETA
import pysam

def usage():
    base = os.path.basename(sys.argv[0])
    print __doc__
    print """
Usage: %s {opts} in.bam {in.bam...}

Options:
    -flanking N      The number of flanking bases on either side to report
                     [default: 12]
                     
    -fasta ref.fa    Ouput in FASTA format (requires reference genome.fa)
                     [default: BED output]
                     
    -cutoff N        Cut-off %% for deletions - if the %% of reads that 
                     include a deletion at a position is higher than this 
                     number, the fragment is reported (0->1.0)
                     [default: 0.1]
                     
    -ns              Don't take the strand of the read into account
    
    -window N        The maximum length of a deletion window
                     [default: 10]
""" % (base)
    sys.exit(1)

def _output_bed(chrom,start,end,num,strand,region_pcts):
    ave_acc = 0.0
    for pct in region_pcts:
        ave_acc += pct
    ave = int(100*ave_acc / len(region_pcts))
    
    
    sys.stdout.write('%s\t%s\t%s\tregion_%s\t%s\t%s\n' % (chrom,start,end,num,ave,strand))

def _output_fasta(chrom,start,end,strand,ref,flanking,region_pcts):
    ave_acc = 0.0
    for pct in region_pcts:
        ave_acc += pct
    ave = ave_acc / len(region_pcts)

    seq = ref.fetch(chrom,start-flanking,end+flanking)
    seq = '%s%s%s' % (seq[:flanking].upper(),seq[flanking:end-start+flanking].lower(),seq[-flanking:].upper())

    if strand == '-':
        rc = []
        for base in seq:
            if base == 'A':
                rc.append('T')
            elif base == 'T':
                rc.append('A')
            elif base == 'C':
                rc.append('G')
            elif base == 'G':
                rc.append('C')
            elif base == 'a':
                rc.append('t')
            elif base == 't':
                rc.append('a')
            elif base == 'c':
                rc.append('g')
            elif base == 'g':
                rc.append('c')
        rc.reverse()
        seq=''.join(rc)
    
    sys.stdout.write('>%s:%s%s%s %s\n%s\n' % (chrom,start-flanking,strand,end+flanking,ave,seq))


def bam_cims_finder(bam_fnames,output='bed',ref_fname=None,flanking=12,cutoff=0.1,stranded=True,window_size=10):
    regions = []
    for bam_fname in bam_fnames:
        sys.stderr.write('%s\n' % bam_fname)
        bam = pysam.Samfile(bam_fname,"rb")
        if output == 'fasta' and ref_fname:
            ref = pysam.Fastafile(ref_fname)
        else:
            ref = None

        printed = False
    
        region_start = None
        region_end = None
        last_chrom = None
        last_strand = None
        region_pcts = []
        region_num = 1
    
        if stranded:
            strands = ['+','-']
        else:
            strands = ['+']
    
        for strand in strands:
            eta = ETA(0,bamfile=bam)
            for pileup in bam.pileup():
                chrom = bam.getrname(pileup.tid)
                eta.print_status(extra='%s:%s' % (chrom,pileup.pos),bam_pos=(pileup.tid,pileup.pos))
        
                if region_start and (chrom != last_chrom or pileup.pos-region_start > window_size):
                    if not (last_chrom,region_start,region_end) in regions:
                        if output == 'bed':
                            _output_bed(last_chrom,region_start,region_end,region_num,strand,region_pcts)
                        else:
                            _output_fasta(last_chrom,region_start,region_end,strand,ref,flanking,region_pcts)

                    regions.append((last_chrom,region_start,region_end))

                    region_num += 1
                    region_start = None
                    region_end = None
                    region_pcts = []
                
                last_chrom = chrom
        
                deletions = 0
                total = 0
        
                for pileupread in pileup.pileups:
                    if not stranded or (strand == '+' and not pileupread.alignment.is_reverse) or (strand == '-' and pileupread.alignment.is_reverse):
                        total += 1
                        if pileupread.is_del:
                            deletions += 1
                    if total > 0 and float(deletions)/total > cutoff:
                        if not region_start:
                            region_start = pileup.pos
                        region_end = pileup.pos+1
                        region_pcts.append(float(deletions)/total)

            if region_start and not (last_chrom,region_start,region_end) in regions:
                if output == 'bed':
                    _output_bed(last_chrom,region_start,region_end,region_num,strand,region_pcts)
                else:
                    _output_fasta(last_chrom,region_start,region_end,strand,ref,flanking,region_pcts)
                regions.append((last_chrom,region_start,region_end))
            eta.done()
        bam.close()

    if ref:
        ref.close()

if __name__ == '__main__':
    bams = []
    ref = None
    output = 'bed'
    cutoff = 0.1
    flanking = 12
    stranded = True
    window = 10
    
    last = None
    for arg in sys.argv[1:]:
        if last == '-flanking':
            flanking = int(arg)
            last = None
        elif last == '-cutoff':
            cutoff = float(arg)
            last = None
        elif last == '-window':
            window = float(arg)
            last = None
        elif last == '-fasta' and not ref and os.path.exists(arg) and os.path.exists('%s.fai' % arg):
            output = 'fasta'
            ref = arg
            last = None
        elif arg in ['-flanking','-fasta','-cutoff','-window']:
            last = arg
        elif arg == '-ns':
            stranded = False
        elif os.path.exists(arg) and os.path.exists('%s.bai' % arg):
            bams.append(arg)
        else:
            print "Unknown option or missing index: %s" % arg
            usage()

    if not bams:
        usage()
    else:
        bam_cims_finder(bams,output,ref,flanking,cutoff,stranded,window)
        
