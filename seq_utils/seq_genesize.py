#!/usr/bin/env python
'''
Extracts the genomic and transcript sizes for genes in RefFlat format
'''

import sys,os

def usage():
    print __doc__
    print 'Usage: sequtils genesize refiso.txt' % os.path.basename(sys.argv[0])
    sys.exit(1)

def seq_genesize(fname):
    genes = {}
    with open(fname) as f:
        for line in f:
            if line[0] == '#':
                continue
            
            cols = line.strip().split('\t')
            
            dna_size = int(cols[5]) - int(cols[4])
            coding_len = 0
            for start,end in zip([int(x) for x in cols[9].split(',') if x],[int(x) for x in cols[10].split(',') if x]):
                coding_len += end-start
            
            if not cols[0] in genes:
                genes[cols[0]] = (dna_size,coding_len)
            else:
                max_dna = max(dna_size,genes[cols[0]][0])
                max_coding = max(coding_len,genes[cols[0]][1])
                genes[cols[0]] = (max_dna,max_coding)
    
    sys.stdout.write('gene\tgenomic size\ttranscript size\n')
    for gene in genes:
        sys.stdout.write('%s\t%s\t%s\n' % (gene,genes[gene][0],genes[gene][1]))

if __name__ == '__main__':
    fname = None
    for arg in sys.argv[1:]:
        if not fname and os.path.exists(arg):
            fname = arg
        else:
            usage()
            
    if not fname:
        usage()
        
    seq_genesize(fname)
