#!/usr/bin/env python
'''
Takes a RefIso gene model and convert each exon to a BED region
'''

import os,sys
import support.refiso

def seq_tobed_gene(fname):
    refiso = support.refiso.RefIso(fname)
    for gene in refiso.genes:
        sys.stdout.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, gene.tx_start, gene.tx_end, gene.name, 0, gene.strand]]))

def seq_tobed_exon(fname):
    refiso = support.refiso.RefIso(fname)
    for gene in refiso.genes:
        for i,start,end,const,names in gene.regions:
            sys.stdout.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, start, end, '%s.%s' % (gene.name, i), 0, gene.strand]]))

def usage():
    print __doc__
    print """Usage: sequtils tobed [-gene|-exon] refiso.txt{.gz}
Options:
  -gene    Convert the whole gene to a BED region (txStart->txEnd)
  -exon    Convert each exon to a BED region
"""
    sys.exit(1)

if __name__ == '__main__':
    model = None
    gene = False
    exon = False
    
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif arg == '-gene':
            gene = True
        elif arg == '-exon':
            exon = True
        elif not model and os.path.exists(arg):
            model = arg
        else:
            sys.stderr.write('Unknown option: %s\n' % arg)
            usage()
    
    if not model:
        usage()

    if gene:
        seq_tobed_gene(model)
    elif exon:
        seq_tobed_exon(model)
    else:
        usage()