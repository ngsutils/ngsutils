#!/usr/bin/env python
'''
Given a tab-delimited file with base positions, this will add data 
from dbSNP, if there is a known SNP at a given position.

dbSNP file must be tab-delimited (from UCSC), and Tabix indexed

Requires a file with the following columns:
chrom
position (1-based)
ref base

We will then insert the following columns (in place):
rs
strand
observed values
frequency (if available)
function

If there is more than one SNP, these will be comma-delimited
'''

import sys,os
sys.path.append(os.path.realpath(os.path.join(os.path.dirname(__file__),'..')))
sys.path.append(os.path.realpath(os.path.join(os.path.dirname(__file__),'..','ext')))
from support.ngs_utils import gzip_opener
from support.eta import ETA
import pysam

def annotate_tab(fname,dbsnp_fname,header=True,insert_at_col=3):
    dbsnp = pysam.Tabixfile(dbsnp_fname)
    
    with gzip_opener(fname) as f:
        eta = ETA(os.stat(fname).st_size, fileobj=f)
        for line in f:
            cols = line.strip().split('\t')
            eta.print_status(extra='%s:%s' % (cols[0],cols[1]))

            if header:
                header = False
                sys.stdout.write('\t'.join(cols[:insert_at_col]))
                sys.stdout.write('\trs\tstrand\tobserved\tfreq\tfunc\t')
                sys.stdout.write('\t'.join(cols[insert_at_col:]))
                sys.stdout.write('\n')
            else:
                rs = []
                strand = []
                observed = []
                freq = []
                func = []
                
                for line in dbsnp.fetch(cols[0],int(cols[1])-1,int(cols[1])):
                    tups = line.split('\t')
                    if int(cols[1]) == int(tups[3]):
                        rs.append(tups[4])
                        strand.append(tups[6])
                        observed.append(tups[9])
                        freq.append(tups[13])
                        func.append(tups[15])
                
                sys.stdout.write('\t'.join(cols[:insert_at_col]))
                sys.stdout.write('\t%s\t%s\t%s\t%s\t%s\t' % (','.join(rs),','.join(strand),','.join(observed),','.join(freq),','.join(func)))
                sys.stdout.write('\t'.join(cols[insert_at_col:]))
                sys.stdout.write('\n')
        eta.done()

def usage():
    print __doc__
    print """\
Usage: %s {-noheader} {-insertcol col} dbsnp.txt filename.tab

Merges dbSNP data into the given tab-delimited file.

Options:

-noheader        -  input file doesn't have a header line.
-insertcol col   -  insert dbSNP information at this column (default 3)

""" % os.path.basename(sys.argv[0])
    sys.exit(1)
    
if __name__ == '__main__':
    insert_col = 3
    tab_fname = None
    dbsnp_fname = None
    header = True
    
    last = None
    for arg in sys.argv[1:]:
        if last == '-insertcol':
            insert_col = int(arg)
            last = None
        elif arg in ['-insertcol']:
            last = arg
        elif arg == '-noheader':
            header = False
        elif not dbsnp_fname and os.path.exists(arg):
            if not os.path.exists('%s.tbi' % arg):
                sys.stderr.write('%s is not Tabix indexed\n' % arg)
                usage()
                
            dbsnp_fname = arg
        elif not tab_fname and os.path.exists(arg):
            tab_fname = arg
        else:
            sys.stderr.write("Unknown option or missing file: %s\n" % arg)
            usage()
            
    if not tab_fname and not dbsnp_fname:
        usage()
        
    annotate_tab(tab_fname,dbsnp_fname,header,insert_col)
