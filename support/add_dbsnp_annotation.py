#!/usr/bin/env python
'''
Given a tab-delimited file with base positions, this will add data 
from dbSNP, if there is a known SNP at a given position.

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
from ngs_utils import gzip_opener

def read_dbsnp(fname):
    '''
    reads a dbSNP flat-file (from UCSC) completely into memory
    returns a dict keyed on (chrom,pos) with values [(name,strand,observed,freq,function),] (list of tuples)
    '''
    dbsnp = {}
    
    with gzip_opener(fname) as f:
        for line in f:
            cols = line.strip().split('\t')
            k = (cols[1],int(cols[3]))
            if not k in dbsnp:
                dbsnp[k] = []
            dbsnp[k].append((cols[4],cols[6],cols[9],cols[13],cols[15]))
            
    return dbsnp
    
def annotate_tab(fname,dbsnp,header=True,insert_at_col=3):
    with gzip_opener(fname) as f:
        for line in f:
            cols = line.strip().split('\t')
            if header:
                header = False
                sys.stdout.write('\t'.join(cols[:insert_at_col]))
                sys.stdout.write('\trs\tstrand\tobserved\tfreq\tfunc\t')
                sys.stdout.write('\t'.join(cols[insert_at_col:]))
                sys.stdout.write('\n')
            else:
                rs = ""
                strand = ""
                observed = ""
                freq = ""
                func = ""
                
                k = (cols[0],int(cols[1]))
                
                if k in dbsnp:
                    rs,strand,observed,freq,func = dbsnp[k]
                
                sys.stdout.write('\t'.join(cols[:insert_at_col]))
                sys.stdout.write('\t%s\t%s\t%s\t%s\t%s\t' % (rs,strand,observed,freq,func))
                sys.stdout.write('\t'.join(cols[insert_at_col:]))
                sys.stdout.write('\n')

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
    for arg in sys.argv[1]:
        if last == '-insertcol':
            insert_col = int(arg)
            last = None
        elif arg in ['-insertcol']:
            last = arg
        elif arg == '-noheader':
            header = False
        elif not dbsnp_fname and os.path.exists(arg):
            dbsnp_fname = arg
        elif not tab_fname and os.path.exists(arg):
            tab_fname = arg
        else:
            sys.stderr.write("Unknown option or missing file: %s\n" % arg)
            usage()
            
    if not tab_fname and not dbsnp_fname:
        usage()
        
    annotate_tab(tab_fname,read_dbsnp(dbsnp_fname),header,insert_col)
    