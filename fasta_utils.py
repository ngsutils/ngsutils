#!/usr/bin/env python

import sys,gzip,os,hashlib
sys.path.append(os.path.join(os.path.dirname(__file__),"..","ext"))
from eta import ETA
import ngs_utils

import pysam

def strip_random(fname):
    if fname == '-':
        f = sys.stdin
    else:
        f = open(fname)
    
    good = False
    for line in f:
        if line[0] == '>':
            sys.stderr.write(line[:-1])
            
            if '_' in line:
                good = False
                sys.stderr.write(" *filtered*")
            else:
                good = True

            sys.stderr.write("\n")
        
        if good:
            sys.stdout.write(line)    
    
    if fname != '-':
        f.close()



def seq_length(fname):
    if fname == '-':
        f = sys.stdin
    elif fname[-3:] == '.gz':
        f = gzip.open(fname)
    else:
        f = open(fname)

    acc = 0
    vals = {}
    for line in f:
        if line[0] == '>':
            if acc:
                if not acc in vals:
                    vals[acc] = 0
                vals[acc] += 1
            acc = 0
        else:
            acc += len(line.strip())
    
    if acc:
        if not acc in vals:
            vals[acc] = 0
        vals[acc] += 1
    
    for k in vals:
        print "%s: %s" % (k,vals[k])

    if f!=sys.stdin:
        f.close()
    
def fasta_tag_parse(*args):
    lastarg = None
    
    prefix = None
    suffix = None
    fname = None
    
    for arg in args:
        if arg == '-prefix':
            lastarg = arg
        elif arg == '-suffix':
            lastarg = arg
        elif lastarg == '-prefix':
            prefix = arg
            lastarg = None
        elif lastarg == '-suffix':
            suffix = arg
            lastarg = None
        elif fname is None:
            fname = arg

    if ( not prefix and not suffix ) or not fname:
        sys.stderr.write("FASTA filename and a prefix or suffix is required (or both)!\n")
    else:
        fasta_tag(fname,prefix,suffix)

def fasta_tag(fname,prefix = None, suffix = None):
    f = open(fname)

    if not prefix:
        prefix = ''

    if not suffix:
        suffix = ''

    for line in f:
        if line[0] == '>':
            sys.stdout.write('>%s%s%s\n' % (prefix,line.strip()[1:],suffix))
        else:
            sys.stdout.write(line)
    f.close()

def filter_bases(fname):
    if fname != '-':
        f = open(fname)
    else:
        f = sys.stdin
        
    fout = sys.stdout
    for line in f:
        if line[0] == '>':
            fout.write(line)
        else:
            fout.write(re.sub('[^atgcnryswkmbdhvATGCNRYSWKMBDHV\n\r \.\-]','N',line))

    if f != sys.stdin:
        f.close()

def single_line(fname):
    if fname != '-':
        f = open(fname)
    else:
        f = sys.stdin
        
    fout = sys.stdout
    seq = ''
    for line in f:
        if line[0] == '>':
            if seq:
                fout.write('%s\n' % seq)
                seq = ''
            fout.write(line)
        else:
            seq += line.strip('\n').strip('\r')

    if seq:
        fout.write('%s\n' % seq)

    if f != sys.stdin:
        f.close()

def truncate(fname,max_len=50):
    if type(max_len) != int:
        max_len = int(max_len)
        
    if fname != '-':
        f = open(fname)
    else:
        f = sys.stdin

    fout = sys.stdout
    seq = ''
    for line in f:
        if line[0] == '>':
            if seq:
                fout.write('%s\n' % seq[:max_len])
                seq = ''
            fout.write(line)
        else:
            seq += line.strip()

    if seq:
        fout.write('%s\n' % seq[:max_len])

    if f != sys.stdin:
        f.close()

def head(fname,bases=50):
    if type(bases) != int:
        bases = int(bases)

    if fname != '-':
        f = open(fname)
    else:
        f = sys.stdin

    fout = sys.stdout
    seq = ''
    for line in f:
        if line[0] == '>':
            if seq:
                fout.write('%s\n' % seq[:bases])
                seq = ''
            fout.write('%s-head-%s\n' % (line.strip(),bases))
        else:
            if len(seq) < bases:
                seq += line.strip()

    if seq:
        fout.write('%s\n' % seq[:bases])

    if f != sys.stdin:
        f.close()

def tail(fname,bases=50):
    if type(bases) != int:
        bases = int(bases)

    if fname != '-':
        f = open(fname)
    else:
        f = sys.stdin

    fout = sys.stdout
    seq = ''
    for line in f:
        if line[0] == '>':
            if seq:
                fout.write('%s\n' % seq[-bases:])
                seq = ''
            fout.write('%s-tail-%s\n' % (line.strip(),bases))
        else:
            seq += line.strip()
            
            if len(seq) > bases:
                seq = seq[-bases:]

    if seq:
        fout.write('%s\n' % seq[-bases:])

    if f != sys.stdin:
        f.close()

def pretty(fname,length=50):
    if type(length) != int:
        length = int(length)

    if fname != '-':
        f = open(fname)
    else:
        f = sys.stdin

    fout = sys.stdout
    seq = ''
    for line in f:
        if line[0] == '>':
            if seq:
                while len(seq) > length:
                    fout.write('%s\n' % seq[:length])
                    seq = seq[length:]
                fout.write('%s\n' % seq)
                seq = ''
            fout.write('%s' % line)
        else:
            seq += line.strip()

            while len(seq) > length:
                fout.write('%s\n' % seq[:length])
                seq = seq[length:]

    if seq:
        while len(seq) > length:
            fout.write('%s\n' % seq[:length])
            seq = seq[length:]

    if f != sys.stdin:
        f.close()
def sort_fasta(fname):
    seqs = {}
    name = ''
    names = []

    with ngs_utils.gzip_opener(fname) as f:
        eta = ETA(os.stat(fname).st_size,fileobj=f)
        for line in f:
            eta.print_status(extra=name)
            if line[0] == '>':
                name = line[1:].split(' ')[0]
                if name in seqs:
                    sys.stderr.write('Error! Duplicate sequence name: %s\n' % name)
                    sys.exit(1)
                names.append(name)
                seqs[name] = line
            else:
                seqs[name]+=line

    eta.done()
    names.sort()
    for name in names:
        sys.stdout.write(seqs[name])

def seq_dict(fname):
    name = None
    size = None
    digester = None
    
    sys.stdout.write('@HD\tVN:1.0\tSO:unsorted\n')
    
    with ngs_utils.gzip_opener(fname) as f:
        eta = ETA(os.stat(fname).st_size,fileobj=f)
        for line in f:
            eta.print_status(extra=name)
            if line[0] == '>':
                if digester:
                    sys.stdout.write('@SQ\tSN:%s\tLN:%s\tUR:file:%s\tM5:%s\n' % (name,size,os.path.abspath(fname),digester.hexdigest()))
                name = line[1:].rstrip().split(' ')[0]
                digester = hashlib.md5()
                size = 0
            else:
                seq = line.strip().upper()
                digester.update(seq)
                size += len(seq)
                
        sys.stdout.write('@SQ\tSN:%s\tLN:%s\tUR:file:%s\tM5:%s\n' % (name,size,fname,digester.hexdigest()))
        eta.done()


if __name__ == "__main__":
    prgs = {
        'tag':(fasta_tag_parse,"Tag each sequence's ID with a prefix or suffix", [('fasta_file','Name of the FASTA file'),('{-prefix prefix}','Prefix to use'),('{-suffix suffix}','Suffix to use')]),
        'sort':(sort_fasta,"Sorts the sequences in a file by name (use only on small files)", [('fasta_file','Name of the FASTA file')]),
        'filter_bases':(filter_bases,"Convert any non-[ATGC] characters to N's", [('fasta_file','Name of the FASTA file')]),
        'strip_random':(strip_random,"Remove random contigs from a FASTA file", [('fasta_file','Name of the FASTA file')]),
        'single_line':(single_line,"Reformat the FASTA file to have the sequences on a single line", [('fasta_file','Name of the FASTA file')]),
        'truncate':(truncate,"Truncate sequences to be at most N long", [('fasta_file','Name of the FASTA file'),('len','Max length')]),
        'head':(head,"First N bases of each sequence", [('fasta_file','Name of the FASTA file'),('bases','Number of bases to return from the start of each sequence')]),
        'tail':(tail,"Last N bases of each sequence", [('fasta_file','Name of the FASTA file'),('bases','Number of bases to return from the end of each sequence')]),
        'pretty':(pretty,"Reformat the sequences to be a common length", [('fasta_file','Name of the FASTA file'),('len','Length of each line (default 50)')]),
        'seq_length':(seq_length,"Tally the length of all of the sequences in a FASTA file", [('fasta_file','Name of the FASTA file')]),
        'dict':(seq_dict,"Produce a Picard-compatible seq dict in SAM format", [('fasta_file','Name of the FASTA file')]),
        }
        
    if len(sys.argv) > 1 and sys.argv[1] in prgs:
        if len(sys.argv) > 2:
            prgs[sys.argv[1]][0](*sys.argv[2:])
        else:
            print "%s - %s" % (sys.argv[1],prgs[sys.argv[1]][1])
            print ""
            print "Parameters:"
            for param, desc in prgs[sys.argv[1]][2]:
                print "%20s\t%s" % (param,desc)
        
    else:
        if len(sys.argv)>1:
            print "Unknown routine: %s" % sys.argv[1]
        print "Available routines:"
        for prg in prgs:
            print " %20s\t%s" % (prg,prgs[prg][1])
        
