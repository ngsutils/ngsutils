#!/usr/bin/env python
"""
Merges a (cs)fasta file and a qual file into a fastq file

This assumes that the fasta and qual files have the same reads in the same 
order (so we don't have to load all the data and merge it).  Unless the 
results are being written to stdout, it creates a new fastq file with the 
same basename as the fasta file.
"""

import sys,os,gzip

try:
    sys.path.append(os.path.normpath(os.path.join(os.path.dirname(__file__),'..')))
    from support.eta import ETA
except:
    ETA = None

def usage():
    print __doc__
    print """
Usage: %s { -tag tag } {-c} {-z} {-q val} {-maxlen val} filename.[cs]fasta [filename.qual]

-q val       Use a constant value for quality (should be char, such as '*'=10)
-maxlen VAL  Trim sequences and quality scores to be a maximum of VAL bases
-c           output to stdout (can only have one fasta file)
-z           output gzip compressed files
-tag tag     add a tag as a suffix to all read names like: readname:suffix

""" % os.path.basename(sys.argv[0])


def qual_to_phred(qual):
    if type(qual) is str:
        qual = qual.strip().split()
    elif type(qual) is not type([]):
        qual = [qual,]
    
    quals = []
    for q in qual:
      if not q:
        q = 0
      else:
        try:
          q = int(q)
        except ValueError, e:
          print e
          print "Line: %s" % qual
          sys.exit(2) 
        if q > 93:
            q = 93
        elif q < 0:
            q = 0
      quals.append(q)
    
    return ''.join(['%c' % (q +33) for q in quals])

def getline(fs):
    line = fs.readline()
    while line and line.startswith('#'):
        line = fs.readline()
    if line:
        line = line.strip()
    return line

def merge_files(fasta,qual,suffix=None,stdout=False, gz=False, common_qual=None,max_len=None):
    sys.stderr.write('Merging %s and %s -> ' % (os.path.basename(fasta),os.path.basename(qual) if qual else common_qual))
    if stdout:
        sys.stderr.write('stdout')
        if gz:
            out = gzip.GzipFile(fileobj=sys.stdout)
        else:
            out = sys.stdout
    else:
      if fasta.lower().endswith('fasta'):
        outfile = fasta[:-1]+'q'
      else:
        outfile = fasta+'.fastq'
      if gz:
        outfile = '%s.gz' % outfile
        out = gzip.open(outfile,'w')
      else:
        out = open(outfile,'w')
      sys.stderr.write(os.path.basename(outfile))

    if gz:
        sys.stderr.write(' [gzip]')
    sys.stderr.write('\n')
    
    f = open(fasta)

    if qual:
        q = open(qual)
    elif common_qual:
        q = None
    else:
        sys.stderr.write('Must specify a common qual value or a qual file!')
        sys.exit(-1)
    
    f_line = getline(f)
    if q:
        q_line = getline(q)
    else:
        q_line = None
    
    if ETA:
        eta = ETA(os.stat(fasta).st_size,fileobj=f)
    else:
        eta = None
    
    colorspace = None
    
    while f_line and (not q or q_line):
        if q:
            assert f_line == q_line
        
        name = trim_name(f_line[1:])

        if eta:
            eta.print_status(extra=name)
        

        if suffix:
            name = "%s:%s" % (name,suffix)
        seq = getline(f)

        if colorspace is None:
            colorspace = False
            for base in seq[1:]:
                if base in "01234":
                    colorspace = True
                    break

        if q:
            qual = qual_to_phred(getline(q))
        elif colorspace:
            qual = common_qual * (len(seq)-1)
        else:
            qual = common_qual * len(seq)
            
        if not maxlen:
            out.write('@%s\n%s\n+\n%s\n' % (name,seq,qual))
        elif colorspace and seq[0].upper() in 'ATGC':
            out.write('@%s\n%s\n+\n%s\n' % (name,seq[:maxlen+1],qual[:maxlen]))
        else:
            out.write('@%s\n%s\n+\n%s\n' % (name,seq[:maxlen],qual[:maxlen]))
        
        f_line = getline(f)
        if q:
            q_line = getline(q)

    if eta:
        eta.done()
    
    f.close()
    if q:
        q.close()
    if out != sys.stdout:
      out.close()

def trim_name(name):
    ''' remove trailing _F3 / _R3 (to match bfast solid2fastq script) '''
    if '_F3' in name:
        return name[:name.rindex('_F3')]
    elif '_R3' in name:
            return name[:name.rindex('_R3')]
    return name

if __name__ == '__main__':
    last = None
    tag = None
    gz = False
    stdout = False
    fasta = None
    qual = None
    common_qual = None
    maxlen = None
    
    for arg in sys.argv[1:]:
        if arg in ['-tag','-q','-maxlen']:
            last = arg
        elif last == '-tag':
            tag = arg
            last = None
        elif last == '-maxlen':
            maxlen = int(arg)
            last = None
        elif last == '-q':
            if len(arg) == 1:
                common_qual = arg
            else:
                sys.stderr.write("ERROR: A common qual value must be only one character\n")
                sys.exit(-1)
            last = None
        elif arg == '-z':
            gz = True
        elif arg == '-c':
            stdout = True
        elif not fasta and os.path.exists(arg):
            fasta = arg
        elif not qual and os.path.exists(arg):
            qual = arg
        else:
            sys.stderr.write("ERROR: Unknown option: %s\n" % arg)
            usage()
            sys.exit(-1)

    if not fasta:
        usage()
        sys.exit(-1)
        
    merge_files(fasta,qual, tag,stdout,gz,common_qual,maxlen)
        
