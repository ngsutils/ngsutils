#!/usr/bin/env python
'''
A collection of utility functions for dealing with FASTQ files.
This is largely a combination of existing functions, brought together to be in
one place.

Marcus Breese <mbreese@iupui.edu>
Dec 16, 2010
'''

import sys,os,gzip
import ngs_utils,localalign
from eta import ETA
def read_fastq(fname,quiet=False,eta_callback=None):
    with ngs_utils.gzip_opener(fname) as f:
        if not quiet:
            eta = ETA(os.stat(fname).st_size,fileobj=f)
        while f:
            try:
                name = f.next().strip()
                seq = f.next().strip()
                f.next()
                qual = f.next().strip()
                
                if eta_callback:
                    extra = eta_callback()
                else:
                    extra = name
                if not quiet:
                    eta.print_status(extra=extra)
                yield (name,seq,qual)
            except Exception, e:
                break
    if not quiet:
        eta.done()

def fastq_truncate(fname,max_len=50):
    '''
    Takes each read in a FASTQ file and truncates it to be a max of {max_len}
    bases.
    '''
    if type(max_len) != int:
        max_len = int(max_len)

    cs = fastq_is_colorspace(fname,True)
    
    for name,seq,qual in read_fastq(fname):
        if cs and seq[0] in "atcgATGC":
            sys.stdout.write('%s\n%s\n+\n%s\n' % (name,seq[:max_len+1],qual[:max_len]))
        else:
            sys.stdout.write('%s\n%s\n+\n%s\n' % (name,seq[:max_len],qual[:max_len]))

def fastq_is_colorspace(fname,quiet=False):
    '''
    Scans reads in a FASTQ file to see if they are colorspace.  This works by
    scanning the first 10 reads that have sequences (aren't Ns or 4s). If there
    are any colorspace values, the entire file is called as colorspace.
    
    It's a bit overkill...
    
    '''
    seqs_checked=0
    is_colorspace = False

    valid_basespace="atcgATCG"
    valid_colorspace="0123456"

    for name,seq,qual in read_fastq(fname,quiet=True):
        if seqs_checked > 10:
            break
        checked = False
        for base in seq:
            if base in valid_colorspace:
                is_colorspace = True
                checked = True
            elif base in valid_basespace:
                checked = True
        if checked:
            seqs_checked += 1

    if not quiet:
        if is_colorspace:
            print "colorspace"
        else:
            print "basespace"
    
    return is_colorspace

def fastq_split(fname,outbase,chunksize=10000000):
    '''
    Splits a FASTQ file into multiple smaller files with a fixed number of
    reads per file.
    
    Output is a set of gzip compressed FASTQ files
    '''
    if type(chunksize) != int:
        chunksize = int(chunksize)

    i=0
    chunk=1

    out = gzip.open('%s.%s.fastq.gz' % (outbase,chunk),'w')

    def _get_chunk():
        return 'File #%s' % chunk

    for name,seq,qual in read_fastq(fname,eta_callback=_get_chunk):
        if i>=chunksize:
            out.close()
            chunk += 1
            out = gzip.open('%s.%s.fastq.gz' % (outbase,chunk),'w')
            i=0
        out.write('%s\n%s\n+\n%s\n' % (name,seq,qual))
        i+=1

    out.close()

def fastq_names(fname):
    '''
    Print the names of all of the reads in a FASTQ file
    '''
    for name,seq,qual in read_fastq(fname):
        sys.stdout.write('%s\n' % name[1:].split()[0]) # remove '@'

def fastq_check(fname):
    '''
    Validates a FASTQ file to confirm that it matches spec
    '''

    cs = fastq_is_colorspace(fname,True);

    line = 0
    for name,seq,qual in read_fastq(fname):
        if not name[0] == '@':
            print 'Invalid formatted record [line %s]' % line
            break
        
        if cs:
            if len(seq) != len(qual)+1:
                print 'Seq / qual mismatch [line %s]' % line
                break
        else:
            if len(seq) != len(qual):
                print 'Seq / qual mismatch [line %s]' % line
                break
        
        line += 4


def fastq_lengths(fname):
    '''
    Finds the length for each read in a FASTQ file and prints a distribution
    '''
    counts = {}
    for name,seq,qual in read_fastq(fname):
        lq=len(qual)
        if not lq in counts:
            counts[lq] = 0

        counts[lq]+=1

    sizes=[x for x in counts]
    sizes.sort()
    for s in sizes:
        sys.stdout.write('%s\t%s\n' % (s,counts[s]))

def fastq_trim_main(*args):
    '''
    Removes linkers from 5' and 3' ends by performing local 
    alignments with the linker sequences
    '''
    linker_5 = None
    linker_3 = None
    fastq = None
    last = None
    min_len = 25
    
    for arg in args:
        if last == '-5':
            linker_5 = arg
            last = None
        elif last == '-3':
            linker_3 = arg
            last = None
        elif last == '-min':
            min_len = int(arg)
            last = None
        elif arg in ['-3','-5','-min']:
            last = arg
        elif not fastq:
            fastq = arg
        else:
            usage('trim')
            
    if not fastq or not os.path.exists(fastq) or (not linker_5 and not linker_3):
        usage('trim')
    else:
        sys.stderr.write("Processing file: %s\n" % fastq)
        if linker_5:
            sys.stderr.write("Removing %s from 5' end\n" % linker_5)
        if linker_3:
            sys.stderr.write("Removing %s from 3' end\n" % linker_3)

        fastq_trim(fastq,linker_5,linker_3,min_len=min_len)

def fastq_trim(fname,linker_5=None,linker_3=None,out=sys.stdout,pct_identity=0.8,min_trim=4,min_len=25):
    '''
    fname - the fastq filename
    linker_5 - the 5' linker to remove
    linker_3 - the 3' linker to remove
    out - an output stream (eg: file, stdout)
    pct_identity - the percentage of matches that must be present in the alignment to strip away linkers
    min_trim - the distance away from the edges that the linkers much match w/in
    '''
    sw = localalign.LocalAlignment(localalign.NucleotideScoringMatrix(2,-1),-1)
    removed = 0
    trimmed = 0
    for name,seq,qual in read_fastq(fname):
        left = 0
        right = len(seq)

        if linker_5:
            aln = sw.align(seq,linker_5)
            if aln.r_pos < min_trim and aln.identity > pct_identity:
                left = aln.r_end

        if linker_3:
            aln = sw.align(seq,linker_3)
            if aln.r_end > len(seq)-min_trim and aln.identity > pct_identity:
                right = aln.r_pos

        s=seq[left:right]
        if len(s) >= min_len:
            if left > 0 or right < len(seq):
                trimmed += 1
            out.write('%s\n%s\n+\n%s\n' % (name,s,qual[left:right]))
        else:
            removed += 1

        # out.write('%s\n%s (%s-%s)\n' % (name,seq,left,right))
        # out.write('x'*left)
        # out.write(seq[left:right])
        # out.write('x' *(len(seq)-right))
        # out.write('\n')
    sys.stderr.write('Trimmed: %s\n' % trimmed)
    sys.stderr.write('Removed: %s (len)\n' % removed)

def fastq_convert_illumina_qual(fname):
    '''
    Converts the quality scaling from Illumina scale to Sanger scale
    (Doesn't adjust for Solexa vs. Phred quality calculations...)
    '''
    for name,seq,qual in read_fastq(fname):
        sys.stdout.write('@%s\n%s\n+\n%s\n' % (name,seq, ''.join([chr(ord(q)-31) for q in qual])))

def fastq_uniq(fname):
    '''
    Takes a FASTQ file and returns only unique sequences.
    The quality scores for each unique sequence are averaged across all
    reads with the same sequence.
    '''
    reads = {}

    for name,seq,qual in read_fastq(fname):
        if not seq in reads:
            reads[seq]={}
            reads[seq]['name'] = 'read_%s' % (len(reads))
            reads[seq]['names'] = []
            reads[seq]['quals'] = []

        reads[seq]['names'].append(name[1:].split()[0])
        reads[seq]['quals'].append(qual)

    sys.stderr.write('Outputting unique reads\n')

    f = open('%s.members.txt' % fname, 'w')
    for seq in reads:
        vals = reads[seq]
        sys.stdout.write("@%s\n%s\n+\n%s\n" % (vals['name'],seq,qual_mean(vals['quals'])))
        f.write('%s\t%s\t%s\n' % (vals['name'],len(vals['names']),','.join(vals['names'])))
    f.close()
    sys.stderr.write('%s uniq reads\n' % len(reads))

def fastq_qualtype(fname, num_to_check=10000):
    '''
    Checks a FASTQ file's quality score to see what encoding/scaling is used:
    Sanger, Solexa, or Illumina
    '''
    reads = {}
    
    # these are the differential values
    sanger = (33,73)
    solexa = (59,104)
    illumina = (64,104)

    sanger_count = 0
    solexa_count = 0
    illumina_count = 0
    
    checked = 0
    for name,seq,qual in read_fastq(fname,quiet=True):
        if checked > num_to_check:
            break
        qmax = None
        qmin = None
        for q in qual:
            if qmin is None or ord(q) < qmin:
                qmin = ord(q)
            if qmax is None or ord(q) > qmax:
                qmax = ord(q)
        if sanger[0] <= qmin <= qmax <= sanger[1]:
            sanger_count += 1
        elif solexa[0] <= qmin <= qmax <= solexa[1]:
            solexa_count += 1
        elif illumina[0] <= qmin <= qmax <= illumina[1]:
            illumina_count += 1
        checked += 1
    
    totals = [(sanger_count,'sanger'),(solexa_count,'solexa'),(illumina_count,'illumina')]
    totals.sort()
    print totals[-1][1]
    

def qual_mean(quals, offset=33):
    q_acc = [0,] * len(quals[0])
    for qual in quals:
        for i in xrange(len(qual)):
            q_acc[i] = ord(qual[i]) - offset
    return ''.join([ chr(int(x/len(quals))+offset) for x in q_acc])


# config: 
# { 'command': (func, description, [(arg_name, arg_desc)...])}

prgs = {
    'colorspace':(fastq_is_colorspace,"Determine if a FASTQ file is in color-space or base-space", [('fastq_file','Name of the FASTQ file'),]),
    'names':(fastq_names,"Print the names of all reads in a FASTQ file", [('fastq_file','Name of the FASTQ file'),]),
    'uniq':(fastq_uniq,"Reduce a FASTQ file to contain only unique sequences", [('fastq_file','Name of the FASTQ file'),]),
    'check':(fastq_check,"Verify a FASTQ formatted file", [('fastq_file','Name of the FASTQ file'),]),
    'qual_type':(fastq_qualtype,"Find out what encoding/scaling the quality scores use", [('fastq_file','Name of the FASTQ file'),]),
    'lengths':(fastq_lengths,"Find out how long reads in a file are", [('fastq_file','Name of the FASTQ file'),]),
    'convert_illumina':(fastq_convert_illumina_qual,"Converts Illumina-scaled quality values to Sanger scaling", [('fastq_file','Name of the FASTQ file'),]),
    'split':(fastq_split,"Split a FASTQ file into multiple smaller files", [
        ('fastq_file','Name of the FASTQ file'),
        ('out_template','The template for the name of new FASTQ files'),
        ('chunk_size','Number of sequences to include in a new file')
        ]),
    'trim':(fastq_trim_main,"Trims 5' and 3' linkers from sequences", [
        ('{-min len[25]}','Minimum length of a read (post-trim)'),
        ('{-5 seq}',"5' linker sequence to remove"),
        ('{-3 seq}',"3' linker sequence to remove"),
        ('fastq_file','Name of the FASTQ file'),
        ]),
    'truncate':(fastq_truncate,"Truncate sequences/quals to be at most N long", [('fastq_file','Name of the FASTQ file'),('len','Max length')]),
    }

def usage(cmd=None):
    if cmd and cmd in prgs:
        print "%s - %s" % (cmd,prgs[cmd][1])
        print prgs[cmd][0].__doc__
        print "Usage: %s %s %s" % (os.path.basename(sys.argv[0]),cmd,' '.join([x[0] for x in prgs[cmd][2]]))
        print ""
        print "Parameters:"
        for param, desc in prgs[cmd][2]:
            print "%20s\t%s" % (param,desc)
    else:
        if cmd:
            print "Unknown command: %s" % cmd
        print "Usage: %s cmd arguments..." % os.path.basename(sys.argv[0])
        print "Available commands:"
        for prg in prgs:
            print " %20s\t%s" % (prg,prgs[prg][1])        
    
    sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] in prgs:
        if len(sys.argv) > 2:
            prgs[sys.argv[1]][0](*sys.argv[2:])
        else:
            usage(sys.argv[1])
    else:
        usage()
