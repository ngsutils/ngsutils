#!/usr/bin/env python
'''
Extracts N-mers from regions surrounding splice junctions (cassettes only at 
this point).

Splice junctions are given in MISO format: UPSTREAM@CASSETTE@DOWNSTREAM
each region is then in ref:start:end:strand format

Example:
chr12:54894009:54894121:-@chr12:54893047:54893189:-@chr12:54889666:54890619:-

Separate counts are given for the following regions:

    -------            -------            -------        
---|       |----------|       |----------|       |-------
    -------            -------            ------- 
       AAAA|BBB    CCC|DDDDDDD|EEE    FFF|GGGG

Up:
  A = min(exon_span or entire exon)
  B = intron_span

Cassette:
  C = intron_span
  D = entire exon
  E = intron_span

Down:
  F = intron_span
  G = min(exon_span or entire exon)

'''

import sys,os,itertools
import pysam
from support.eta import ETA

def combinations_with_replacement(iterable, r):
    # copied from Python 2.7
    # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
    pool = tuple(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield tuple(pool[i] for i in indices)

### TODO: FIXME - Doesn't generate all nmers correctly...
def generate_nmers(size,seed=None):
    l = []
    for tup in combinations_with_replacement('ACGT',size):
        l.append(''.join(tup))
    return l

def output_mers(fileobj,name,nmers,known):
    fileobj.write(name)
    for nmer in known:
        fileobj.write('\t')
        if nmer in nmers:
            fileobj.write(str(nmers[nmer]))
        else:
            fileobj.write('0')

    fileobj.write('\n')

__revcomp = {
'A':'T',
'C':'G',
'G':'C',
'T':'A',
}
    
def revcomp(seq):
    return ''.join([__revcomp[x] for x in seq[::-1]])
            
def find_nmers(ref,chrom,start,end,strand,size):
    seq = ref.fetch(chrom,start,end).upper()

    if strand == '-':
        seq = revcomp(seq)
    
    nmers = {}
    
    while len(seq) > size:
        sub = seq[:size]
        if not sub in nmers:
            nmers[sub] = 1
        else:
            nmers[sub] += 1
        seq = seq[1:]
    
    return nmers
    
def extract_nmers(ref_fname,cassette_fname,output_template,nmer_size=6,intron_span=200,exon_span=400):
    ref = pysam.Fastafile(ref_fname)
    known_nmers = generate_nmers(nmer_size)
    
    outfiles = {}
    regions = ['up_exon','up_intron','alt_5_intron','alt_exon','alt_3_intron','down_intron','down_exon',]
    for reg in regions:
        outfiles[reg] = open('%s.%s.txt' % (output_template,reg),'w')
        outfiles[reg].write('event')
        for nmer in known_nmers:
            outfiles[reg].write('\t%s' % nmer)
        outfiles[reg].write('\n')
    
    with open(cassette_fname) as f:
        eta = ETA(os.stat(cassette_fname).st_size,fileobj=f)
        f.next()
        for line in f:
            event = line.split('\t')[0]
            eta.print_status(extra=event)
            up,alt,down = [x.split(':') for x in event.split('@')]
            
            up_start = int(up[1])
            up_end = int(up[2])

            alt_start = int(alt[1])
            alt_end = int(alt[2])

            down_start = int(down[1])
            down_end = int(down[2])
            
            if up[3] == '+':
                # forward strand
                
                #up_exon
                if (up_end - up_start) > exon_span:
                    up_start = up_end - exon_span
                output_mers(outfiles['up_exon'],event,find_nmers(ref,up[0],up_start,up_end,up[3],nmer_size),known_nmers)

                #up_intron
                output_mers(outfiles['up_intron'],event,find_nmers(ref,up[0],up_end,up_end+intron_span,up[3],nmer_size),known_nmers)

                #alt
                output_mers(outfiles['alt_5_intron'],event,find_nmers(ref,alt[0],alt_start-intron_span,alt_start,alt[3],nmer_size),known_nmers)
                output_mers(outfiles['alt_exon'],event,find_nmers(ref,alt[0],alt_start,alt_end,alt[3],nmer_size),known_nmers)
                output_mers(outfiles['alt_3_intron'],event,find_nmers(ref,alt[0],alt_end,alt_end+intron_span,alt[3],nmer_size),known_nmers)

                #down_intron
                output_mers(outfiles['down_intron'],event,find_nmers(ref,down[0],down_start-intron_span,down_start,down[3],nmer_size),known_nmers)

                #down_intron
                if (down_end - down_start) > exon_span:
                    down_end = down_start + exon_span
                
                output_mers(outfiles['down_exon'],event,find_nmers(ref,down[0],down_start,down_end,down[3],nmer_size),known_nmers)
            else:
                # reverse strand
                # MISO event format is strange.  It has the order for events correct (5->3, with strand correct)
                # but the start/end of the regions is relative to the + strand always (start < end)

                #up_exon
                if (up_end - up_start) > exon_span:
                    up_end = up_start + exon_span
                output_mers(outfiles['up_exon'],event,find_nmers(ref,up[0],up_start,up_end,up[3],nmer_size),known_nmers)

                #up_intron
                output_mers(outfiles['up_intron'],event,find_nmers(ref,up[0],up_start-intron_span,up_start,up[3],nmer_size),known_nmers)

                #alt (switches 5/3)
                output_mers(outfiles['alt_3_intron'],event,find_nmers(ref,alt[0],alt_start-intron_span,alt_start,alt[3],nmer_size),known_nmers)
                output_mers(outfiles['alt_exon'],event,find_nmers(ref,alt[0],alt_start,alt_end,alt[3],nmer_size),known_nmers)
                output_mers(outfiles['alt_5_intron'],event,find_nmers(ref,alt[0],alt_end,alt_end+intron_span,alt[3],nmer_size),known_nmers)

                #down_intron
                output_mers(outfiles['down_intron'],event,find_nmers(ref,down[0],down_end,down_end+intron_span,down[3],nmer_size),known_nmers)

                #down_intron
                if (down_end - down_start) > exon_span:
                    down_start = down_end - exon_span
                
                output_mers(outfiles['down_exon'],event,find_nmers(ref,down[0],down_start,down_end,down[3],nmer_size),known_nmers)
        eta.done()
            


def test():
    print generate_nmers(1)
    print generate_nmers(2)
    print generate_nmers(3)
    print len(generate_nmers(4))
    print len(generate_nmers(5))
    print len(generate_nmers(6))

def usage():
    print __doc__
    print """Usage: %s {opts} reference.fa cassettes.txt output_template

Required:
  reference.fa    - indexed FASTA file
  cassettes.txt   - tab-delimited file with cassettes in first column 
                    (ex: MISO bayes-factor format)
  output_template - template for output filenames
                    The following files will be created:
                        template.up_exon.txt
                        template.up_intron.txt
                        template.alt_5_intron.txt
                        template.alt_exon.txt
                        template.alt_3_intron.txt
                        template.down_intron.txt
                        template.down_exon.txt
Options:
  -intron_span  N   How far into the introns to extract 
                    (default: 200)
  -exon_span    N   How far into the up/down exons to extract 
                    (default: 400)
  -size         N   Size of the fragments to extract 
                    (default: 6)
"""
    sys.exit(1)

if __name__ == '__main__':
    ref_fname = None
    cassette_fname = None
    out_template = None
    intron_span = 200
    exon_span = 400
    nmer_size = 6
    
    last = None
    
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if arg == 'test':
            test()
            sys.exit()
        if last == '-intron_span':
            intron_span = int(arg)
            last = None
        elif last == '-exon_span':
            exon_span = int(arg)
            last = None
        elif last == '-size':
            nmer_size = int(arg)
            last = None
        elif arg in ['-intron_span','-exon_span','-size']:
            last = arg
        elif not ref_fname and os.path.exists(arg) and os.path.exists('%s.fai' % arg):
            ref_fname = arg
        elif not cassette_fname and os.path.exists(arg):
            cassette_fname = arg
        elif not out_template:
            out_template = arg
        else:
            print "Unknown option: %s" % arg
            usage()
    
    if not ref_fname or not cassette_fname or not out_template:
        usage()
    
    
    
    extract_nmers(ref_fname,cassette_fname,out_template,nmer_size,intron_span,exon_span)