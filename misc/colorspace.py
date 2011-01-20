#!/usr/bin/env python
'''
Convert colorspace sequences
'''

import sys,os

def usage():
    print __doc__
    print """
Usage: %s [encode|decode] [type] {opts} filename

Where type is:
-fa     FASTA
-fq     FASTQ

Options:
decode 
    -pseudo    Don't to di-base conversions, just convert 0=A, 1=C, 2=G, 3=T.  
               This is useful for using tools that expect straight base-space
               sequences.
encode
   -prefix N   The prefix/primer base to include in the sequence
   -pseudo     Use pseudo-conversions (as above)

Encode will encode a base-space file into color-space
Decode will decode a color-space file into base-space

""" % os.path.basename(sys.argv[0])
    sys.exit(1)

__encode = {
    'AA':'0',
    'AC':'1',
    'AG':'2',
    'AT':'3',
    'A.':'4',
    '.A':'4',
    'AN':'4',
    'CA':'1',
    'CC':'0',
    'CG':'3',
    'CT':'2',
    'C.':'4',
    '.C':'4',
    'CN':'4',
    'GA':'2',
    'GC':'3',
    'GG':'0',
    'GT':'1',
    'G.':'4',
    '.G':'4',
    'GN':'4',
    'TA':'3',
    'TC':'2',
    'TG':'1',
    'TT':'0',
    'T.':'4',
    '.T':'4',
    'TN':'4',
    'N.':'4',
    '.N':'4',
    'NA':'4',
    'NC':'4',
    'NG':'4',
    'NT':'4',
    'NN':'4',
    'N.':'4',
    '..':'4'
    # 'NA':5,
    # 'NC':5,
    # 'NG':5,
    # 'NT':5,
    # 'NN':6,
    # 'N.':6,
    # '..':6
}

__decode = {
    'A0': 'A',
    'A1': 'C',
    'A2': 'G',
    'A3': 'T',
    'A4': 'N',
    'A.': 'N',
    'C0': 'C',
    'C1': 'A',
    'C2': 'T',
    'C3': 'G',
    'C4': 'N',
    'C.': 'N',
    'G0': 'G',
    'G1': 'T',
    'G2': 'A',
    'G3': 'C',
    'G4': 'N',
    'G.': 'N',
    'T0': 'T',
    'T1': 'G',
    'T2': 'C',
    'T3': 'A',
    'T4': 'N',
    'T.': 'N',
    'N0': 'N',
    'N1': 'N',
    'N2': 'N',
    'N3': 'N',
    'N4': 'N',
    'N.': 'N',
}


__pseudo_encode = {
    'A' : '0',
    'C' : '1',
    'G' : '2',
    'T' : '3',
}

__pseudo_decode = {
    '0':'A',
    '1':'C',
    '2':'G',
    '3':'T',
}

def encode_cs_seq(seq,last=None,pseudo=False):
    ret = []
    seq = seq.upper()
    if last:
        last = last.upper()
    if pseudo:
        if last:
            ret.append(last)
        for base in seq:
            ret.append(__pseudo_encode[base])
    else:
        if last:
            ret.append(last)
            ret.append(__encode['%s%s' % (last,seq[0])])
        else:
            ret.append(seq[0])
        last = seq[0]
        for base in seq[1:]:
            ret.append(__encode['%s%s' % (last,base)])
            last = base
    return ''.join(ret)

def decode_cs_seq(seq,last=None,pseudo=False):
    ret = []
    if last:
        last = last.upper()

    if pseudo:
        for base in seq:
            if base in __pseudo_decode:
                ret.append(__pseudo_decode[base])
    else:
        if last:
            ret.append(__decode['%s%s' % (last,seq[0])])
            last = ret[-1]
        else:
            last = seq[0]
        for base in seq[1:]:
            ret.append(__decode['%s%s' % (last,base)])
            last = ret[-1]
    return ''.join(ret)


def encode_cs_fasta(fname,pseudo=False,primer=''):
    with open(fname) as f:
        seq = ''
        for line in f:
            if line[0] == '>':
                if seq:
                    sys.stdout.write(encode_cs_seq(seq,primer,pseudo))
                    sys.stdout.write('\n')
                sys.stdout.write(line)
                seq = ''
            else:
                seq += line.strip()
        if seq:
            sys.stdout.write(encode_cs_seq(seq,primer,pseudo))
            sys.stdout.write('\n')

def decode_cs_fasta(fname,pseudo=False,primer=''):
    with open(fname) as f:
        seq = ''
        for line in f:
            if line[0] == '>':
                if seq:
                    sys.stdout.write(decode_cs_seq(seq,primer,pseudo))
                    sys.stdout.write('\n')
                sys.stdout.write(line)
                seq = ''
            else:
                seq += line.strip()
        if seq:
            sys.stdout.write(decode_cs_seq(seq,pseudo))
            sys.stdout.write('\n')

def encode_cs_fastq(fname,pseudo=False,primer=''):
    with open(fname) as f:
        while True:
            try:
                name = f.next()
                seq = f.next()
                plus = f.next()
                qual = f.next()
                sys.stdout.write('%s%s%s\n+\n%s' % (name,primer,encode_cs_seq(seq,primer,pseudo),qual))
            except:
                break

def decode_cs_fastq(fname,pseudo=False,primer=''):
    with open(fname) as f:
        while True:
            try:
                name = f.next()
                seq = f.next()
                plus = f.next()
                qual = f.next()
                sys.stdout.write('%s%s\n+\n%s' % (name,decode_cs_seq(seq,primer,pseudo),qual))
            except:
                break
                
if __name__ == '__main__':
    op = None
    ftype = None
    pseudo = False
    primer = None
    fname = None
    last = None
    
    for arg in sys.argv[1:]:
        if last == '-primer':
            primer = arg
        elif arg in ['encode', 'decode'] and not op:
            op = arg
        elif arg in ['-fq','-fa','-seq'] and not ftype:
            ftype = arg
        elif arg == '-pseudo':
            pseudo = True
        elif arg == '-primer':
            last = arg
        elif not fname and os.path.exists(arg):
            fname = arg
        else:
            sys.stderr.write('Unknown argument: %s\n' % arg)
            usage()
    if not op or not ftype or not fname:
        usage()
    
    if op == 'encode' and ftype == '-fa':
        encode_cs_fasta(fname,pseudo,primer)
    elif op == 'encode' and ftype == '-fq':
        encode_cs_fastq(fname,pseudo,primer)
    elif op == 'decode' and ftype == '-fa':
        decode_cs_fasta(fname,pseudo,primer)
    elif op == 'decode' and ftype == '-fq':
        decode_cs_fastq(fname,pseudo,primer)
        