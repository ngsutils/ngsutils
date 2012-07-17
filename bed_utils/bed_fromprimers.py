#!/usr/bin/env python
## category Conversion
## desc Converts a list of PCR primer pairs to BED regions
'''
Converts a list of PCR primer pairs to BED regions

Given a genomic FASTA file and a list of primers,  this will generate a new
FASTA file for the targeted region.  This targeted FASTA file can then be
used for targeted resequencing.

PCR input is expected to be paired FASTA entries with names like:
>primername/1
atatcgtgctacgatc
>primername/2
ttgatcggcatataaa

Tab-delimited input should be:
fwd-primer {tab} rev-primer

Uses http://genome.ucsc.edu/cgi-bin/hgPcr to perform ePCR.

Note: This operates by screen scraping,  so it may break unexpectedly in the
future.
'''

import os
import sys
import urllib2
import re
import urllib
import time


def usage():
    print __doc__
    print """\
Usage: bedutils fromprimers {opts} -fasta file.fa
       bedutils fromprimers {opts} -tab file.txt
       bedutils fromprimers {opts} fwd_primer rev_primer

Options:
-db         name    DB to use
-perfect    bases   Minimum perfect bases (default: 15)
-good       bases   Minimum good bases (default: 15)
-size       bases   Max product size (default: 4000)
-flip               Flip the reverse primer (default: False)

"""
    sys.exit(1)


def insilico_pcr_tab(fasta,  db,  perfect=15,  good=15,  size=4000,  flip=False):
    pairs = set()
    with open(fasta) as f:
        for line in f:
            fwd,  rev = [x.strip().upper() for x in line.strip().split('\t')[0:2]]
            if not (fwd,  rev) in pairs:
                pairs.add((fwd,  rev))
                insilico_pcr_pairs(fwd,  rev,  len(pairs),  db,  perfect,  good,  size,  flip)


def insilico_pcr_fasta(tab,  db,  perfect=15,  good=15,  size=4000,  flip=False):
    primers = {}
    names = []
    with open(tab) as f:
        seq = ''
        name = None
        for line in f:
            if line[0] == '>':
                if seq:
                    if not name in primers:
                        primers[name] = [seq, ]
                        names.append(name)
                    else:
                        primers[name].append(seq)

                name = line[1:].split('/')[0]
                seq = ''
            else:
                seq += line.strip()

        if not name in primers:
            primers[name] = [seq, ]
        else:
            primers[name].append(seq)

    for name in names:
        insilico_pcr_pairs(primers[name][0],  primers[name][1],  name,  db,  perfect,  good,  size,  flip)


def insilico_pcr_pairs(fwd,  rev,  name,  db,  perfect=15,  good=15,  size=4000,  flip=False):
        params = [
            ('db', db),
            ('wp_perfect', perfect),
            ('wp_good', good),
            ('boolshad.wp_flipReverse', '0'),
            ('wp_f', fwd),
            ('wp_r', rev),
            ('Submit', 'submit'),
            ('wp_size', size)
            ]
        if flip:
            params.append(('wp_flipReverse',  'on'))

        attempt = 0
        result = None
        while attempt < 5 and not result:
            try:
                sys.stderr.write('>%s %s->%s (#%s) ' % (name,  fwd,  rev,  attempt + 1))
                f = urllib2.urlopen('http://genome.ucsc.edu/cgi-bin/hgPcr?%s' % urllib.urlencode(params),  timeout=180)
                if f:
                    result = f.read()
            except:
                pass

            if not result:
                attempt += 1
                sys.stderr.write("Timeout... waiting 30 sec\n")
                time.sleep(30)

        if result:
            count = 0
            for m in re.finditer('><A HREF=".*">(.*):(\d+)([\+\-])(\d+)</A>', result):
                if count == 0:
                    print "%s\t%s\t%s\t%s\t0\t%s" % (m.group(1), int(m.group(2)) - 1, m.group(4), name, m.group(3))
                else:
                    print "%s\t%s\t%s\t%s.%s\t0\t%s" % (m.group(1), int(m.group(2)) - 1, m.group(4), name, count, m.group(3))
                m = re.search('<PRE>><A HREF=".*">(.*):(\d+)([\+\-])(\d+)</A>', result)
                count += 1
            sys.stderr.write(" %s match%s\n" % (count,  'es' if count != 1 else ''))
        else:
            sys.stderr.write("Error - skipping\n")


if __name__ == '__main__':
    db = None
    perfect = 15
    good = 15
    size = 4000
    flip = False
    fasta = None
    tab = None
    fwd = None
    rev = None

    last = None
    for arg in sys.argv[1:]:
        if last == '-perfect':
            perfect = int(arg)
            last = None
        elif last == '-good':
            good = int(arg)
            last = None
        elif last == '-size':
            size = int(arg)
            last = None
        elif last == '-db':
            db = arg
            last = None
        elif last == '-fasta' and os.path.exists(arg):
            fasta = arg
            last = None
        elif last == '-tab' and os.path.exists(arg):
            tab = arg
            last = None
        elif arg in ['-db',  '-good',  '-perfect',  '-size',  '-fasta',  '-tab']:
            last = arg
        elif arg == '-h':
            usage()
        elif arg == '-flip':
            flip = True
        else:
            ok = True
            for base in arg:
                if base not in 'ATGCatgc':
                    ok = False
                    break
            if ok and not fwd:
                fwd = arg
            elif ok:
                rev = arg
            else:
                print "Unknown option: %s" % arg
                usage()
    if not db:
        usage()

    if not fasta and not tab and (not fwd and not rev):
        usage()

    if fasta:
        insilico_pcr_fasta(fasta, db, perfect, good, size, flip)
    elif tab:
        insilico_pcr_tab(tab, db, perfect, good, size, flip)
    else:
        insilico_pcr_pairs(fwd, rev, 'stdin_primers', db, perfect, good, size, flip)
