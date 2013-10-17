#!/usr/bin/env python
## category Misc
## desc Strips sequences from reference FASTA files
'''
Strips alternative assemblies from reference FASTA files

This removes all sequences by looking for specific substrings in their name
(defaults to '_')

'''

import sys
import os
from eta import eta_open_iter

def strip_fasta(fname, substr='_'):
    good = False
    name = ''
    rmlist = []
    def callback():
        if good:
            return '%s (removed %s)' % (name, len(rmlist))
        else:
            return '%s [x] (removed %s)' % (name, len(rmlist))

    for line in eta_open_iter(fname, callback=callback):
        if line[0] == '>':
            name = line[1:].strip().split(' ')[0]
            if substr in name:
                good = False
                rmlist.append(name)
            else:
                good = True
        if good:
            sys.stdout.write(line)

    sys.stderr.write('Removed:\n')
    for n in rmlist:
        sys.stderr.write('%s\n' % n)


def usage(msg=None):
    if msg:
        sys.stdout.write('%s\n' % msg)
    sys.stdout.write(__doc__)
    sys.stdout.write('Usage: ngsutils strip_fasta filename.fa [string-to-search-for]\n\n')
    sys.exit(1)


if __name__ == '__main__':
    fname = None
    substr = '_'
    for arg in sys.argv[1:]:
        if not fname:
            if os.path.exists(arg):
                fname = arg
            else:
                usage('%s missing!' % arg)
        else:
            substr = arg

    if not fname:
        usage()

    strip_fasta(fname, substr)