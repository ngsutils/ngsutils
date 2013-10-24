#!/usr/bin/env python
## category Misc
## desc Tag FASTA sequence names with a prefix or suffix
'''
Tag FASTA sequence names with a prefix or suffix

'''

import sys
import os
from eta import eta_open_iter

def tag_fasta(fname, prefix='', suffix=''):
    name = ''
    for line in eta_open_iter(fname, callback=lambda: name):
        if line[0] == '>':
            spl = line[1:].strip().split(' ',1)
            name = '%s%s%s' % (prefix, spl[0], suffix)
            if len(spl) > 1:
                comment = ' %s' % spl[1]
            else:
                comment = ''
            sys.stdout.write('>%s%s\n' % (name, comment))
        else:
            sys.stdout.write(line)


def usage(msg=None):
    if msg:
        sys.stdout.write('%s\n' % msg)
    sys.stdout.write(__doc__)
    sys.stdout.write('Usage: ngsutils tag_fasta {-prefix val} {-suffix val} filename.fa\n\n')
    sys.exit(1)


if __name__ == '__main__':
    fname = None
    prefix = ''
    suffix = ''
    last = None

    for arg in sys.argv[1:]:
        if last == '-prefix':
            prefix = arg
            last = None
        elif last == '-suffix':
            suffix = arg
            last = None
        elif arg in ['-prefix', '-suffix']:
            last = arg
        elif not fname:
            if os.path.exists(arg):
                fname = arg
            else:
                usage('%s missing!' % arg)

    if not fname or (not prefix and not suffix):
        usage()

    tag_fasta(fname, prefix, suffix)
