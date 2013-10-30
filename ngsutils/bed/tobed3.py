#!/usr/bin/env python
## category Conversion
## desc Removes extra columns from a BED (or BED compatible) file
'''
Removes extra columns from a BED (or BED compatible) file.

'''

import os
import sys


def usage():
    print __doc__
    print """\
Usage: bedutils tobed3 bedfile

"""
    sys.exit(1)


def bed_tobed3(fname, out=sys.stdout):
    if fname == '-':
        f = sys.stdin
    else:
        f = open(fname)

    for line in f:
        cols = line.strip().split('\t')
        out.write('%s\t%s\t%s\n' % (cols[0], cols[1], cols[2]))

    if f != sys.stdin:
        f.close()

if __name__ == '__main__':
    fname = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg
        else:
            print "Unknown option: %s" % arg
            usage()

    if not fname:
        usage()

    bed_tobed3(fname)
