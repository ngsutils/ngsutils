#!/usr/bin/env python
## category Misc
## desc Checks a BAM file for corruption
'''
Checks a BAM file for corruption
'''
import sys
import os

import ngsutils.bam


def bam_check(fname, quiet=False):
    if not quiet:
        sys.stdout.write('%s: ' % fname)
        sys.stdout.flush()
    fail = False
    try:
        bamfile = ngsutils.bam.bam_open(fname)
        for read in ngsutils.bam.bam_iter(bamfile):
            pass
        bamfile.close()
    except KeyboardInterrupt:
        if not quiet:
            sys.stdout.write('\n')
        sys.exit(-1)
    except:
        fail = True
        pass

    if fail:
        if not quiet:
            sys.stdout.write('ERROR\n')
        return False

    if not quiet:
        sys.stdout.write('OK\n')
    return True


def usage():
    print __doc__
    print "Usage: bamutils check bamfile..."
    sys.exit(-1)

if __name__ == "__main__":
    fnames = []

    for arg in sys.argv[1:]:
        if arg == "-h":
            usage()
        elif os.path.exists(arg):
            fnames.append(arg)
        else:
            sys.stderr.write("File: %s not found!\n" % arg)
            usage()

    if not fnames:
        usage()

    fail = False
    for f in fnames:
        if not bam_check(f):
            fail = True

    if fail:
        sys.exit(1)
