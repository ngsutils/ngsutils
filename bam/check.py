#!/usr/bin/env python
## category Misc
## desc Checks a BAM file for corruption
'''
Checks a BAM file for corruption
'''
import sys
import os
import pysam


def bam_check(fname):
    sys.stdout.write('%s: ' % fname)
    sys.stdout.flush()
    fail = False
    try:
        bamfile = pysam.Samfile(fname, "rb")
        for read in bamfile:
            pass
        bamfile.close()
    except KeyboardInterrupt:
        sys.stdout.write('\n')
        sys.exit(-1)
    except:
        fail = True
        pass

    if fail:
        sys.stdout.write('ERROR\n')
        return False

    sys.stdout.write('OK\n')


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
