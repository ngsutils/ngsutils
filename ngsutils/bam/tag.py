#!/usr/bin/env python
## category General
## desc Update read names with a suffix (for merging)
'''
Tags each read in a BAM file

Currently supported tags:
  -suffix
  -xs
'''
import sys
import os
from ngsutils.bam import bam_iter
import pysam


class BamWriter(object):
    def __init__(self, outname, template):
        self.outname = outname
        self.template = template

    def run_chain(self, chain):
        tmp = os.path.join(os.path.dirname(self.outname), '.tmp.%s' % os.path.basename(self.outname))
        outbam = pysam.Samfile(tmp, 'wb', template=self.template)

        for read in chain.filter():
            outbam.write(read)

        outbam.close()

        if os.path.exists(outfname):
            os.unlink(outfname)  # Not really needed on *nix
        os.rename(tmp, outfname)


class BamReader(object):
    def __init__(self, fname):
        self.bamfile = pysam.Samfile(fname, "rb")

    def filter(self):
        for read in bam_iter(self.bamfile):
            yield read

        self.bamfile.close()


class Suffix(object):
    def __init__(self, parent, suffix):
        self.parent = parent
        self.suffix = suffix

    def filter(self):
        for read in self.parent.filter():
            read.qname = "%s%s" % (read.qname, self.suffix)
            yield read


class Tag(object):
    def __init__(self, parent, tag):
        self.parent = parent

        spl = tag.rsplit(':', 1)
        self.key = spl[0]

        if self.key[-2:].lower() == ':i':
            self.value = int(spl[1])
        elif self.key[-2:].lower() == ':f':
            self.value = float(spl[1])
        else:
            self.value = spl[1]

        if ':' in self.key:
            self.key = self.key.split(':')[0]

    def filter(self):
        for read in self.parent.filter():
            read.tags = read.tags + [(self.key, self.value)]
            yield read


class CufflinksXS(object):
    def __init__(self, parent):
        self.parent = parent

    def filter(self):
        for read in self.parent.filter():
            if not read.is_unmapped:
                if read.is_reverse:
                    read.tags = read.tags + [('XS', '-')]
                else:
                    read.tags = read.tags + [('XS', '+')]

            yield read


def usage():
    print __doc__
    print """\
Usage: bamutils tag {opts} in.bamfile out.bamfile

Arguments:
  in.bamfile    The input BAM file
  out.bamfile   The name of the new output BAM file

Options:
  -suffix suff  A suffix to add to each read name
  -xs           Add the XS:A tag for +/- strandedness (req'd by Cufflinks)
  -tag tag      Add an arbitrary tag (ex: -tag XX:Z:test)
  -f            Force overwriting the output BAM file if it exists
"""
    sys.exit(1)

if __name__ == "__main__":
    infname = None
    outfname = None
    force = False
    last = None

    args = []

    for arg in sys.argv[1:]:
        if arg == '-f':
            force = True
        elif last == '-suffix':
            args.append([Suffix, arg])
            last = None
        elif last == '-tag':
            args.append([Tag, arg])
            last = None
        elif arg in ['-suffix', '-tag']:
            last = arg
        elif arg == '-xs':
            args.append([CufflinksXS, ])
        elif not infname and os.path.exists(arg):
            infname = arg
        elif not outfname:
            outfname = arg

    if not infname or not outfname or not args:
        usage()

    if not force and os.path.exists(outfname):
        sys.stderr.write('ERROR: %s already exists! Not overwriting without force (-f)\n\n' % outfname)
        sys.exit(1)

    chain = BamReader(infname)
    writer = BamWriter(outfname, chain.bamfile)

    for arg in args:
        chain = arg[0](chain, *arg[1:])

    writer.run_chain(chain)
