#!/usr/bin/env python
## category General
## desc Update read names with a suffix (for merging)
'''
Tags each read in a BAM file

Currently supported tags:
  -suffix
  -xs
  -orig-ref
  -orig-pos
  -orig-cigar
'''
import sys
import os
from ngsutils.bam import bam_iter, cigar_tostr
import pysam


class BamWriter(object):
    def __init__(self, outname, infname):
        self.outname = outname
        self.inbam = pysam.Samfile(infname, "rb")

    def run_chain(self, chain):
        tmp = os.path.join(os.path.dirname(self.outname), '.tmp.%s' % os.path.basename(self.outname))
        outbam = pysam.Samfile(tmp, 'wb', template=self.inbam)

        for read in chain.filter(self.inbam):
            outbam.write(read)

        self.inbam.close()
        outbam.close()

        if os.path.exists(outfname):
            os.unlink(outfname)  # Not really needed on *nix
        os.rename(tmp, outfname)



class BamReader(object):
    def filter(self, bam):
        for read in bam_iter(bam):
            yield read


class Suffix(object):
    def __init__(self, parent, suffix):
        self.parent = parent
        self.suffix = suffix

    def filter(self, bam):
        for read in self.parent.filter(bam):
            read.qname = "%s%s" % (read.qname, self.suffix)
            yield read


class OrigRef(object):
    def __init__(self, parent, tag):
        self.parent = parent
        self.tag = tag

    def filter(self, bam):
        refs = list(bam.references)
        for read in self.parent.filter(bam):
            if not read.is_unmapped:
                read.tags = read.tags + [(self.tag, refs[read.tid])]
            else:
                read.tags = read.tags + [(self.tag, '*')]
            yield read


class OrigPos(object):
    def __init__(self, parent, tag):
        self.parent = parent
        self.tag = tag

    def filter(self, bam):
        for read in self.parent.filter(bam):
            if not read.is_unmapped:
                read.tags = read.tags + [(self.tag, read.pos)]
            else:
                read.tags = read.tags + [(self.tag, -1)]
            yield read


class OrigCIGAR(object):
    def __init__(self, parent, tag):
        self.parent = parent
        self.tag = tag

    def filter(self, bam):
        for read in self.parent.filter(bam):
            if not read.is_unmapped:
                read.tags = read.tags + [(self.tag, cigar_tostr(read.cigar))]
            else:
                read.tags = read.tags + [(self.tag, '*')]
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

    def filter(self, bam):
        for read in self.parent.filter(bam):
            read.tags = read.tags + [(self.key, self.value)]
            yield read


class CufflinksXS(object):
    def __init__(self, parent):
        self.parent = parent

    def filter(self, bam):
        for read in self.parent.filter(bam):
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
  -suffix suff     A suffix to add to each read name

  -xs              Add the XS:A tag for +/- strandedness (req'd by Cufflinks)

  -tag tag         Add an arbitrary tag (ex: -tag XX:Z:test)

  -orig-ref tag    Add a new tag with the original reference name (For
                   example, in a region-based BAM will be converted to
                   standard coordinates)

  -orig-pos tag    Add a new tag with the original reference pos

  -orig-cigar tag  Add a new tag with the original CIGAR alignment

  -f               Force overwriting the output BAM file if it exists
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
        elif last == '-orig-ref':
            args.append([OrigRef, arg])
            last = None
        elif last == '-orig-pos':
            args.append([OrigPos, arg])
            last = None
        elif last == '-orig-cigar':
            args.append([OrigCIGAR, arg])
            last = None
        elif arg in ['-suffix', '-tag', '-orig-ref', '-orig-pos', '-orig-cigar']:
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

    writer = BamWriter(outfname, infname)

    chain = BamReader()
    for arg in args:
        chain = arg[0](chain, *arg[1:])

    writer.run_chain(chain)
