#!/usr/bin/env python
## category General
## desc Count the number of reads containing a certain tag
## experimental
"""
Counts the number of reads containing unique tag values.

For example, if reads were annotated with a certain tag value, this
will find the unique tag values and count how many reads matched the
unique tag.
"""

import os
import sys
import pysam
import ngsutils.bam


def usage(msg=None):
    if msg:
        print msg
    print __doc__
    print """
Usage: bamutils tallytags {opts} input.bam

Options
  -tag VAL         The primary tag to tally. The number of reads that match
                   each unique value for this tag will be counted.

  -secondary VAL   If a secondary tag is given, the average value of this will
                   be determined for each read that matches the primary. (In this
                   case, the primary tag is similar to the SQL GROUP BY clause).

  -fragment        If each read of a fragment have the same tag value, then if
                   this is set, both reads will be counted as one.
"""
    sys.exit(1)


class SecondaryCounter(object):
    def __init__(self):
        self.acc = 0
        self.count = 0

    def add(self, val):
        self.acc += val
        self.count += 1

    def mean(self):
        return float(self.acc) / self.count


class TagTally(object):
    def __init__(self, secondary_tag_names):
        self.counts = {}
        self.secondary_tag_names = secondary_tag_names

    def add(self, primary_tag, secondary_tags):
        if not primary_tag in self.counts:
            self.counts[primary_tag] = { 'count': 0 }
            for tag in self.secondary_tag_names:
                self.counts[primary_tag][tag] = SecondaryCounter()

        self.counts[primary_tag]['count'] += 1

        for tag, val in zip(self.secondary_tag_names, secondary_tags):
            self.counts[val][tag].add(val)


def bam_tallytag(infile, primary_tag, secondary_tags=None, fragment_count=False, quiet=False):
    counts = {}
    inbam = pysam.Samfile(infile, "rb")

    paircache = {}

    for read in ngsutils.bam.bam_iter(inbam, quiet=quiet):
        try:
            val = read.opt(primary_tag)
        except KeyError:
            continue

        if not val in counts:
            counts[val] = { 'count': 0 }
            for tag in secondary_tags:
                counts[val][tag] = SecondaryCounter()

        count = 1
        if fragment_count:
            if read.qname in paircache:
                if paircache[read.qname] != val:
                    counts[val]['count'] += 1

                del paircache[read.qname]
            else:
                counts[val]['count'] += 1
                tag_vals = []
                for tag in secondary_tags:
                    try:
                        tag_val = read.opt(tag)
                        tag_vals.append(tag_val)
                        counts[val][tag].add(tag_val)
                    except KeyError:
                        continue
                paircache[read.qname] = val
        else:
            counts[val]['count'] += 1
            for tag in secondary_tags:
                try:
                    counts[val][tag].add(read.opt(tag))
                except KeyError:
                    continue

    sys.stdout.write('%s\tcount' % primary_tag)
    for sec_tag in secondary_tags:
        sys.stdout.write('\t%s ave' % sec_tag)
    sys.stdout.write('\n')

    for pri_tag in counts:
        sys.stdout.write('%s\t%s' % (pri_tag, counts[pri_tag]['count']))
        for sec_tag in secondary_tags:
            sys.stdout.write('\t%s' % counts[pri_tag][sec_tag].mean())
        sys.stdout.write('\n')



if __name__ == '__main__':
    infile = None
    fragment_count = False
    primary_tag = None
    secondary_tags = []
    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-tag':
            primary_tag = arg
            last = None
        elif last == '-secondary':
            secondary_tags.append(arg)
            last = None
        elif arg in ['-tag', '-secondary']:
            last = arg
        elif arg == '-fragment':
            fragment_count = True
        elif not infile and os.path.exists(arg):
            infile = arg
        else:
            usage('Unknown option: %s' % arg)

    if not primary_tag or not infile:
        usage()

    bam_tallytag(infile, primary_tag, secondary_tags, fragment_count)
