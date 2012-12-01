import collections
import gzip
import os
import sys
import re
from eta import ETA


class FASTARecord(collections.namedtuple('FASTARecord', 'name comment seq')):
    def __repr__(self):
        if self.comment:
            return '>%s %s\n%s\n' % (self.name, self.comment, self.seq)
        return '>%s\n%s\n' % (self.name, self.seq)


class FASTAFile(object):
    def __init__(self, fname=None, fileobj=None):
        self.fname = fname
        if fileobj:
            self.fileobj = fileobj
        else:
            if self.fname == '-':
                self.fileobj = sys.stdin
            elif self.fname[-3:] == '.gz' or self.fname[-4:] == '.bgz':
                self.fileobj = gzip.open(os.path.expanduser(self.fname))
            else:
                self.fileobj = open(os.path.expanduser(self.fname))

        if not self.fileobj:
            raise ValueError("Missing valid filename or fileobj")

    def close(self):
        if self.fileobj != sys.stdout:
            self.fileobj.close()

    def read(self, quiet=False):
        name = ''
        comment = ''
        seq = ''

        if not quiet and self.fname and self.fname != '-':
            eta = ETA(os.stat(self.fname).st_size, fileobj=self.fileobj)
        else:
            eta = None

        for line in self.fileobj:
            line = line.strip()
            if line[0] == '#':
                continue

            if line[0] == '>':
                if name and seq:
                    if eta:
                        eta.print_status(name)
                    yield FASTARecord(name, comment, seq)

                spl = re.split(r'[ \t]', line[1:], maxsplit=1)
                name = spl[0]
                if len(spl) > 1:
                    comment = spl[1]
                else:
                    comment = ''
                seq = ''

            else:
                seq += line

        if name and seq:
            if eta:
                eta.print_status(name)
            yield FASTARecord(name, comment, seq)

        if eta:
            eta.done()


def gzip_reader(fname, quiet=False):
    if fname == '-':
        f = sys.stdin
    elif fname[-3:] == '.gz' or fname[-4:] == '.bgz':
        f = gzip.open(os.path.expanduser(fname))
    else:
        f = open(os.path.expanduser(fname))

    if quiet or fname == '-':
        eta = None
    else:
        eta = ETA(os.stat(fname).st_size, fileobj=f)

    for line in f:
        if eta:
            eta.print_status()
        yield line

    if f != sys.stdout:
        f.close()

    if eta:
        eta.done()


class Symbolize(object):
    'Converts strings to symbols - basically a cache of strings'
    def __init__(self):
        self.__cache = {}

    def __getitem__(self, k):
        if not k in self.__cache:
            self.__cache[k] = k

        return self.__cache[k]

symbols = Symbolize()

_compliments = {
'a': 't',
'A': 'T',
'c': 'g',
'C': 'G',
'g': 'c',
'G': 'C',
't': 'a',
'T': 'A',
'n': 'n',
'N': 'N'
}


def revcomp(seq):
    '''
    >>> revcomp('ATCGatcg')
    'cgatCGAT'
    '''
    ret = []

    for s in seq:
        ret.append(_compliments[s])

    ret.reverse()
    return ''.join(ret)


class Counts(object):
    '''
    Setup simple binning.  Bins are continuous 0->max.  Values are added to
    bins and then means / distributions can be calculated.
    '''
    def __init__(self):
        self.bins = []

    def add(self, val):
        while len(self.bins) <= val:
            self.bins.append(0)
        self.bins[val] += 1

    def mean(self):
        acc = 0
        count = 0

        for i, val in enumerate(self.bins):
            acc += (i * val)
            count += val

        return float(acc) / count

    def max(self):
        return len(self.bins) - 1
