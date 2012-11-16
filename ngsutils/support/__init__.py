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
