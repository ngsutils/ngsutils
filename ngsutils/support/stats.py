'''
various statistical tests and methods...
'''
import math


def median(vals):
    '''
    >>> median([1,2,3])
    2
    >>> median([1,2,3,4])
    2.5
    '''
    vals.sort()

    if len(vals) % 2 == 1:
        return vals[len(vals) / 2]
    else:
        a = vals[(len(vals) / 2) - 1]
        b = vals[(len(vals) / 2)]
        return float(a + b) / 2


def mean_stdev(l):
    '''
    >>> mean_stdev([1,2,2,2])
    (1.75, 0.5)
    >>> mean_stdev([2,2,2,2])
    (2.0, 0.0)
    '''

    acc = 0
    for el in l:
        acc += el

    mean = float(acc) / len(l)
    acc = 0
    for el in l:
        acc += (el - mean) ** 2

    if len(l) > 2:
        stdev = math.sqrt(float(acc) / (len(l) - 1))
    else:
        stdev = 0.0

    return (mean, stdev)


def counts_median(d):
    '''
    Calculate the median from counts stored in a dictionary
    >>> counts_median({ 1: 4, 2: 1, 3: 4 })
    2
    >>> counts_median({ 1: 4, 3: 4 })
    2

    '''
    count = 0
    for k in d:
        count += d[k]

    if count == 0:
        return 0

    acc = 0.0
    last = 0
    for k in sorted(d):
        if last:
            return (last + k) / 2
        acc += d[k]
        if acc / count == 0.5:
            last = k
        elif acc / count > 0.5:
            return k


def counts_mean_stdev(d):
    '''

    calc mean / stdev when data is stored as counts in a dictionary

    Ex:
        { 1: 4, 2: 1, 3: 4 } = [1, 1, 1, 1, 2, 3, 3, 3, 3]

    >>> counts_mean_stdev({ 1: 4, 2: 1, 3: 4 })
    (2.0, 1.0)

    '''

    acc = 0
    count = 0
    for k in d:
        acc += k * d[k]
        count += d[k]

    mean = float(acc) / count

    acc = 0
    for k in d:
        acc += (((k - mean) ** 2) * d[k])

    if count > 2:
        stdev = math.sqrt(float(acc) / (count - 1))
    else:
        stdev = 0.0

    return (mean, stdev)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
