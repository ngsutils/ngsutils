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

    stdev = math.sqrt(float(acc) / (len(l) - 1))
    return (mean, stdev)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
