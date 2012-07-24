'''
various statistical tests and methods...
'''

import sys
import math

# RPy version of fisher test - call out to R
# def _fisher_test(vals):
#     '''
#     vals is a list: [sample1_pres,sample1_notpres,sample2_pres,sample2_notpres]
#     Computes a Fisher Exact test
#
#     The actual test is performed in R
#
#                      Group1   Group2    total
#                      ------   ------    -----
#        Present         a        b        a+b
#        Not Present     c        d        c+d
#                        a+c      b+d      n
#
#                (a+b)!(c+d)!(a+c)!(b+d)!
#        p   =   -----------------------
#                    n!a!b!c!d!
#     '''
#     if robjects:
#         table=robjects.r.matrix(robjects.FloatVector([vals[0],vals[2],vals[1],vals[3]]),nr=2)
#         p=robjects.r['fisher.test'](table)[0][0]
#         return p
#     return None
#


def fisher_test(a, b, c, d):
    '''
    vals is a list: [sample1_pres,sample1_notpres,sample2_pres,sample2_notpres]
    Computes a Fisher Exact test

    The actual test is performed in R

                     Group1   Group2    total
                     ------   ------    -----
       Present         a        b        a+b
       Not Present     c        d        c+d
                       a+c      b+d      n

               (a+b)!(c+d)!(a+c)!(b+d)!
       p   =   -----------------------
                   n!a!b!c!d!

    Note: This is a very simple calculation and will only work on small counts.
    For larger counts, use scipy.stats.fisher_exact
    '''
    n = a + b + c + d
    ab = math.factorial(a + b)
    cd = math.factorial(c + d)
    ac = math.factorial(a + c)
    bd = math.factorial(b + d)
    nf = math.factorial(n)
    af = math.factorial(a)
    bf = math.factorial(b)
    cf = math.factorial(c)
    df = math.factorial(d)
    try:
        p = float(ab * cd * ac * bd) / (nf * af * bf * cf * df)
    except Exception, e:
        sys.stderr.write('Error calculating p-value: a,b,c,d -> (%s, %s, %s, %s)\n' % (a, b, c, d))
        sys.stderr.write(str(e))
        sys.stderr.write('\n')
        return 1.0

    return p


def mean_stdev(l):
    acc = 0
    for el in l:
        acc += el

    mean = float(acc) / len(l)
    acc = 0
    for el in l:
        acc += (el - mean) ** 2

    stdev = math.sqrt(float(acc) / (len(l) - 1))
    return (mean, stdev)
