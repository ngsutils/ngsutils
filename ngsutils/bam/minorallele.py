#!/usr/bin/env python
## category DNA-seq
## desc Find potential minor allele frequency
## experimental
"""
Calculates a minor allele frequency in (pooled) genomic sequencing.

Given a BAM file and a genomic reference, for each position covered in the
BAM file, show the reference base, the potential minor allele, and probable
background.

This assumes that if a SNP exists, there is likely only one possible variation.
So, this calculation will fail if more than one minor allele is present.  This
also ignores indels.

If R is installed and -alleles is given, this will also calculate a 95%
confidence interval. If rpy2 is installed, that bridge will be used to call R,
otherwise a process is forked for each allele.
"""

import os
import sys
import math
import subprocess
from ngsutils.bam import bam_pileup_iter
import pysam


class Blackhole(object):
    def write(self, string):
        pass

__sink = Blackhole()

try:
    __rsrc = os.path.join(os.path.dirname(__file__), 'support', 'minorallele_cpci.R')
    import rpy2.robjects as robjects
    rscript = True
    with open(__rsrc) as f:
        robjects.r(f.read())
except Exception:
    robjects = None
    rscript = os.path.join(os.path.dirname(__file__), 'support', 'minorallele_cpci.rsh')
    if not os.path.exists(rscript):
        sys.stderr.write('Missing R script: %s\n' % rscript)
        sys.exit(-1)

    stdout = sys.stdout
    sys.stdout = __sink
    retval = subprocess.Popen([rscript], stdout=subprocess.PIPE).wait()
    sys.stdout = stdout
    if retval != 0:
        sys.stderr.write('Error calling R script: %s\n' % rscript)
        rscript = None


def usage():
    print __doc__
    print """
Usage: bamutils minorallele {opts} in.bam ref.fa

Arguments:
  in.bam        BAM files to import
  ref.fa        Genomic reference sequence (indexed FASTA)

Options:
  -name    name  Sample name (default to filename)
  -qual    val   Minimum quality level to use in calculations
                 (numeric, Sanger scale) (default 0)
  -count   val   Only report bases with a minimum coverage of {val}
                 (default 0)
  -ci-low  val   Only report bases where the 95% CI low is greater than {val}
                 (default: show all)
  -alleles val   The number of alleles included in this sample
                 If given, a Clopper-Pearson style confidence interval will
                 be calculated. (requires rpy2 or R)
"""
    if robjects:
        print "rpy2 detected!"
    else:
        print "rpy2 not detected! Consider installing rpy2 for faster processing!"

    sys.exit(1)


def bam_minorallele(bam_fname, ref_fname, min_qual=0, min_count=0, num_alleles=0, name=None, min_ci_low=None):
    bam = pysam.Samfile(bam_fname, "rb")
    ref = pysam.Fastafile(ref_fname)

    if not name:
        name = os.path.basename(bam_fname)

    if num_alleles:
        print "# %s" % num_alleles

    sys.stdout.write('\t'.join("chrom pos refbase altbase total refcount altcount background refback altback".split()))
    if num_alleles and rscript:
        sys.stdout.write("\tci_low\tci_high\tallele_lowt\tallele_high")
    sys.stdout.write('\n')

    for pileup in bam_pileup_iter(bam, mask=1540):
        chrom = bam.getrname(pileup.tid)

        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        total = 0

        for pileupread in pileup.pileups:
            if not pileupread.is_del:
                if min_qual:
                    if pileupread.alignment.qual[pileupread.qpos] < min_qual:
                        continue
                if pileupread.indel == 0:
                    base = pileupread.alignment.seq[pileupread.qpos].upper()
                    if base != 'N':
                        counts[base] += 1
                        total += 1

        if total > min_count:
            refbase = ref.fetch(chrom, pileup.pos, pileup.pos + 1).upper()
            if not refbase in counts:
                continue

            refcount = counts[refbase]

            # sort non-ref counts.  first is alt, next is background

            scounts = []
            for c in counts:
                if c != refbase:
                    scounts.append((counts[c], c))

            scounts.sort()
            scounts.reverse()

            altbase = scounts[0][1]
            altcount = scounts[0][0]
            background = scounts[1][0]

            refback = refcount - background
            altback = altcount - background

            if (altback + refback) == 0:
                altfreq = 0
            else:
                altfreq = float(altback) / (altback + refback)

            cols = [chrom, pileup.pos + 1, refbase, altbase, total, refcount, altcount, background, refback, altback, altfreq]
            if num_alleles and rscript:
                ci_low, ci_high = calc_cp_ci(refback + altback, altback, num_alleles)
                allele_low = ci_low * num_alleles
                allele_high = ci_high * num_alleles

                cols.append(ci_low)
                cols.append(ci_high)
                cols.append(allele_low)
                cols.append(allele_high)
            else:
                ci_low = 0

            if not math.isnan(ci_low) and (min_ci_low is None or ci_low > min_ci_low):
                print '\t'.join([str(x) for x in cols])

    bam.close()
    ref.close()


__ci_cache = {}


def calc_cp_ci(N, count, num_alleles):
    if (N, count, num_alleles) in __ci_cache:
        return __ci_cache[(N, count, num_alleles)]

    vals = (float('nan'), float('nan'))

    if robjects:
        stdout = sys.stdout
        sys.stdout = __sink
        stderr = sys.stderr
        sys.stderr = __sink
        vals = robjects.r['CP.CI'](N, count, num_alleles)
        sys.stdout = stdout
        sys.stderr = stderr
    else:
        vals = [float(x) for x in subprocess.Popen([str(x) for x in [rscript, N, count, num_alleles]], stdout=subprocess.PIPE).communicate()[0].split()]

    __ci_cache[(N, count, num_alleles)] = vals
    return vals

if __name__ == '__main__':
    bam = None
    ref = None

    min_qual = 0
    min_count = 0
    min_ci = None
    num_alleles = 0
    name = None

    last = None
    for arg in sys.argv[1:]:
        if last == '-qual':
            min_qual = int(arg)
            last = None
        elif last == '-count':
            min_count = int(arg)
            last = None
        elif last == '-alleles':
            num_alleles = int(arg)
            last = None
        elif last == '-ci-low':
            min_ci = int(arg)
            last = None
        elif last == '-name':
            name = arg
            last = None
        elif arg == '-h':
            usage()
        elif arg in ['-qual', '-count', '-alleles', '-name', '-ci-low']:
            last = arg
        elif not bam and os.path.exists(arg) and os.path.exists('%s.bai' % arg):
            bam = arg
        elif not ref and os.path.exists(arg) and os.path.exists('%s.fai' % arg):
            ref = arg
        else:
            print "Unknown option or missing index: %s" % arg
            usage()

    if min_ci and not num_alleles:
        print "Can't specify a minimum CI level without specifying the number of alleles!"
        usage()

    if not bam or not ref:
        usage()

    bam_minorallele(bam, ref, min_qual, min_count, num_alleles, name, min_ci)
