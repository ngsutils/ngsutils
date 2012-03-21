
import os
import support.ngs_utils
import support.localalign
from support.eta import ETA


def read_fastq(fname, quiet=False, eta_callback=None):
    with support.ngs_utils.gzip_opener(fname) as f:
        if not quiet:
            eta = ETA(os.stat(fname).st_size, fileobj=f)
        while f:
            try:
                name = f.next().strip()
                seq = f.next().strip()
                f.next()
                qual = f.next().strip()

                if eta_callback:
                    extra = eta_callback()
                else:
                    extra = name
                if not quiet:
                    eta.print_status(extra=extra)
                yield (name, seq, qual)
            except:
                break
    if not quiet:
        eta.done()


def fastq_qualtype(fname, num_to_check=10000):
    '''
    Checks a FASTQ file's quality score to see what encoding/scaling is used:
    Sanger, Solexa, or Illumina
    '''

    # these are the differential values
    sanger = (33, 73)
    solexa = (59, 104)
    illumina = (64, 104)

    sanger_count = 0
    solexa_count = 0
    illumina_count = 0
    unknown_count = 0

    checked = 0
    for name, seq, qual in read_fastq(fname, quiet=True):
        if checked > num_to_check:
            break
        qmax = None
        qmin = None
        for q in qual:
            if qmin is None or ord(q) < qmin:
                qmin = ord(q)
            if qmax is None or ord(q) > qmax:
                qmax = ord(q)
        if sanger[0] <= qmin <= qmax <= sanger[1]:
            sanger_count += 1
        elif illumina[0] <= qmin <= qmax <= illumina[1]:
            illumina_count += 1
        elif solexa[0] <= qmin <= qmax <= solexa[1]:
            solexa_count += 1
        else:
            unknown_count += 1
        checked += 1

    totals = [(sanger_count, 'Sanger'), (solexa_count, 'Solexa'), (illumina_count, 'Illumina'), (unknown_count, 'Unknown')]
    totals.sort()
    return totals


def is_colorspace_fastq(fname):
    '''
    This works by scanning the first 10 reads that have sequences (aren't Ns
    or 4s). If there are any colorspace values, the entire file is called as
    colorspace.

    It's a bit overkill...
    '''

    seqs_checked = 0
    is_colorspace = False

    valid_basespace = "atcgATCG"
    valid_colorspace = "0123456"

    for name, seq, qual in read_fastq(fname, quiet=True):
        if seqs_checked > 10:
            break
        checked = False
        for base in seq:
            if base in valid_colorspace:
                is_colorspace = True
                checked = True
            elif base in valid_basespace:
                checked = True
        if checked:
            seqs_checked += 1

    return is_colorspace


def is_paired_fastq(fname):
    '''
    Determines if a FASTQ file has paired reads. This returns True is the file has
    paired reads with the same name consecutively.
    '''

    last_name = None
    last_count = 0
    count = 0

    frag_counts = []

    for name, seq, qual in read_fastq(fname, quiet=True):
        name = name.split()[0]
        if name != last_name:
            if last_count == 1:
                return 0
            else:
                if last_name:
                    frag_counts.append(last_count)
                last_count = 1
                last_name = name
        else:
            last_count += 1

        count += 1
        if count > 10:
            return max(frag_counts)
