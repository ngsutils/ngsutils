
import sys
import os
import gzip
import math
import collections
from eta import ETA


FASTQRead = collections.namedtuple('FASTQRead', 'name seq qual')


class FASTQ(object):
    def __init__(self, fname=None, fileobj=None):
        self.fname = fname

        if fileobj:
            self.fileobj = fileobj
        elif fname:
            if fname == '-':
                self.fileobj = sys.stdin
            elif fname[-3:] == '.gz' or fname[-4:] == '.bgz':
                self.fileobj = gzip.open(os.path.expanduser(fname))
            else:
                self.fileobj = open(os.path.expanduser(fname))
        else:
            raise ValueError("Must pass either a fileobj or fname!")

    def seek(self, pos, whence=0):
        self.fileobj.seek(pos, whence)

    def fetch(self, quiet=False):
        if self.fname and not quiet:
            eta = ETA(os.stat(self.fname).st_size, fileobj=self.fileobj)
        else:
            eta = None

        while True:
            try:
                name = self.fileobj.next().strip()[1:]
                seq = self.fileobj.next().strip()
                self.fileobj.next()
                qual = self.fileobj.next().strip()

                if eta:
                    eta.print_status(name)
                yield FASTQRead(name, seq, qual)

            except:
                break

        if eta:
            eta.done()

    def close(self):
        if self.fileobj != sys.stdout:
            self.fileobj.close()

    def check_qualtype(self, num_to_check=10000):
        '''
        Checks a FASTQ file's quality score to see what encoding/scaling is used:
        Sanger, Solexa, or Illumina

        returns "Sanger", "Solexa", "Illumina", or "Unknown"
        '''

        # these are the differential values, unscaled from chr()
        sanger = (33, 73)
        solexa = (59, 104)
        illumina = (64, 104)

        sanger_count = 0
        solexa_count = 0
        illumina_count = 0
        unknown_count = 0

        checked = 0
        for read in self.fetch(quiet=True):
            if checked > num_to_check:
                break
            qmax = None
            qmin = None
            for q in [ord(x) for x in read.qual]:
                if qmin is None or q < qmin:
                    qmin = q
                if qmax is None or q > qmax:
                    qmax = q

            if sanger[0] <= qmin <= qmax <= sanger[1]:
                sanger_count += 1
            elif illumina[0] <= qmin <= qmax <= illumina[1]:
                illumina_count += 1
            elif solexa[0] <= qmin <= qmax <= solexa[1]:
                solexa_count += 1
            else:
                unknown_count += 1
            checked += 1

        self.seek(0)

        if unknown_count > 0:
            return 'Unknown'  # We don't have any idea about at least one of these reads

        if solexa_count > 0:
            # If there are any reads that fall in the Solexa range,
            # this must be a Solexa scale file. This should be rare.
            return 'Solexa'

        if sanger_count > illumina_count:
            return 'Sanger'
        return 'Illumina'

    def is_colorspace(self, num_to_check=10):
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

        for read in self.fetch(quiet=True):
            if seqs_checked > num_to_check:
                break
            checked = False
            for base in read.seq:
                if base in valid_colorspace:
                    is_colorspace = True
                    checked = True
                elif base in valid_basespace:
                    checked = True
            if checked:
                seqs_checked += 1

        return is_colorspace

    def is_paired(self, num_to_check=10):
        '''
        Determines if a FASTQ file has paired reads. This returns True is the file has
        paired reads with the same name in consecutive order.
        '''

        last_name = None
        last_count = 0
        count = 0

        frag_counts = []

        for read in self.fetch(quiet=True):
            name = read.name.split()[0]
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

            if count > num_to_check:
                return max(frag_counts) > 1


# def read_fastq(fname, quiet=False, eta_callback=None):
#     with ngsutils.support.ngs_utils.gzip_opener(fname) as f:
#         if fname == '-':
#             quiet = True
#         if not quiet:
#             eta = ETA(os.stat(fname).st_size, fileobj=f)
#         while f:
#             try:
#                 name = f.next().strip()
#                 seq = f.next().strip()
#                 f.next()
#                 qual = f.next().strip()

#                 if eta_callback:
#                     extra = eta_callback()
#                 else:
#                     extra = name
#                 if not quiet:
#                     eta.print_status(extra=extra)
#                 yield (name, seq, qual)
#             except:
#                 break
#     if not quiet:
#         eta.done()

def convert_illumina_qual(qual):
    '''
    Illumina char: QPhred + 64
    Phred char: QPhred + 33
    '''

    return ''.join([chr(ord(q) - 31) for q in qual])


def convert_solexa_qual(qual):
    '''
    Illumina char: QSolexa + 64  (note: this is for very old samples)
    Phred char: QPhred + 33

    QPhred = -10 * log10 (1/error)
    QSolexa = -10 * log10 (error/(1-error))

    QPhred = 10 * log10 (10 ^ (QSolexa/10) + 1)

    '''

    rv = []
    for q in qual:
        val = ord(q) - 64
        qp = 10 * math.log10(10 ** (val / 10) + 1)
        rv.append(chr(qp + 33))
    return ''.join(rv)
