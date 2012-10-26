import sys
import os
import pysam
from ngsutils.support.eta import ETA


def bam_open(fname, mode='r'):
    if fname.lower()[-4:] == '.bam':
        return pysam.Samfile(fname, '%sb' % mode)
    return pysam.Samfile(fname, '%s' % mode)


def bam_iter(bam, quiet=False, show_ref_pos=False):
    if not quiet:
        eta = ETA(os.stat(bam.filename).st_size)

    if os.path.exists('%s.bai' % bam.filename):
        # This is an indexed file, so it is ref sorted...
        # Meaning that we should show chrom:pos, instead of read names
        show_ref_pos = True

    for read in bam:
        pos = bam.tell()
        bgz_offset = pos >> 16

        if not quiet:
            if (show_ref_pos):
                eta.print_status(bgz_offset, extra='%s:%s %s' % (bam.getrname(read.tid), read.pos, read.qname))
            else:
                eta.print_status(bgz_offset, extra='%s' % read.qname)
        yield read

    if not quiet:
        eta.done()


def read_calc_mismatches(read):
    inserts = 0
    deletions = 0
    indels = 0
    edits = int(read.opt('NM'))
    #
    # NM counts the length of indels
    # We really just care about *if* there is an indel, not the size
    #

    for op, length in read.cigar:
        if op == 1:
            inserts += length
            indels += 1
        elif op == 2:
            deletions += length
            indels += 1

    return edits - inserts - deletions + indels


def _extract_md_matches(md, maxlength):
    md_pos = 0

    while md and md_pos < maxlength:
        tmp = '0'  # preload a zero so that immediate mismatches will be caught
                   # the zero will have no affect otherwise...

        # look for matches
        while md and md[0] in '0123456789':
            tmp += md[0]
            md = md[1:]

        pos = int(tmp)
        if pos > maxlength:
            return (maxlength, '%s%s' % (pos - maxlength, md))
        return (pos, md)


def read_calc_variations(read):
    'see _read_calc_variations'
    for tup in _read_calc_variations(read.pos, read.cigar, read.opt('MD'), read.seq):
        yield tup


def _read_calc_variations(start_pos, cigar, md, seq):
    '''
    For each variation, outputs a tuple: (op, pos, seq)

    op  - operation (0 = mismatch, 1 = insert, 2 = deletion) (like CIGAR)
    pos - 0-based position of the variation (relative to reference)
    seq - the base (or bases) involved in the variation
          for mismatch or insert, this is the sequence inserted
          for deletions, this is the reference sequence that was removed

    MD is the mismatch string. Not all aligners include the tag. If your aligner
    doesn't include this, then you'll need to add it, or use a different function
    (see: read_calc_mismatches_gen).

    Special care must be used to handle RNAseq reads that cross
    an exon-exon junction.

    Also: MD is a *really* dumb format that can't be read correctly with
          a regex. It must be processed in concert with the CIGAR alignment
          in order to catch all edge cases. Some implementations insert 0's
          at the end of inserts / deltions / variations to make parsing easier
          but not everyone follows this. Look at the complex examples: the
          CIGAR alignment may show an insert, but the MD just shows all matches.

    Examples: See: http://davetang.org/muse/2011/01/28/perl-and-sam/
              Also from CCBB actual mappings and manual altered (shortened,
              made more complex)
              (doctests included)

    Match/mismatch
    CIGAR: 36M
    MD:Z:  1A0C0C0C1T0C0T27
    MD:Z:  1ACCC1TCT27 (alternative)
                   1         2
          123456789012345678901234567890123456
    ref:  CGATACGGGGACATCCGGCCTGCTCCTTCTCACATG
           XXXX XXX
    read: CACCCCTCTGACATCCGGCCTGCTCCTTCTCACATG
          MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
          -ACCC-TCT---------------------------
    >>> list(_read_calc_variations(1, [(0,36)], '1A0C0C0C1T0C0T27', 'CACCCCTCTGACATCCGGCCTGCTCCTTCTCACATG'))
    [(0, 2, 'A'), (0, 3, 'C'), (0, 4, 'C'), (0, 5, 'C'), (0, 7, 'T'), (0, 8, 'C'), (0, 9, 'T')]

    Insert
    CIGAR: 6M1I29M
    MD:Z: 0C1C0C1C0T0C27
          C1CC1CTC27 (alt)
                    1         2
          123456^789012345678901234567890123456
    ref:  CACCCC^TCTGACATCCGGCCTGCTCCTTCTCACAT
          X XX X|XX
    read: GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT
          MMMMMMIMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
          G-GA-GGGG---------------------------
    >>> list(_read_calc_variations(1, [(0,6), (1,1), (0, 29)], '0C1C0C1C0T0C27', 'GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT'))
    [(0, 1, 'G'), (0, 3, 'G'), (0, 4, 'A'), (0, 6, 'G'), (1, 7, 'G'), (0, 7, 'G'), (0, 8, 'G')]
    >>> list(_read_calc_variations(1, [(0,6), (1,1), (0, 29)], 'C1CC1CTC27', 'GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT'))
    [(0, 1, 'G'), (0, 3, 'G'), (0, 4, 'A'), (0, 6, 'G'), (1, 7, 'G'), (0, 7, 'G'), (0, 8, 'G')]


    Deletion
    CIGAR: 9M9D27M
    MD:Z: 2G0A5^ATGATGTCA27
          2GA5^ATGATGTCA27 (alt)
    ref:  AGGAATGGGATGATGTCAGGGGTTCCAGGTGGAGACGAGGACTCC
            XX     ^^^^^^^^^
    read: AGTGATGGG^^^^^^^^^GGGGTTCCAGGTGGAGACGAGGACTCC
          MMMMMMMMMDDDDDDDDDMMMMMMMMMMMMMMMMMMMMMMMMMMM
          --TG-----ATGATGTCA---------------------------
    >>> list(_read_calc_variations(1, [(0,9), (2,9), (0, 27)], '2G0A5^ATGATGTCA27', 'AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC'))
    [(0, 3, 'T'), (0, 4, 'G'), (2, 10, 'ATGATGTCA')]


    Complex
    CIGAR: 9M9D11M1I15M
    MD:Z: 2G0A5^ATGATGTCAA26
    MD:Z: 2G0A5^ATGATGTCA0G26 (alt)
                   1         2         3         4
    pos:  123456789012345678901234567890123456789012345
    ref:  AGGAATGGGATGATGTCAGGGGTTCCAGG^GGAGACGAGGACTCC
            XX     ^^^^^^^^^X          |
    read: AGTGATGGG^^^^^^^^^AGGGTTCCAGGTGGAGACGAGGACTCC
          MMMMMMMMMDDDDDDDDDMMMMMMMMMMMMMMMMMMMMMMMMMMM
          --TG-----ATGATGTCAG----------T---------------
    >>> list(_read_calc_variations(1, [(0,9), (2,9), (0,11), (1,1), (0,15)], '2G0A5^ATGATGTCAA26', 'AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC'))
    [(0, 3, 'T'), (0, 4, 'G'), (2, 10, 'ATGATGTCA'), (0, 19, 'G'), (1, 30, 'T')]


    Complex example - inserts aren't separately handled by MD, only visible in CIGAR
    CIGAR: 14M2D16M3I42M
    MD:Z:  14^TC58
                   1         2         3            4         5         6         7
    pos:  12345678901234567890123456789012^^^345678901234567890123456789012345678901234567
    ref:  caagtatcaccatgtcaggcatttttttcatt^^^tttgtagagagagaagacttgctatgttgcccaagctggcct
                        ^^                |||
    read: CAAGTATCACCATG^^AGGCATTTTTTTCATTTGGTTTGTAGAGAGAGAAGACTTGCTATGTTGCCCAAGCTGGCCT
          MMMMMMMMMMMMMMDDMMMMMMMMMMMMMMMMIIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
          --------------tc----------------TGG------------------------------------------
    >>> list(_read_calc_variations(1, [(0,14), (2,2), (0,16), (1,3), (0,42)], '14^TC58', 'CAAGTATCACCATGAGGCATTTTTTTCATTTGGTTTGTAGAGAGAGAAGACTTGCTATGTTGCCCAAGCTGGCCT'))
    [(2, 15, 'TC'), (1, 33, 'TGG')]


    Complex example 2:
    CIGAR: 41M3I10M1I5M1I2M2I10M
    MD:Z:  44C2C6T6T6
                   1         2         3         4            5             6
    pos:  12345678901234567890123456789012345678901^^^2345678901^23456^78^^9012345678
    ref:  AGGGTGGCGAGATCGATGACGGCATTGGCGATGGTGATCTT^^^GAGCCACATG^CGGTC^GC^^GGATCTCCAG
                                                   |||   X  X   |   X |  ||   X
    read: AGGGTGGCGAGATCGATGACGGCATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG
          MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMIIIMMMMMMMMMMIMMMMMIMMIIMMMMMMMMMM
          -----------------------------------------tta---A--T---c---A----gt---G------


    13M 28M 3I 10M 1I 5M 1I 2M 2I 10M
    >>> list(_read_calc_variations(1, [(0, 41), (1, 3), (0, 10), (1, 1), (0, 5), (1, 1), (0, 2), (1, 2), (0, 10)], '44C2C6T6T6', 'AGGGTGGCGAGATCGATGACGGCATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG'))
    [(1, 42, 'TTA'), (0, 45, 'A'), (0, 48, 'T'), (1, 52, 'C'), (0, 55, 'A'), (1, 57, 'G'), (1, 59, 'GT'), (0, 62, 'G')]


    Splice junction example:
    CIGAR: 62M100N13M
    MD:Z:  2T27C44
                                                                                 1      1
                   1         2         3         4         5         6           6      7
    pos:  12345678901234567890123456789012345678901234567890123456789012| [100] |3456789012345
    ref:  CCTCATGACCAGCTTGTTGAAGAGATCCGACATCAAGTGCCCACCTTGGCTCGTGGCTCTCA|-------|CTTGCTCCTGCTC
            X                           X
    read: CCGCATGACCAGCTTGTTGAAGAGATCCGATATCAAGTGCCCACCTTGGCTCGTGGCTCTCA|-------|CTTGCTCCTGCTC
          --G---------------------------T-----------------------------------------------------

    >>> list(_read_calc_variations(1, [(0,62), (4,100), (0,13)], '2T27C44', 'CCGCATGACCAGCTTGTTGAAGAGATCCGATATCAAGTGCCCACCTTGGCTCGTGGCTCTCACTTGCTCCTGCTC'))
    [(0, 3, 'G'), (0, 31, 'T')]


    Splice junction example 2:
    CIGAR: 13M100N28M3I10M1I5M1I2M2I10M
    MD:Z:  44C2C6T6T6
                                      1         1         1            1             1
                   1                  2         3         4            5             6
    pos:  1234567890123| [100] |4567890123456789012345678901^^^2345678901^23456^78^^9012345678
    ref:  AGGGTGGCGAGAT|-------|CGATGACGGCATTGGCGATGGTGATCTT^^^GAGCCACATG^CGGTC^GC^^GGATCTCCAG
                                                            |||   X  X   |   X |  ||   X
    read: AGGGTGGCGAGAT|-------|CGATGACGGCATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG
          MMMMMMMMMMMMM         MMMMMMMMMMMMMMMMMMMMMMMMMMMMIIIMMMMMMMMMMIMMMMMIMMIIMMMMMMMMMM
          -------------         ----------------------------tta---A--T---c---A----gt---G------

    13M 100N 28M 3I 10M 1I 5M 1I 2M 2I 10M
    >>> list(_read_calc_variations(1, [(0, 13), (3, 100), (0, 28), (1, 3), (0, 10), (1, 1), (0, 5), (1, 1), (0, 2), (1, 2), (0, 10)], '44C2C6T6T6', 'AGGGTGGCGAGATCGATGACGGCATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG'))
    [(1, 142, 'TTA'), (0, 145, 'A'), (0, 148, 'T'), (1, 152, 'C'), (0, 155, 'A'), (1, 157, 'G'), (1, 159, 'GT'), (0, 162, 'G')]


    Splice junction example 2A:
    CIGAR: 13M100N7M2D19M3I10M1I5M1I2M2I10M
    MD:Z:  9A10^GG22C2C6T6T6
                                      1         1         1            1             1
                   1                  2         3         4            5             6
    pos:  1234567890123| [100] |4567890123456789012345678901^^^2345678901^23456^78^^9012345678
    ref:  AGGGTGGCGAGAT|-------|CGATGACGGCATTGGCGATGGTGATCTT^^^GAGCCACATG^CGGTC^GC^^GGATCTCCAG
                                       ^^                   |||   X  X   |   X |  ||   X
    read: AGGGTGGCGCGAT|-------|CGATGAC^^CATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG
          MMMMMMMMMMMMM         MMMMMMMDDMMMMMMMMMMMMMMMMMMMIIIMMMMMMMMMMIMMMMMIMMIIMMMMMMMMMM
          ---------C---         ----------------------------tta---A--T---c---A----gt---G------
          .........A...         .......GG...................   ...C..C... ...T. ..  ...T......
              9    A        10        ^GG             22          C 2C   6   T     6   T   6

    >>> list(_read_calc_variations(1, [(0, 13), (3, 100), (0, 7), (2, 2), (0, 19), (1, 3), (0, 10), (1, 1), (0, 5), (1, 1), (0, 2), (1, 2), (0, 10)], '9A10^GG22C2C6T6T6', 'AGGGTGGCGCGATCGATGACCATTGGCGATGGTGATCTTTTAGAGACATATGCCGGACGGCGTGGAGCTCCAG'))
    [(0, 10, 'C'), (2, 121, 'GG'), (1, 142, 'TTA'), (0, 145, 'A'), (0, 148, 'T'), (1, 152, 'C'), (0, 155, 'A'), (1, 157, 'G'), (1, 159, 'GT'), (0, 162, 'G')]

    Real Example
    242_1071_1799_B1
    CIGAR: 42M10I3M1D9M1D11M
    MD:Z:  27G16A0^T6C2^T1C9
                   1         2         3         4                   5         6         7
    pos:  123456789012345678901234567890123456789012          345678901234567890123456789012345
    ref:  ACTGAGAAACCCAACCCTCTGAGACCAGCACACCCCTTTCAA^^^^^^^^^^GCATGTTCCTCCCTCCCCTTCTTTG
                                     X                          X^      X  ^ X
    read: ACTGAGAAACCCAACCCTCTGAGACCAACACACCCCTTTCAACACATTTTTGGCC^GTTCCTGCC^CGCCTTCTTTG
          MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMIIIIIIIIIIMMMDMMMMMMMMMDMMMMMMMMMMM
          ---------------------------A--------------^^^^^^^^^^--CT------G--T-G---------

    >>> list(_read_calc_variations(1, [(0,42), (1,10), (0, 3), (2, 1), (0, 9), (2, 1), (0, 11)], '27G16A0^T6C2^T1C9', 'ACTGAGAAACCCAACCCTCTGAGACCAACACACCCCTTTCAACACATTTTTGGCCGTTCCTGCCCGCCTTCTTTG',  ))
    [(0, 28, 'A'), (1, 43, 'CACATTTTTG'), (0, 45, 'C'), (2, 46, 'T'), (0, 53, 'G'), (2, 56, 'T'), (0, 58, 'G')]


    Real example 2
    577_1692_891_A1
    CIGAR: 34M100N39M2I
    MD:Z:  3T69
                                                          1         1         1         1
                   1         2         3                  4         5         6         7
    pos:  1234567890123456789012345678901234| [100] |567890123456789012345678901234567890123
    ref:  GGATTCTTCCCACTGGGTCGATGTTGTTTGTGAT|-------|CTGAGAGAGAGTTGCATCTGCACATGCTTTCCTGGCGTC^^

    read: GGAATCTTCCCACTGGGTCGATGTTGTTTGTGAT|-------|CTGAGAGAGAGTTGCATCTGCACATGCTTTCCTGGCGTCTC
          MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM  NNNNN  MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMII
          ---A------------------------------         ---------------------------------------TC

    >>> list(_read_calc_variations(1, [(0,34), (3,100), (0, 39), (1, 2)], '3T69', 'GGAATCTTCCCACTGGGTCGATGTTGTTTGTGATCTGAGAGAGAGTTGCATCTGCACATGCTTTCCTGGCGTCTC',  ))
    [(0, 4, 'A'), (1, 174, 'TC')]

    '''

    ref_pos = start_pos
    read_pos = 0

    for op, length in cigar:
        if md and md[0] == '0':
            md = md[1:]
        # sys.stderr.write('%s, %s, %s\n' %(op, length, md))
        if op == 0:  # M
            # how far in the chunk are we? (do *not* update ref_pos until end)
            md_pos = 0
            last = None
            while md and md_pos < length:
                if last == (op, length, md):
                    sys.stderr.write('\nInfinite loop in variant finding!\nPos: %s\nCIGAR: (%s, %s)\n' % (ref_pos, op, length))
                    sys.exit(1)
                last = (op, length, md)
                # sys.stderr.write('%s, %s, %s\n' %(op, length, md))
                chunk_size, md = _extract_md_matches(md, length - md_pos)
                # sys.stderr.write('   -> %s, %s\n' %(chunk_size, md))
                md_pos += chunk_size

                # look for mismatches
                while md_pos < length and md and md[0] not in '0123456789^':
                    yield (op, ref_pos + md_pos, seq[read_pos + md_pos])
                    md = md[1:]

                    md_pos += 1

            ref_pos += length
            read_pos += length

        elif op == 1:  # I
            # nothing in MD about inserts...
            yield (op, ref_pos, seq[read_pos:read_pos + length])
            read_pos += length

        elif op == 2:  # D
            # prefixed with '^' and includes all of the removed bases
            if md[0] == '^':
                md = md[1:]
            yield (op, ref_pos, md[:length])
            md = md[length:]
            ref_pos += length

        elif op == 3:  # N
            ref_pos += length


def read_calc_mismatches_gen(ref, read, chrom):
    start = read.pos
    ref_pos = 0
    read_pos = 0

    for op, length in read.cigar:
        if op == 1:
            yield ref_pos, op, None
            read_pos += length
        elif op == 2:
            yield ref_pos, op, None
            ref_pos += length
        elif op == 3:
            ref_pos += length
        elif op == 0:
            refseq = ref.fetch(chrom, start + ref_pos, start + ref_pos + length)
            cur_pos = start + ref_pos
            for refbase, readbase in zip(refseq.upper(), read.seq[read_pos:read_pos + length].upper()):
                if refbase != readbase:
                    yield op, cur_pos, readbase
                cur_pos += 1
            ref_pos += length
            read_pos += length


def read_calc_mismatches_ref(ref, read, chrom):
    edits = 0

    for op, pos, base in read_calc_mismatches_gen(ref, read, chrom):
        edits += 1

    return edits

if __name__ == '__main__':
    import doctest
    doctest.testmod()
