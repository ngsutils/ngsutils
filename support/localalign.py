#!/usr/bin/env python
'''
Simple Smith-Waterman aligner
'''

import sys
import StringIO


class ScoringMatrix(object):
    def __init__(self, filename=None, text=None):
        assert filename or text

        if filename:
            fs = open(filename)
        else:
            fs = StringIO.StringIO(text)

        self.scores = []
        self.bases = None

        for line in fs:
            if not self.bases:
                self.bases = line.split()
                self.base_count = len(self.bases)
            else:
                cols = line.split()
                self.scores.extend([int(x) for x in cols[1:]])
        fs.close()

    def score(self, one, two):
        one_idx = 0
        two_idx = 0
        for i, b in enumerate(self.bases):
            if b == one:
                one_idx = i
            if b == two:
                two_idx = i

        return self.scores[(one_idx * self.base_count) + two_idx]


class NucleotideScoringMatrix(object):
    def __init__(self, match=1, mismatch=-1):
        self.match = match
        self.mismatch = mismatch

    def score(self, one, two):
        if one == two:
            return self.match
        return self.mismatch


class Matrix(object):
    def __init__(self, rows, cols, init=None):
        self.rows = rows
        self.cols = cols
        self.values = [init, ] * rows * cols

    def get(self, row, col):
        return self.values[(row * self.cols) + col]

    def set(self, row, col, val):
        self.values[(row * self.cols) + col] = val


class LocalAlignment(object):
    def __init__(self, scoring_matrix, gap_penalty=-1):
        self.scoring_matrix = scoring_matrix
        self.gap_penalty = gap_penalty

    def align(self, ref, query):
        ref = ref.upper()
        query = query.upper()

        matrix = Matrix(len(query) + 1, len(ref) + 1, (0, ' '))

        max_val = 0
        max_row = 0
        max_col = 0

        # calculate matrix
        for row in xrange(1, matrix.rows):
            for col in xrange(1, matrix.cols):
                mm = matrix.get(row - 1, col - 1)[0] + self.scoring_matrix.score(query[row - 1], ref[col - 1])
                delete = matrix.get(row - 1, col)[0] + self.gap_penalty
                insert = matrix.get(row, col - 1)[0] + self.gap_penalty
                m = max(0, mm, delete, insert)

                if m == mm:
                    val = (m, 'm')
                elif m == delete:
                    val = (m, 'd')
                elif m == insert:
                    val = (m, 'i')
                else:
                    val = (0, 'x')

                if val[0] >= max_val:
                    max_val = val[0]
                    max_row = row
                    max_col = col

                matrix.set(row, col, val)

        # backtrack
        row = max_row
        col = max_col
        val = max_val
        op = 'm'

        aln = []

        while val > 0:
            aln.append(op)

            op = matrix.get(row, col)[1]

            if op == 'm':
                row -= 1
                col -= 1
            elif op == 'd':
                row -= 1
            elif op == 'i':
                col -= 1
            else:
                break

            val, op = matrix.get(row, col)
        aln.reverse()

        return Alignment(query, ref, row, col, self._reduce_cigar(aln), max_val)

    def _reduce_cigar(self, cigar):
        count = 1
        last = None
        ret = []
        for op in cigar:
            if last and op == last:
                count += 1
            elif last:
                ret.append((count, last.upper()))
                count = 1
            last = op

        if last:
            ret.append((count,last.upper()))
        return ret

    def dump_matrix(self, ref, query, matrix, show_row=-1, show_col=-1):
        sys.stdout.write('   -   ')
        sys.stdout.write('   '.join(ref))
        sys.stdout.write('\n')
        for row in xrange(matrix.rows):
            if row == 0:
                sys.stdout.write('-')
            else:
                sys.stdout.write(query[row - 1])

            for col in xrange(matrix.cols):
                if show_row == row and show_col == col:
                    sys.stdout.write('   *')
                else:
                    sys.stdout.write(' %2s%s' % (matrix.get(row, col)[0], matrix.get(row, col)[1]))
            sys.stdout.write('\n')


class Alignment(object):
    def __init__(self, query, ref, q_pos, r_pos, cigar, score):
        self.query = query
        self.ref = ref
        self.q_pos = q_pos
        self.r_pos = r_pos
        self.cigar = cigar
        self.score = score

        q_len = 0
        r_len = 0

        self.matches = 0
        self.mismatches = 0

        i = self.q_pos
        j = self.r_pos

        for count, op in self.cigar:
            if op == 'M':
                q_len += count
                r_len += count
                for k in xrange(count):
                    if self.query[i] == self.ref[j]:
                        self.matches += 1
                    else:
                        self.mismatches += 1
                    i += 1
                    j += 1

            elif op == 'I':
                r_len += count
                for k in xrange(count):
                    j += 1
                    self.mismatches += 1
            elif op == 'D':
                q_len += count
                for k in xrange(count):
                    i += 1
                    self.mismatches += 1

        self.q_end = q_pos + q_len
        self.r_end = r_pos + r_len
        if self.mismatches + self.matches > 0:
            self.identity = float(self.matches) / (self.mismatches + self.matches)
        else:
            self.identity = 0

    def dump(self, out=sys.stdout):
        i = self.q_pos
        j = self.r_pos

        q = '%4s ' % self.q_pos
        r = '%4s ' % self.r_pos
        m = '     '

        for count, op in self.cigar:
            if op == 'M':
                for k in xrange(count):
                    q += self.query[i]
                    r += self.ref[j]
                    if self.query[i] == self.ref[j]:
                        m += '|'
                    else:
                        m += ' '

                    i += 1
                    j += 1
            elif op == 'I':
                for k in xrange(count):
                    q += '-'
                    r += self.ref[j]
                    m += ' '
                    j += 1
            elif op == 'D':
                for k in xrange(count):
                    q += self.query[i]
                    r += '-'
                    m += ' '
                    i += 1

            elif op == 'N':
                q += '-//-'
                r += '-//-'
                m += '    '

        out.write("%s %s\n" % (q, self.q_end))
        out.write("%s\n" % m)
        out.write("%s %s\n" % (r, self.r_end))
        out.write("Score: %s\n" % self.score)
        out.write("Matches: %s (%.1f%%)\n" % (self.matches, self.identity * 100))
        out.write("Mismatches: %s\n" % (self.mismatches,))

if __name__ == '__main__':
    sw = LocalAlignment(NucleotideScoringMatrix(2, -1))
    sw.align(sys.argv[1], sys.argv[2]).dump()
#     sw.align('ACACACTA','AGCACACA').dump()
#     aln=sw.align("AAGGGGAGGACGATGCGGATGTTC","AGGGAGGACGATGCGG")
#     aln.dump()
