#!/usr/bin/env python
## category General
## desc Splits a FASTQ/FASTA file based on sequence barcodes
'''
Given a FASTA/FASTQ file and a list of barcodes (name and seq),
the FASTA/FASTQ file is split into N number of new files representing
the reads for each barcode. The format of the input file is determined
automatically. The barcode will be removed from read sequence (and quality
score if FASTQ).

barcode input format:
tagname\tsequence\torientation (5 or 3)\tstrip tag? (y/n, optional - default y)

The output files will be named: out_template_tagname.fast[qa]

Reads missing a barcode (or if the tag is too degenerate) will be written to:
out_template_missing.fast[qa]

Note: This is a slow method using a Smith-Waterman alignment algorithm. It is
      primarily useful with data that is inherently noisy, such as PacBio
      sequencing reads. For less-noisy sequences, a faster, non-alignment
      approach will be better.

Note 2: This isn't appropriate for color-space FASTQ files with a prefix base
        included in the read sequence, since it trims an equal number of bases
        from the sequence and quality FASTQ lines.
'''

import sys
import os
import gzip

from ngsutils.support import revcomp, FASTA
from ngsutils.fastq import FASTQ

import swalign

sw = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(2, -1))


def fastx_barcode_split(reader, outtempl, barcodes, edits=0, pos=0, allow_revcomp=False, gzip_output=False, stats_fname=None):
    '''
    Split FAST[QA] reads from {fname} using {barcodes} (hash) to write them to
    output files named like {templ}.
    '''

    outs = {}
    if gzip_output:
        outtempl += '.gz'
        outs[''] = gzip.open(outtempl % 'missing', 'w')
    else:
        outs[''] = open(outtempl % 'missing', 'w')

    tag_count = {}
    for tag in barcodes:
        if gzip_output:
            outs[tag] = gzip.open(outtempl % tag, 'w')
        else:
            outs[tag] = open(outtempl % tag, 'w')
        tag_count[tag] = 0

    matched = 0
    perfect = 0
    mismatched = 0
    mispositioned = 0
    missing = 0

    for record in reader.fetch():
        ismatch, (tag, aln, is_forward, reason) = check_tags(barcodes, record.seq, edits, pos, allow_revcomp)

        if not ismatch:
            missing += 1
            comment = '#guess: %s, %s%s (%s, %s, %smm) %s' % (tag, barcodes[tag][0], '' if is_forward else '[rc]', aln.q_pos, aln.extended_cigar_str, aln.mismatches, reason)
            record.subseq(None, None, comment).write(outs[''])
        else:
            matched += 1
            if aln.mismatches == 0:
                if aln.q_pos == 0 or aln.q_end == len(record.seq):
                    perfect += 1
                else:
                    mispositioned += 1
            else:
                mismatched += 1

            comment = '#%s%s (%s, %s)' % (tag, '' if is_forward else '[rc]', aln.q_pos, aln.mismatches,)

            if not barcodes[tag][2]:
                record.write(outs[tag])

            else:
                if is_forward:
                    if barcodes[tag][1] == '5':
                        newrec = record.subseq(aln.q_end, None, comment)
                    else:
                        subseqlen = len(barcodes[tag][0]) + edits + pos
                        newrec = record.subseq(None, -1 * (subseqlen - aln.q_pos), comment)
                else:
                    if barcodes[tag][1] == '5':
                        subseqlen = len(barcodes[tag][0]) + edits + pos
                        newrec = record.subseq(None, -1 * (subseqlen - aln.q_pos), comment)
                    else:
                        newrec = record.subseq(aln.q_end, None, comment)

                newrec.write(outs[tag])

            tag_count[tag] += 1

    for tag in outs:
        outs[tag].close()

    if stats_fname:
        with open(stats_fname, 'w') as f:
            f.write("Unmatched\t%s\n" % missing)
            f.write("Matched\t%s\n" % matched)
            f.write(" Perfect\t%s\n" % perfect)
            f.write(" Mismatched\t%s\n" % mismatched)
            f.write(" Mispositioned\t%s\n" % mispositioned)
            f.write('\n')
            for tag in tag_count:
                f.write("Tag\t%s\t%s\n" % (tag, tag_count[tag]))


def _tag_aln_check(aln, seqlen, barcodelen, orientation, edit=0, pos=0):
    r_len = aln.r_end - aln.r_pos
    mismatches = barcodelen - r_len + aln.mismatches

    if mismatches <= edit:
        if orientation == '5' and aln.q_pos <= pos:  # allow at most {pos} extra bases at start
            return True, ''
        elif orientation == '3' and aln.q_end >= seqlen - pos:  # allow at most {pos} extra bases at end
            return True, ''
        else:
            return False, 'Not at right location'
    else:
        return False, 'too many mismatches'

    return False, '?'


def check_tags(barcodes, seq, edit, pos, allow_revcomp=False, verbose=False):
    '''
    For each barcode, pull out the appropriate 5' or 3' sub sequence from {seq}. Then
    run a local alignment of the barcode to the subseq. If a good match is found, return
    it; otherwise, find the best match and return that.

    returns a multi-tuple:
        valid, (tag, valid_seq, edits, reason_if_fail)

    For the alignments, the reference is the barcode, the query is the subset of the read
    that is possibly the barcode (5'/3' subseq)
    '''

    best = None

    # check perfect matches first...
    # for tag in barcodes:
    #     barcodeseq, orientation = barcodes[tag]
    #     if orientation == '5':
    #         if seq[:len(barcodeseq)] == barcodeseq:
    #             return True, (tag, seq[len(barcodeseq):], 0, '')
    #     else:
    #         if seq[-len(barcodeseq):] == barcodeseq:
    #             return True, (tag, seq[:-len(barcodeseq)], 0, '')

    for tag in barcodes:
        barcodeseq, orientation, strip = barcodes[tag]
        if orientation == '5':
            testseq = seq[:len(barcodeseq) + edit + pos]
        else:
            testseq = seq[-1 * (len(barcodeseq) + edit + pos):]

        aln = sw.align(barcodeseq, testseq)
        valid, reason = _tag_aln_check(aln, len(testseq), len(barcodeseq), orientation, edit, pos)
        if verbose:
            print 'Testing tag: %s vs %s' % (str(barcodes[tag]), testseq)
            aln.dump()
            print valid, reason
        if valid:
            return True, (tag, aln, True, '')

        if not best or aln.score > best[1].score:
            best = (tag, aln, True, reason)

        if allow_revcomp:
            if orientation == '5':
                testseq = seq[-1 * (len(barcodeseq) + edit + pos):]
            else:
                testseq = seq[:len(barcodeseq) + edit + pos]

            aln = sw.align(revcomp(barcodeseq), testseq)
            valid, reason = _tag_aln_check(aln, len(testseq), len(barcodeseq), '5' if orientation == '3' else '3', edit, pos)
            if verbose:
                print 'Testing tag: %s [rc] vs %s' % (str(barcodes[tag]), testseq)
                aln.dump()
                print valid, reason
            if valid:
                return True, (tag, aln, False, '')

            if not best or aln.score > best[1].score:
                best = (tag, aln, False, reason)

    if verbose:
        print 'BEST: ', best
        best[1].dump()

    return False, best


def usage():
    sys.stdout.write('%s\n' % __doc__)
    sys.stdout.write('''
Usage: %s {options} barcodes.txt input.fastx{.gz} output_template

Options:
  -edit num         Number of mismatches/indels to allow in barcode discovery.
                    (default: 0)

  -pos num          Allow barcode to be within [num] bases from the start.
                    (or end for 3' barcodes) (default: 0)

  -allow-revcomp    Allow barcodes to match the reverse-compliment of a read.
                    (non-strand specific sequencing) The read's orientation
                    will *not* be changed in the output file.

  -gz               GZip compress the output files

  -stats            Output stats file (output_template.stats.txt)

''' % os.path.basename(sys.argv[0]))
    sys.exit(1)


def is_fasta_file(fname):
    '''
    Open a fileobj and read the first line
        If it starts with '>', it's FASTA
        If it starts with '@', it's FASTQ
    Return the actual FASTX object
    '''
    if fname[-3:] == '.gz':
        f = gzip.open(fname)
    else:
        f = open(fname)

    try:
        while True:
            line = f.next()
            if line[0] == '>':
                f.close()
                return True
            elif line[0] == '@':
                f.close()
                return False
    except:
        f.close()
        raise ValueError("Unknown type of file!")

if __name__ == '__main__':
    barcodes = None
    fname = None
    templ = None
    edit = 0
    pos = 0
    allow_revcomp = False
    gz = False
    stats = False

    last = None

    for arg in sys.argv[1:]:
        if last == '-edit':
            edit = int(arg)
            last = None
        elif last == '-pos':
            pos = int(arg)
            last = None
        elif arg in ['-edit', '-pos']:
            last = arg
        elif arg == '-allow-revcomp':
            allow_revcomp = True
        elif arg == '-gz':
            gz = True
        elif arg == '-stats':
            stats = True
        elif not barcodes and os.path.exists(arg):
            barcodes = {}
            with open(arg) as f:
                first = True
                for line in f:
                    if first:
                        first = False
                        continue
                    cols = line.strip().split('\t')
                    if len(cols) < 4:
                        cols.append('y')
                    barcodes[cols[0]] = (cols[1], cols[2], True if cols[3] in 'TtyY' else False)
        elif not fname and os.path.exists(arg):
            fname = arg
        elif not templ:
            templ = arg

    if not barcodes or not fname or not templ:
        usage()

    if is_fasta_file(fname):
        outtempl = '%s_%%s.fasta' % templ
        reader = FASTA(fname)
    else:
        outtempl = '%s_%%s.fastq' % templ
        reader = FASTQ(fname)

    if stats:
        stats_fname = '%s.stats.txt' % templ
    else:
        stats_fname = None

    fastx_barcode_split(reader, outtempl, barcodes, edit, pos, allow_revcomp, gz, stats_fname)
