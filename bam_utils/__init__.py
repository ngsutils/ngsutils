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


def read_calc_mismatches_ref(ref, read, chrom):
    start = read.pos
    edits = 0
    ref_pos = 0
    read_pos = 0

    for op, length in read.cigar:
        if op == 1:
            edits += 1
            read_pos += length
        elif op == 2:
            edits += 1
            ref_pos += length
        elif op == 3:
            ref_pos += length
        elif op == 0:
            refseq = ref.fetch(chrom, start + ref_pos, start + ref_pos + length)
            for refbase, readbase in zip(refseq, read.seq[read_pos:read_pos + length]):
                if refbase.upper() != readbase.upper():
                    edits += 1
            ref_pos += length
            read_pos += length

    return edits
