#!/usr/bin/python
## category General
## desc Converts to/from BFQ (binary-compressed FASTQ format)
## experimental
'''
Converts to/from BFQ (binary-compressed FASTQ format)

This stores the name/seq/qual information for a FASTQ file in a highly compact
binary format. This stores the sequence + qual in an 8-bit encoded format. It
can store more than one sequence+qual fragments for each read, meaning that it
efficiently deals with paired-end data (with support for upto 256 fragments).

To save paired end data, include one file for each fragment, or use a single
FASTQ file with one record for each fragment (consecuatively using the same
name for each fragment).

BFQ also compresses the read data with zlib, so it is very compact. (This can
be turned off for faster processing).

Finally, the file includes an MD5 checksum at the end so that file integrity
be confirmed at any point.
'''

'''
8 bit-encoding:
bits: SSQQ QQQQ

Sequence is a 2-bit format (0 -> A, 1-> C, 2-> G, 3-> T for base-space)
Seq can be obtained by: byte >> 6

Quality is a value from 0-62. This can be obtained by: byte & 0x3F. Anything
higher is truncated to 62.

If the entire byte is 255 (0xFF), then the base is a wildcard (N) and the
quality value is zero.

All values are little-endian, regardless of platform.

Limits:
    length of description: 4 billion (uint)
    length of a fragment name: 64K (ushort)
    number of fragments per read: 256 (uchar)
    length of read name: 64K (ushort)
    length of sequence: 64K (ushort)
    quality score: 0-62

File format:

[bfq]
{header}
{body}
16byte  md5 hash of entire file

[header]
uint    magic_byte    (always uncompressed)
ushort  file version  (always uncompressed)
ushort  flags bitmap  (always uncompressed)
(compression starts here - if enabled)
uchar   number of fragments per read
uint    description length
char*   description
{fragment_header}+

    Valid flags:
        0x1     The file is zlib compressed
        0x2     The description contains a filename followed by binary data
                The format is:
                  filename\0binary data


[fragment_header]
ushort  fragment flags bitmap
ushort  name length
char*   name

    Valid fragment flags:
        0x1     Fragment is in colorspace
        0x2     If fragment is in colorspace, the last base of the adapter
                is included (in base-space)

[body]
{read}+

[read]
ushort  length of the name (must be > 0)
char*   read name
{seq_qual}+ (for each fragment, if a read is missing a fragment, the length is 0)

[seq_qual]
ushort  length of the sequence
uchar*  seq+qual in 8bit encoding

'''

import sys
import os
import gzip
import struct
import zlib
import hashlib
import collections

_BFQ_fragment_flags = collections.namedtuple('_BFQ_fragment_flags', 'flags colorspace cs_include_prefix')


class BFQ(object):
    _magic = 0xE1EEBEA4
    __cs_encode = {'0': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 4, '6': 4, '.': 4}
    __cs_decode = '0123'
    __nt_encode = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    __nt_decode = 'ACGT'

    def __init__(self, filename, mode='r', description='', fragment_count=1, fragment_names=None, fragment_flags=None, description_binary=False):
        assert mode in ['r', 'w', 'wz']
        self.filename = filename
        self.description = description
        self.fragment_count = fragment_count
        self.fragment_flags = fragment_flags
        self.fragment_names = fragment_names

        self.description_binary = description_binary
        self.compress = False
        self.version = 1
        self.flags = 0

        self.__fragment_colorspace = []
        self.__fragment_cs_include_prefix = []

        self.__errors_printed = set()
        self.__last_name = None
        self.__seqnum = 0

        mode = mode.lower()

        if mode == 'wz':
            self.mode = 'w'
            self.compress = True
            self.__compressor = zlib.compressobj()
        else:
            self.mode = mode

        if filename == '-':
            if self.mode == 'r':
                self.fobj = sys.stdin
            else:
                self.fobj = sys.stdout
        else:
            self.fobj = open(filename, '%sb' % self.mode)

        if self.mode == 'w':
            self.__wrote_header = False
            self.__md5 = hashlib.md5()

            if not self.fragment_flags:
                self.fragment_flags = [0, ] * fragment_count
            if not self.fragment_names:
                self.fragment_names = ['', ] * fragment_count

        else:
            self._read_header()

    def __iter__(self):
        return self

    def next(self):
        try:
            lenname, = self.__read_struct('<H')
            name, = self.__read_struct('<%ss' % lenname)

            seqs = []
            quals = []

            for i in xrange(self.fragment_count):
                lenseq = self.__read_struct('<H')

                if lenseq == 0:
                    seqs.append(None)
                    quals.append(None)
                    continue

                seqqual = self.__read_struct('<%sB' % lenseq)

                seq = []
                qual = []
                if self.__fragment_cs_include_prefix[i]:
                        seq.append(BFQ.__nt_decode[seqqual[0] >> 6])
                        seqqual = seqqual[1:]

                for sq in seqqual:
                    if sq == 0xFF:
                        if self.__fragment_colorspace[i]:
                            seq.append('.')
                        else:
                            seq.append('N')
                        qual.append(0)
                    else:
                        if self.__fragment_colorspace[i]:
                            seq.append(BFQ.__cs_decode[sq >> 6])
                        else:
                            seq.append(BFQ.__nt_decode[sq >> 6])

                        qual.append(sq & 0x3F)

                seqs.append(''.join(seq))
                quals.append(qual)
            return (name, seqs, quals)
        except Exception:
            raise StopIteration

    def write_read(self, name, seq, qual):
        if not self.__wrote_header:
            self._write_header()

        if self.fragment_count == 1 or name != self.__last_name:
            if self.__last_name:
                while self.__seqnum < self.fragment_count:
                    self.__write_data(struct.pack('<H', 0))
                    self.__seqnum += 1
            self.__write_data(struct.pack('<H%ss' % len(name), len(name), name))
            self.__seqnum = 0

        self.__last_name = name
        self.__write_data(struct.pack('<H', len(seq)))

        if self.__fragment_cs_include_prefix[self.__seqnum]:
                self.__write_data(struct.pack('<B', (BFQ.__nt_encode[seq[0]] << 6)))
                seq = seq[1:]

        for s, q in zip(seq, qual):
            self.__write_base_qual(s, q)
        self.__seqnum += 1

    def tell(self):
        return self.fobj.tell()

    def close(self):
        if self.mode == 'w':
            if self.compress:
                data = self.__compressor.flush()
                if data:
                    self.fobj.write(data)
                    self.__md5_update(data)
            self.fobj.write(self.__md5.digest())
            self.fobj.flush()
        self.fobj.close()

    def _read_header(self):
        magic, = struct.unpack('<I', self.fobj.read(4))

        if magic != BFQ._magic:
            raise IOError("File isn't in BFQ format")

        self.version, self.flags, self.fragment_count, lendesc = struct.unpack('<HHBI', self.fobj.read(9))
        self.compress = self.flags & 0x1 == 0x1
        self.description_binary = self.flags & 0x2 == 0x2

        if lendesc > 0:
            self.description, = struct.unpack('<%ss' % lendesc, self.fobj.read(lendesc))

        self.fragment_names = []
        self.fragment_flags = []

        for i in xrange(self.fragment_count):
            flags, namelen = struct.unpack('<HH', self.fobj.read(4))

            if namelen > 0:
                name, = struct.unpack('<%ss' % namelen, self.fobj.read(namelen))
            else:
                name = ''

            self.fragment_names.append(name)
            self.fragment_flags.append(flags)

            self.__fragment_colorspace.append(flags & 0x1 == 0x1)
            self.__fragment_cs_include_prefix.append(flags & 0x2 == 0x2)

        if self.compress:
            self.__decompressor = zlib.decompressobj()
            self.__read_buffer = ''

    def _write_header(self):
        flags = 0

        if self.compress:
            flags = flags | 0x1
        if self.description_binary:
            flags = flags | 0x2

        data = struct.pack('<IHH', BFQ._magic, self.version, flags)
        self.__write_data(data, True)

        data = struct.pack('<BI%ss' % len(self.description), self.fragment_count, len(self.description), self.description)
        self.__write_data(data)

        while len(self.fragment_names) < len(self.fragment_flags):
            self.fragment_names.append('')

        for name, ff in zip(self.fragment_names, self.fragment_flags):
            data = struct.pack('<HH%ss' % len(name), ff, len(name), name)
            self.__write_data(data)

            self.__fragment_colorspace.append(ff & 0x1 == 0x1)
            self.__fragment_cs_include_prefix.append(ff & 0x2 == 0x2)

        self.__wrote_header = True

    def __write_base_qual(self, base, qual):
        if self.__fragment_colorspace[self.__seqnum]:
            b = BFQ.__cs_encode[base]
        else:
            b = BFQ.__nt_encode[base]

        if b == 4:
            b = 0xFF
        else:
            if qual >= 0x3F:
                self.__errmsg("Warning: Quality value truncated")
                qual = 0x3E
            b = (b << 6) | (qual & 0x3F)

        self.__write_data(struct.pack('<B', b))

    def __errmsg(self, msg):
        if not msg in self.__errors_printed:
            self.__errors_printed.add(msg)
            sys.stderr.write(msg)

    def __read_struct(self, fmt):
        size = struct.calcsize(fmt)

        if self.compress:
            while len(self.__read_buffer) < size:
                b = self.fobj.read(16384)
                if not b:
                    return
                self.__read_buffer += self.__decompressor.decompress(b)
            data = self.__read_buffer[:size]
            self.__read_buffer = self.__read_buffer[size:]
        else:
            data = self.fobj.read(size)

        if len(data) != size:
            self.__errmsg('Error: File truncated')

        return struct.unpack(fmt, data)

    def __write_data(self, data, force_uncompressed=False):
        if self.compress and not force_uncompressed:
            chunk = self.__compressor.compress(data)
            if chunk:
                self.fobj.write(chunk)
                self.__md5_update(chunk)
        else:
            self.fobj.write(data)
            self.__md5_update(data)

    def __md5_update(self, data):
        self.__md5.update(data)


def _check_fragment_flags(seq):
    flags = 0

    if '0' in seq or '1' in seq or '2' in seq or '3' in seq or '4' in seq or '.' in seq:
        flags = flags | 0x1

        if seq[0].upper() in 'ATCG':
            flags = flags | 0x2

    return flags


def _read_fastq(fobjs, fileno=-1):
    try:
        i = fileno + 1
        while True:
            if i >= len(fobjs):
                i = 0
            name = fobjs[i].next()
            if not name:
                return
            seq = fobjs[i].next()
            fobjs[i].next()
            qual = fobjs[i].next()
            # sys.stderr.write('*%s*' % name.strip())
            yield (name, seq, qual, i)

            if len(fobjs) > 1:
                i += 1

    except Exception:
        return


def bfq_decode(fname):
    bfq = BFQ(fname)

    for name, seqs, quals in bfq:
        for seq, qual in zip(seqs, quals):
            sys.stdout.write('@%s\n%s\n+\n%s\n' % (name, seq, ''.join([chr(x + 33) for x in qual])))

    bfq.close()


def bfq_encode_fastq(fnames, stdout=False, compress=False, description='', outname=None, force=False, quiet=False, fragment_names=None, description_binary=False):
    fobjs = []
    for fname in fnames:
        if fname == '-':
            if not os.isatty(0):
                fsize = 0
                f = sys.stdin
                quiet = True
            else:
                sys.stderr.write('Missing data on stdin\n')
                usage()
        elif fname[-3:].lower() == '.gz':
            f = gzip.open(fname)
            fsize = os.stat(fname).st_size
            teller = f.fileobj
        else:
            f = open(fname)
            fsize = os.stat(fname).st_size
            teller = f
        fobjs.append(f)

    buf = []
    flags = {}
    last_fileno = 0
    for name, seq, qual, fileno in _read_fastq(fobjs):
        buf.append((name, seq, qual))

        if not name in flags:
            flags[name] = []

        flags[name].append(_check_fragment_flags(seq))

        if len(flags) > 10:
            last_fileno = fileno
            break

    nc = []
    for name in flags:
        nc.append((len(flags[name]), flags[name]))

    nc.sort()
    fragment_count = nc[-1][0]
    fragment_flags = nc[-1][1]

    if compress:
        mode = 'wz'
    else:
        mode = 'w'

    if stdout:
        outname = '-'
    elif not outname:
        if fnames[0] == '-':
            outname = '-'
        else:
            spl = fnames[0].split('.')
            base = []
            skipped = False
            for ext in spl[::-1]:
                if ext.lower() not in ['gz', 'fastq', 'fq', 'txt'] or skipped:
                    base.insert(0, ext)
                    skipped = True

            base.append('bfq')
            outname = '.'.join(base)

    if outname and os.path.exists(outname) and not force:
        sys.stderr.write('%s exists!\n' % outname)
        sys.exit(1)

    if not fragment_names:
        fragment_names = []
        for f in fnames:
            fragment_names.append(os.path.basename(fname))

    bfq = BFQ(outname, mode, description=description, fragment_count=fragment_count, fragment_names=fragment_names, fragment_flags=fragment_flags, description_binary=description_binary)

    last_pos = 0
    diff = .05 * fsize
    pos = 0

    if not quiet:
        sys.stderr.write('[...................]\r[>')

    for name, seq, qual in buf:
        bfq.write_read(name[1:].strip(), seq.strip(), [ord(q) - 33 for q in qual.strip()])

    for name, seq, qual, fileno in _read_fastq(fobjs, last_fileno):
        bfq.write_read(name[1:].strip(), seq.strip(), [ord(q) - 33 for q in qual.strip()])

        if fsize and not quiet:
            pos = teller.tell()
            if pos - last_pos > diff:
                sys.stderr.write('>')
                sys.stderr.flush()
                last_pos = pos

    bfq.close()

    if not quiet:
        sys.stderr.write('\r                     \r')

    for fobj in fobjs:
        if fobj != sys.stdin:
            fobj.close()


def bfq_info(fnames):
    for fname in fnames:
        bfq = BFQ(fname)
        sys.stdout.write('[%s]\n' % fname if fname != '-' else 'stdin')
        if bfq.description:
            sys.stdout.write('[Description]\n')
            if bfq.description_binary:
                for c in bfq.description:
                    if c == '\0':
                        sys.stdout.write(' *binary*\n')
                        break
                    else:
                        sys.stdout.write(c)
            else:
                sys.stdout.write('%s\n' % bfq.description)
        sys.stdout.write('[Flags]\n0x%02x ' % bfq.flags)
        f = []
        if bfq.flags & 0x01 == 0x01:
            f.append('compressed')
        else:
            f.append('not-compressed')

        sys.stdout.write('(%s)\n' % (','.join(f)))
        sys.stdout.write('[Fragments]\n%s fragment(s)\n' % bfq.fragment_count)
        i = 0
        for n, ff in zip(bfq.fragment_names, bfq.fragment_flags):
            i += 1
            sys.stdout.write('  #%-4s: %s\n  Flags: 0x%02x ' % (i, n, ff))
            f = []
            if ff & 0x01 == 0x01:
                f.append('colorspace')
            else:
                f.append('basespace')
            if ff & 0x02 == 0x02:
                f.append('cs_prefix')

            sys.stdout.write('(%s)\n' % (','.join(f)))
        sys.stdout.write('\n')


def bfq_description(fname):
    bfq = BFQ(fname)
    if bfq.description:
        i = 0
        if bfq.description_binary:
            for c in bfq.description:
                i += 1
                if c == '\0':
                    break
        sys.stdout.write(bfq.description[i:])


def bfq_test(fnames):
    if type(fnames) == str:
        fnames = [fnames, ]

    allgood = True
    for fname in fnames:
        bfq = BFQ(fname)
        bfq.close()
        sys.stdout.write('%s\t' % fname if fname != '-' else 'stdin')

        f = open(fname)
        f.seek(-16, 2)
        hash_pos = f.tell()
        known_md5 = f.read(16)
        sys.stdout.write(''.join(['%02x' % ord(b) for b in known_md5]))

        md5 = hashlib.md5()
        f.seek(0)
        buf = ''

        bufsize = 16 * 1024
        pos = 0
        while pos + bufsize < hash_pos:
            buf = f.read(bufsize)
            md5.update(buf)
            pos += len(buf)

        buf = f.read(hash_pos - pos)
        md5.update(buf)
        f.close()

        digest = md5.digest()

        if digest == known_md5:
            sys.stdout.write('\tOK\n')
        else:
            sys.stdout.write('\tERROR (')
            sys.stdout.write(''.join(['%02x' % ord(b) for b in digest]))
            sys.stdout.write(')\n')
            allgood = False

    if not allgood:
        sys.exit(1)


def usage():
    print __doc__
    print """\
Usage: fastqutils bfq {opts} infile1.fastq {infile2.fastq ... }
       fastqutils bfq -d {opts} filename.bfq
       fastqutils bfq -t filename.bfq

Options:
    -d              Decode BFQ file to stdout in FASTQ format (one-file)
    -t              Test file integrity
    -i              Display information about the file and fragments

    -h              Display this message

Encoding options:
    -c              Output to stdout
    -f              Force overwriting existing files
    -fn             Include a name for a fragment (there can be as
                    many of these are there are fragments)
    -desc text      Include a description of the file
    -descf fname    The contents of the file fname will be used as the
                    description (this can be a binary file)
    -nc             Disable compression
    -o fname        Name of the output file
                    (defaults to input filename.bfq)
    -q              Quiet (no progress bar)

Decoding options:
    -desc           Extract the description only (to stdout)

"""
    sys.exit(1)

if __name__ == '__main__':
    action = 'encode'
    stdout = False
    desc = ''
    desc_bin = False
    compress = True
    force = False
    quiet = False
    outname = None

    fragment_names = []
    fnames = []

    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()

        if last == '-o':
            outname = arg
            last = None
        elif last == '-descf':
            if os.path.exists(arg):
                with open(arg) as f:
                    desc = '%s\0%s' % (os.path.basename(arg), f.read())
                    desc_bin = True
            else:
                sys.stderr.write('Description file: %s is missing!\n' % arg)
                sys.exit(1)
            last = None
        elif last == '-desc':
            desc = arg
            last = None
        elif last == '-fn':
            fragment_names.append(arg)
            last = None
        elif arg == '-desc' and action == 'decode':
            action = 'desc'
        elif arg in ['-desc', '-o', '-fn', '-descf']:
            last = arg
        elif arg == '-nc':
            compress = False
        elif arg == '-f':
            force = True
        elif arg == '-c':
            stdout = True
        elif arg == '-q':
            quiet = True
        elif arg == '-i':
            action = 'info'
        elif arg == '-t':
            action = 'test'
        elif arg == '-d':
            action = 'decode'
        elif arg == '-' or os.path.exists(arg):
            fnames.append(arg)

    # try:
    if True:
        if action == 'test':
            if not fnames:
                usage()
            bfq_test(fnames)
        elif action == 'decode':
            if not fnames:
                fnames.append('-')
            else:
                bfq_decode(fnames[0])
        elif action == 'info':
            if not fnames:
                fnames.append('-')
            bfq_info(fnames)
        elif action == 'desc':
            if not fnames:
                fnames.append('-')
            bfq_description(fnames[0])
        elif action == 'encode':
            if not fnames:
                fnames = ['-']
            bfq_encode_fastq(fnames, stdout=stdout, compress=compress, description=desc, outname=outname, force=force, fragment_names=fragment_names, quiet=quiet, description_binary=desc_bin)
    # except Exception, e:
       # sys.stderr.write('%s: %s' % (type(e).__name__,e))
       # sys.stderr.write('\n')
       # sys.exit(1)
