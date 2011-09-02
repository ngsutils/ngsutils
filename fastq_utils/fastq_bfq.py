#!/usr/bin/python
'''
BFQ - A Binary FastQ format. This stores the name/seq/qual information for a 
FASTQ file in a highly compact binary format. This stores the name as a char*,
and the sequence + qual in an 8-bit encoded format. It can store more than one 
sequence+qual fragments for each read, meaning that it efficiently deals with 
paired-end data (with support for upto 256 fragments).

BFQ also compresses the read data with zlib, so it is very compact. (This can 
be turned off for faster processing).

Finally, the file includes an MD5 checksum at the end so that file integrity 
be confirmed at any point.

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
    length of name: 64K (ushort)
    number of fragments per read: 256 (uchar)
    length of sequence: 64K (ushort)
    quality score: 0-62

File format:

[bfq]
{header}
{body}
16byte  md5 hash of entire file

[header]
uint    magic_byte
ushort  file version
ushort  flags bitmap
uchar   number of sequences per read (paired-end/multiple fragments)
uint    description length
char*   description

        Valid flags:
        0x1     File is in colorspace
        0x2     If file is in colorspace, the last base of the adapter is 
                included (in base-space)
        0x4     The body is zlib compressed

[body]
{read}+

[read]
ushort  length of the name (must be > 0)
char*   read name
{seq_qual}+

[seq_qual]
ushort  length of the sequence
uchar*  seq+qual

'''

import sys,os,gzip,struct,zlib,hashlib

class BFQ(object):
    _magic = 0xE1EEBEA4
    __cs_encode = { '0': 0,'1':1,'2':2,'3':3,'4':4,'5':4,'6':4,'.':4}
    __nt_encode = { 'A': 0,'C':1,'G':2,'T':3,'N':4}
    __nt_decode = 'ACGT'
    
    def __init__(self, filename, mode='r', description='', seq_count = 1):
        assert mode in ['r','w','wz']
        self.filename = filename
        self.description = description
        self.sequences_per_read = seq_count

        self.colorspace = False
        self.cs_include_prefix = False
        self.compress = False
        self.version = 1
        self.__errors_printed = set()
        self.__last_name = None
        self.__seqnum = 1

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
            self.fobj = open(filename,'%sb' % self.mode)
        
        if self.mode == 'w':
            self.__wrote_header = False
            self.__md5 = hashlib.md5()
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
            
            for i in xrange(self.sequences_per_read):
                lenseq = self.__read_struct('<H')
                
                if lenseq == 0:
                    seqs.append(None)
                    quals.append(None)
                    continue

                seqqual = self.__read_struct('<%sB' % lenseq)
                
                seq = []
                qual = []
                if self.colorspace:
                    if self.cs_include_prefix:
                        seq.append(BFQ.__nt_decode[seqqual[0] >> 6])
                        seqqual = seqqual[1:]
                
                for sq in seqqual:
                    if sq == 0xFF:
                        if self.colorspace:
                            seq.append('.')
                        else:
                            seq.append('N')
                        qual.append(0)
                    else:
                        if self.colorspace:
                            seq.append(str(sq >> 6))
                        else:
                            seq.append(BFQ.__nt_decode[sq >> 6])
            
                        qual.append(sq & 0x3F)
                
                seqs.append(''.join(seq))
                quals.append(qual)
            return (name,seqs,quals)
        except Exception, e:
            raise StopIteration
        
    
    def write_read(self,name,seq,qual):
        if not self.__wrote_header:
            self.__check_format(name,seq,qual)
            self._write_header()

        if self.sequences_per_read == 1 or name != self.__last_name:
            while self.__seqnum < self.sequences_per_read:
                self.__write_struct(struct.pack('<H',0))
                self.__seqnum += 1
            self.__write_struct(struct.pack('<H%ss' % len(name),len(name),name))
            self.__seqnum = 1
        else:
            self.__seqnum += 1
        
        self.__write_struct(struct.pack('<H',len(seq)))
            
        if self.colorspace:
            if self.cs_include_prefix:
                self.__write_struct(struct.pack('<B',(BFQ.__nt_encode[seq[0]] << 6)))
                seq = seq[1:]

        for s,q in zip(seq,qual):
            self.__write_base_qual(s,q)

    def tell(self):
        return self.fobj.tell()
    
    def close(self):
        if self.mode == 'w':
            data = self.__compressor.flush()
            if data:
                self.fobj.write(data)
                self.__md5_update(data)
            self.fobj.write(self.__md5.digest())
            self.fobj.flush()
        self.fobj.close()
        

    def _read_header(self):
        header = self.fobj.read(13)
        magic, self.version, flags, self.sequences_per_read, lendesc = struct.unpack('<IHHBI',header)

        if magic != BFQ._magic:
            raise IOError,"File isn't in BFQ format"

        if lendesc > 0:
            desc = self.fobj.read(lendesc)
            self.description, = struct.unpack('<%ss' % lendesc,desc)


        self.colorspace = flags & 0x1 == 0x1
        self.cs_include_prefix = flags & 0x2 == 0x2
        self.compress = flags & 0x4 == 0x4
        
        if self.compress:
            self.__decompressor = zlib.decompressobj()
            self.__read_buffer = ''

    def _write_header(self):
        flags = 0
        
        if self.colorspace:
            flags = flags | 0x1
        if self.cs_include_prefix:
            flags = flags | 0x2
        if self.compress:
            flags = flags | 0x4
        
        data = struct.pack('<IHHBI%ss' % len(self.description),BFQ._magic,self.version,flags,self.sequences_per_read,len(self.description),self.description)
        self.fobj.write(data)
        self.__md5_update(data)
        self.__wrote_header = True
    
    def __write_base_qual(self,base,qual):
        if self.colorspace:
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
        
        self.__write_struct(struct.pack('<B',b))
    
    def __errmsg(self,msg):
        if not msg in self.__errors_printed:
            self.__errors_printed.add(msg)
            sys.stderr.write(msg)

    def __read_struct(self,fmt):
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
        
        return struct.unpack(fmt,data)
    
    def __write_struct(self,data):
        if self.compress:
            chunk = self.__compressor.compress(data)
            if chunk:
                self.fobj.write(chunk)
                self.__md5_update(chunk)
        else:
            self.fobj.write(data)
            self.__md5_update(data)

    def __md5_update(self,data):
        self.__md5.update(data)

    def __check_format(self,name,seq,qual):
        if '0' in seq or '1' in seq or '2' in seq or '3' in seq:
            self.colorspace = True
            self.cs_include_prefix = seq[0].upper() in 'ATCG'
        

def _read_fastq(fobjs):
    try:
        i = 0
        while True:
            name = fobjs[i].next()
            if not name:
                return
            seq = fobjs[i].next()
            fobjs[i].next()
            qual = fobjs[i].next()
            yield (name,seq,qual)
            
            if len(fobjs) > 1:
                i += 1
                if i > len(fobjs):
                    i=0

    except:
        return


def bfq_decode(fname):
    bfq = BFQ(fname)

    for name,seqs,quals in bfq:
        for seq,qual in zip(seqs,quals):
            sys.stdout.write('@%s\n%s\n+\n%s\n' % (name,seq,''.join([chr(x+33) for x in qual])))

    bfq.close()


def bfq_encode_fastq(fnames, stdout=False, compress=False, description='', outname=None, force = False, quiet = False):
    fobjs = []
    for fname in fnames:
        if fname == '-':
            if not os.isatty(0):
                fsize = 0
                f = sys.stdin
                quiet = True
            else:
                sys.stderr.write('Missing data on stdin\n')
                sys.exit(1)
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

    if len(fnames) > 1:
        seq_count = len(fnames)
    else:
        if fobjs[0] != sys.stdin:
            name_count = {}
            for name,seq,qual in _read_fastq(fobjs):
                buf.append((name,seq,qual))
                if not name in name_count:
                    name_count[name] = 1
                else:
                    name_count[name] += 1
                
                if len(name_count) > 10:
                    break
            nc = []
            for name in name_count:
                nc.append(name_count[name])
            
            nc.sort()
            seq_count = nc[-1]
        else:
            seq_count = 1
    
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
                if ext.lower() not in ['gz','fastq','fq','txt'] or skipped:
                    base.insert(0,ext)
                    skipped = True
                    
            base.append('bfq')
            outname = '.'.join(base)

            if os.path.exists(outname) and not force:
                sys.stderr.write('%s exists!' % outname)
                sys.exit(1)
            

    bfq = BFQ(outname,mode,description=description,seq_count=seq_count)
    
    last_pos = 0
    diff = .05 * fsize
    
    if not quiet:
        sys.stderr.write('[...................]\r[>')
    
    for name,seq,qual in buf:
        bfq.write_read(name[1:].strip(),seq.strip(),[ord(q)-33 for q in qual.strip()])
    
    for name,seq,qual in _read_fastq(fobjs):
        if fsize:
            pos = teller.tell()
            if pos - last_pos > diff:
                sys.stderr.write('>')
                sys.stderr.flush()
#                sys.stderr.write('%s%%\n' % (pos* 100 / fsize))
                last_pos = pos
        bfq.write_read(name[1:].strip(),seq.strip(),[ord(q)-33 for q in qual.strip()])
    
    bfq.close()
    sys.stderr.write('\r                     \r')

    for fobj in fobjs:
        if fobj != sys.stdin:
            fobj.close()

def bfq_show_description(fnames):
    for fname in fnames:
        bfq = BFQ(fname)
        if len(fnames) > 1:
            sys.stdout.write('[%s]\n' % fname if fname != '-' else 'stdin')
        sys.stdout.write('%s\n' % bfq.description)
    
def bfq_test(fnames):
    if type(fnames) == str:
        fnames = [fnames,]
        
    allgood = True
    for fname in fnames:
        bfq = BFQ(fname)
        sys.stdout.write('%s\t' % fname if fname != '-' else 'stdin')
    
        f = open(fname)
        f.seek(-16,2)
        hash_pos = f.tell()
        known_md5 = f.read(16)
        sys.stdout.write(''.join(['%02x' % ord(b) for b in known_md5]))

        md5 = hashlib.md5()
        f.seek(0)
        buf = ''
    
        bufsize = 16*1024
        pos = 0
        while pos+bufsize < hash_pos:
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
    print """\
Usage: fastqutils bfq {opts} infile1.fastq {infile2.fastq ... }
       fastqutils bfq -d {opts} filename.bfq
       fastqutils bfq -t filename.bfq
       
Options:
    -d              Decode BFQ file to stdout in FASTQ format
    -t              Test file integrity

    -h              Display this message
    
Encoding options:
    -c              Output to stdout
    -f              Force overwriting existing files
    -i description  Include a description of the file
    -noz            Disable compression
    -o fname        Name of the output file
                    (defaults to input filename.bfq)

Decoding options:
    -s              Only display the description

"""
    sys.exit(1)

if __name__ == '__main__':
    action = 'encode'
    stdout = False
    desc = ''
    compress = True
    show = False
    force = False
    outname = None
    
    fnames = []
    
    last = None
    
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
            
        if last == '-o':
            outname = arg
            last = None
        elif last == '-i':
            desc = arg
            last = None
        elif arg in ['-i','-o']:
            last = arg
        elif arg == '-nz':
            compress = False
        elif arg == '-f':
            force = True
        elif arg == '-c':
            stdout = True
        elif arg == '-s':
            show = True
        elif arg == '-t':
            action = 'test'
        elif arg == '-d':
            action = 'decode'
        elif arg == '-' or os.path.exists(arg):
            fnames.append(arg)

    try:
        if action == 'test':
            if not fnames:
                usage()
            bfq_test(fnames)
        elif action == 'decode':
            if not fnames:
                fnames.append('-')
        
            if show:
                bfq_show_description(fnames)
            else:
                bfq_decode(fnames[0])
        elif action == 'encode':
            if not fnames:
                fnames = ['-']
            bfq_encode_fastq(fnames, stdout=stdout, compress=compress, description=desc, outname=outname, force=force)
    except Exception, e:
       sys.stderr.write('%s: %s' % (type(e).__name__,e))
       sys.stderr.write('\n')
       sys.exit(1)
