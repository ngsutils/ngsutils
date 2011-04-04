
import os
import support.ngs_utils,support.localalign
from support.eta import ETA

def read_fastq(fname,quiet=False,eta_callback=None):
    with support.ngs_utils.gzip_opener(fname) as f:
        if not quiet:
            eta = ETA(os.stat(fname).st_size,fileobj=f)
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
                yield (name,seq,qual)
            except Exception, e:
                break
    if not quiet:
        eta.done()
