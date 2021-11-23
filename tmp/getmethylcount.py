#!/usr/bin/env python

import sys
import os
import gzip

def getcount(infh, outfh):
    covdict = dict()
    methdict = dict()

    with open(outfh, 'w') as fo:
        with gzip.open(infh, 'r') as fh:
            for line in fh:
                if not line.startswith('Bismark'):
                    (reads, methmark, chrom, pos, methstate) = line.strip().split('\t')
                    fetchid = '%s:%s' % (chrom, pos)
                    if fetchid in covdict:
                        covdict[fetchid] += 1
                        if (methstate == 'Z'):
                            methdict[fetchid] += 1
                        else:
                            pass
                    else:
                        covdict[fetchid] = 1
                        methdict[fetchid] = 0
                        if (methstate == 'Z'):
                            methdict[fetchid] += 1

        for getid in covdict:
            (gchrom, gpos) = getid.split(':')
            pfline = '%s\t%s\t%s\t%s\n' % (gchrom, gpos, covdict[getid], methdict[getid])
            fo.write(pfline)

if __name__ == '__main__':
    getcount(sys.argv[1], sys.argv[2])
