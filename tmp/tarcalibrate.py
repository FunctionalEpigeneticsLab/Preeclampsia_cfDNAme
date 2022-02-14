#!/usr/bin/env python

import sys
import os

def target_calib(inbed, outbed):
    #inbed needs to be sorted by genomic coordinates
    memdict = dict()

    with open(inbed, 'r') as fh:
        for line in fh:
            linfo = line.strip().split('\t')
            (tchrom, tstart, tend, tindex, tname, locstart, locend) = (linfo[0], linfo[1], linfo[2], linfo[3], linfo[4], linfo[6], linfo[7])
            if tindex in memdict:
                (schrom, sstart, send, sname, slocstart, slocend) = memdict[tindex].split(':')
                memval = f'{schrom}:{sstart}:{send}:{sname}:{slocstart}:{locend}'
                memdict[tindex] = memval
            else:
                memval = f'{tchrom}:{tstart}:{tend}:{tindex}:{tname}:{locstart}:{locend}'
                memkey =
                memdict[tindex] = memval

    with open(outbed, 'w') as fo:
        for myindex in memdict:
            (schrom, sstart, send, sname, slocstart, locend) = memdict[myindex].split(':')
            pfline = f'{schrom}\t{sstart}\t{send}\t{myindex}\t{sname}\t{slocstart}\t{locend}\n'
            fo.write(pfline)

if __name__ == '__main__':
    target_calib(sys.argv[1], sys.argv[2])
