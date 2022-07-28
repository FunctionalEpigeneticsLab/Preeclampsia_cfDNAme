#!/usr/bin/env python

import sys
import os
import gzip
import subprocess
from collections import defaultdict

def format_extract(fraginfofh, outfh):
    tmpfh = f'{outfh}.tmp'
    with open(tmpfh, 'w') as fo:
        with open(fraginfofh, 'r') as fh:
            for line in fh:
                tarhapdict = defaultdict(int)
                cstringl = []
                cposl = []
                (chrom, fleft, fright, flen, mcrate, mpos, umpos, fname, fxmstring, tarchr, tarstart, tarend, tartype, overlaps) = line.strip().split('\t')
                if tarchr != "." and mpos != "NULL" and umpos != "NULL":
                    if mpos != "NULL":
                        mposlist = mpos.split(",")
                        for Z in mposlist:
                            #if int(Z) > int(tarstart) and int(Z) < int(tarend):
                            if int(Z) > (int(tarstart)+20) and int(Z) < (int(tarend)-20):
                                tarhapdict[Z] = "C"
                            else:
                                pass

                    if umpos != "NULL":
                        umposlist = umpos.split(",")
                        for z in umposlist:
                            if int(z) > (int(tarstart)+20) and int(z) < (int(tarend)-20):
                                tarhapdict[z] = "T"
                            else:
                                pass

                    for cpos, cstring in sorted(tarhapdict.items(), key=lambda kv: kv[0]):
                        cposl.append(str(cpos))
                        cstringl.append(cstring)

                    if cstringl:
                        cstringf = ''.join(cstringl)
                        cposf = ','.join(cposl)
                        pfline = f'{tarchr}:{tarstart}-{tarend}\t{cstringf}\t1\t{cposf}\t{tartype}\n'
                        fo.write(pfline)

    recordict = defaultdict()

    with open(outfh, 'w') as fo1:
        with open(tmpfh, 'r') as fh1:
            for line1 in fh1:
                if line1 in recordict:
                    recordict[line1] += 1
                else:
                    recordict[line1] = 1

        for record, count in sorted(recordict.items(), key=lambda item: item[1]):
            (tarid, cstringf, fixcnt, cposf, tartype) = record.strip().split("\t")
            pfline = f'{tarid}\t{cstringf}\t{count}\t{cposf}\t{tartype}\n'
            fo1.write(pfline)

    shcmd = f'rm {tmpfh}'
    subprocess.call(shcmd, shell=True)


if __name__ == '__main__':
    format_extract(sys.argv[1], sys.argv[2])
