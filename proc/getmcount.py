#!/usr/bin/env python

from __future__ import division
import sys
import os
import time
import gzip
import subprocess
from collections import defaultdict


def check_fexist(f):
    if os.path.isfile(f) == True:
        return(1)
    else:
        return(0)

def getmavgcountfromcov(samplelist, capindexfh, outdir):

    with open(samplelist, 'r') as fh0:
        for l0 in fh0:
            (infh, fhprefix) = l0.strip().split('\t')
            sortfh = f'{outdir}/{fhprefix}.sorted_by_name.deduplicated.bismark.cov.sorted.bed.gz'
            shcmd0 = f'bedtools sort -i {infh} | gzip - > {sortfh}'
            subprocess.call(shcmd0, shell=True)

            ontarfh = f'{outdir}/{fhprefix}.sorted_by_name.deduplicated.bismark.cov.sorted.ontar.bed.gz'
            shcmd1 = f'bedtools intersect -a {capindexfh} -b {sortfh} -wao | gzip - > {ontarfh}'
            print(shcmd1)
            subprocess.call(shcmd1, shell=True)

            if check_fexist(ontarfh):
                mcounter = defaultdict()
                outfh = f'{outdir}/{fhprefix}.mavg.count.tsv'
                with open(outfh, 'w') as fo:
                    with gzip.open(ontarfh, 'r') as fh1:
                        for l1 in fh1:
                            l1d = l1.decode('utf-8')
                            l1info = l1d.strip().split('\t')

                            (chrom, capstart, capend, capindex, capsite, mcount, umcount) = (l1info[0], l1info[1], l1info[2], l1info[3], l1info[4], l1info[9], l1info[10])

                            if capindex in mcounter:
                                (totmcnt, totumcnt, dchrom, dcapstart, dcapend, dcapsite) = mcounter[capindex].split(':')
                                totmcnt = int(totmcnt)+int(mcount)
                                totumcnt = int(totumcnt)+int(umcount)
                                storeval = f'{mcount}:{umcount}:{dchrom}:{dcapstart}:{dcapend}:{dcapsite}'
                                mcounter[capindex] = storeval

                            else:
                                if mcount == "." and umcount == ".":
                                    storeval = f'0:0:{chrom}:{capstart}:{capend}:{capsite}'
                                else:
                                    storeval = f'{mcount}:{umcount}:{chrom}:{capstart}:{capend}:{capsite}'
                                mcounter[capindex] = storeval


                    for capindex, storeval in sorted(mcounter.items(), key=lambda item: int(item[0][1:])):
                        (totmcnt, totumcnt, chrom, capstart, capend, capsite) = mcounter[capindex].split(':')
                        pfline = f'{chrom}\t{capstart}\t{capend}\t{capindex}\t{totmcnt}\t{totumcnt}\t{capsite}\n'
                        fo.write(pfline)

if __name__ == '__main__':
    getmavgcountfromcov(sys.argv[1], sys.argv[2], sys.argv[3])
