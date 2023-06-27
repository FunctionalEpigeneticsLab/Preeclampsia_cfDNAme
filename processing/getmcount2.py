#!/usr/bin/env python

from __future__ import division
import sys
import os
import time
import gzip
import subprocess
from collections import defaultdict
import statistics

'''
   extract on-target methylated/unmethylated count from bismark coverage output
   - samplelist: a list of samples to be extracted with one sample per line
   - capindexfh: capture index file
   - outdir: output directory
'''


def check_fexist(f):
    if os.path.isfile(f) == True:
        return(1)
    else:
        return(0)

def getmavgindcountfromcov(samplelist, capindexfh, outdir):

    with open(samplelist, 'r') as fh0:
        for l0 in fh0:
            (infh, fhprefix) = l0.strip().split('\t')
            sortfh = f'{outdir}/{fhprefix}.sorted_by_name.deduplicated.bismark.cov.sorted.bed.gz'
            if check_fexist(sortfh) == 0:
                shcmd0 = f'bedtools sort -i {infh} | gzip - > {sortfh}'
                print(f'processing {infh}')
                subprocess.call(shcmd0, shell=True)

            ontarfh = f'{outdir}/{fhprefix}.sorted_by_name.deduplicated.bismark.cov.sorted.ontar.bed.gz'
            shcmd1 = f'bedtools intersect -a {capindexfh} -b {sortfh} -wao | gzip - > {ontarfh}'
            print("extracting on-target counts")
            subprocess.call(shcmd1, shell=True)

            if check_fexist(ontarfh):
                mcounter = defaultdict()
                outfh = f'{outdir}/{fhprefix}.ind.mavg.count.tsv'
                with open(outfh, 'w') as fo:
                    fo.write("Chromosome\tStart\tEnd\tIndex\tMePer\tProbe\n")
                    with gzip.open(ontarfh, 'r') as fh1:
                        for l1 in fh1:
                            l1d = l1.decode('utf-8')
                            l1info = l1d.strip().split('\t')

                            chrom, capstart, capend, capindex, capsite, mper = l1info[0], l1info[1], l1info[2], l1info[3], l1info[4], l1info[8]
                            keyval = f'{capindex}:{chrom}:{capstart}:{capend}:{capsite}'

                            if keyval in mcounter:
                                mper = float(mper)
                                mcounter[keyval].append(mper)

                            else:
                                if mper == ".":
                                    storeval = [0]
                                else:
                                    mper = float(mper)
                                    storeval = [mper]
                                mcounter[keyval] = storeval


                    #for capindex, storeval in sorted(mcounter.items(), key=lambda item: int(item[0][1:])):
                    for keyval, storeval in sorted(mcounter.items(), key=lambda item: int(item[0].split(':')[0])):
                        (capindex, chrom, capstart, capend, capsite) = keyval.split(':')
                        meanmper = statistics.mean(storeval)
                        pfline = f'{chrom}\t{capstart}\t{capend}\t{capindex}\t{meanmper}\t{capsite}\n'
                        fo.write(pfline)

if __name__ == '__main__':
    getmavgindcountfromcov(sys.argv[1], sys.argv[2], sys.argv[3])
