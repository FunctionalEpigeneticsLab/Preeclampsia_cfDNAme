#!/usr/bin/env python

import sys
import os
import re
import subprocess
import gzip
from collections import defaultdict
import statistics

def mergerepcount(samplemetafh, countdir, outdir):
    with open(samplemetafh, 'r') as fh0:
        for l0 in fh0:
            if not l0.startswith('SubjectID'):
                l0info = l0.strip().split('\t')
                (subjid, samids) = (l0info[0], l0info[1])
                subjoutfh = f'{outdir}/{subjid}.ind.mavg.count.merge.tsv'

                samlist = samids.split(':')
                if len(samlist) == 1:
                    samid = samids
                    saminfh = f'{countdir}/{samid}.ind.mavg.count.tsv'
                    shcmd0 = f'cp {saminfh} {subjoutfh}'
                    print(shcmd0)
                    subprocess.call(shcmd0, shell=True)

                else:
                    mergedict = defaultdict()
                    (totm, totum) = (0, 0)
                    print(subjid)
                    for i in range(0,len(samlist)):
                        samid = samlist[i]

                        saminfh = f'{countdir}/{samid}.sorted_by_name.deduplicated.bismark.cov.sorted.ontar.bed.gz'
                        with gzip.open(saminfh, 'r') as subfh:
                            print(f'processing {saminfh}')
                            for subl in subfh:
                                subl = subl.decode('utf-8')
                                if not subl.startswith('Chromosome'):
                                    (capchr,capstart,capend,capindex,capname,chrom,cstart,cend,meper,mcnt,umcnt,overbase) = subl.strip().split('\t')
                                    if mcnt != ".":
                                        mcnt = int(mcnt)
                                    else:
                                        mcnt = 0
                                    if umcnt != ".":
                                        umcnt = int(umcnt)
                                    else:
                                        umcnt = 0
                                    mergekey = f'{capindex}:{capchr}:{capstart}:{capend}:{capname}:{chrom}:{cstart}'

                                    if mergekey in mergedict:
                                        smcnt,sumcnt = mergedict[mergekey].split(':')
                                        mergemcnt = int(smcnt) + mcnt
                                        mergeumcnt = int(sumcnt) + umcnt
                                        subval = f'{mergemcnt}:{mergeumcnt}'
                                        mergedict[mergekey] = subval
                                    else:
                                        subval = f'{mcnt}:{umcnt}'
                                        mergedict[mergekey] = subval

                    mcounter = defaultdict()
                    for mergekey in mergedict:
                        capindex, capchr, capstart, capend, capname, chrom, cstart = mergekey.split(':')
                        mcnt, umcnt = mergedict[mergekey].split(':')
                        if (int(mcnt)+int(umcnt)) == 0:
                            mper = 0
                        else:
                            mper = 100*int(mcnt)/(int(mcnt)+int(umcnt))
                        keyval = f'{capindex}:{capchr}:{capstart}:{capend}:{capname}'

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


                    with open(subjoutfh, 'w') as fo:
                        fo.write("Chromosome\tStart\tEnd\tIndex\tMePer\tProbe\n")
                        for keyval, storeval in sorted(mcounter.items(), key=lambda item: int(item[0].split(':')[0])):
                            (capindex, chrom, capstart, capend, capsite) = keyval.split(':')
                            meanmper = statistics.mean(storeval)
                            pfline = f'{chrom}\t{capstart}\t{capend}\t{capindex}\t{meanmper}\t{capsite}\n'
                            fo.write(pfline)


if __name__ == '__main__':
    mergerepcount(sys.argv[1], sys.argv[2], sys.argv[3])
