#!/usr/bin/env python

import sys
import os
import statistics

def summarize_cv(workdir, metricsumout):
    devlist = []
    mselist = []
    auclist = []
    with open(metricsumout, 'w') as fo:
        for msefile in os.listdir(workdir):
            if msefile.endswith('devmse.out'):
                paramprefix = msefile.rsplit('.',2)[0]
                aucfile = f'{paramprefix}.predres.out'
                with open(msefile, 'r') as fh0:
                    for l0 in fh0:
                        if l0.startswith('DEV'):
                            curdev = l0.strip().split('\t')[1:]
                            devlist = devlist + curdev
                        elif l0.startswith('MSE'):
                            curmse = l0.strip().split('\t')[1:]
                            mselist = mselist + curmse

                with open(aucfile, 'r') as fh1:
                    for l1 in fh1:
                        if l1.startswith('AUC'):
                            curauc = l1.strip().split(': ')[1]
                            auclist.append(curauc)


                devlist = list(map(float, devlist))
                mselist = list(map(float, mselist))
                auclist = list(map(float, auclist))
                meandev = statistics.mean(devlist)
                sddev = statistics.stdev(devlist)
                meanmse = statistics.mean(mselist)
                sdmse = statistics.stdev(mselist)
                meanauc = statistics.mean(auclist)
                sdauc = statistics.stdev(auclist)
                pfline = f'{msefile}\t{meandev}\t{sddev}\t{meanmse}\t{sdmse}\t{meanauc}\t{sdauc}\n'
                fo.write(pfline)


if __name__ == '__main__':
    summarize_cv(sys.argv[1], sys.argv[2])
