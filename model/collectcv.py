#!/usr/bin/env python

import sys
import os
import statistics

def summarize_cv(workdir, metricsumout):
    with open(metricsumout, 'w') as fo:
        fo.write('Output\tmeanDEV\tsdDEV\tmeanMSE\tsdMSE\tmeanAUC\tsdAUC\tmeanThresh1\tsdThresh1\tmeanThresh2\tsdThresh2\tmeanSpec1\tsdSpec1\tmeanSens1\tsdSens1\tmeanSpec2\tsdSpec2\tmeanSens2\tsdSens2\tmeanAcc1\tmeanPrec1\tmeanRecall1\tmeanAcc2\tmeanPrec2\tmeanRecall2\n')
        for msefile in os.listdir(workdir):
            devlist = []
            mselist = []
            auclist = []
            thresh1list = []
            thresh2list = []
            spec1list = []
            sens1list = []
            spec2list = []
            sens2list = []
            acc1list = []
            prec1list = []
            recall1list = []
            acc2list = []
            prec2list = []
            recall2list = []
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
                    print(aucfile)
                    for l1 in fh1:
                        if l1.startswith('AUC'):
                            curauc = l1.strip().split(': ')[1]
                            auclist.append(curauc)
                        elif l1.startswith('Youden'):
                            (curmet, curthresh1, curspec1, cursens1, curacc1, curprec1, currecall1) = l1.strip().split('\t')
                            thresh1list.append(curthresh1)
                            spec1list.append(curspec1)
                            sens1list.append(cursens1)
                            acc1list.append(curacc1)
                            if (curprec1 != "NA"):
                                prec1list.append(curprec1)
                            recall1list.append(currecall1)
                        elif l1.startswith('Topleft'):
                            (curmet, curthresh2, curspec2, cursens2, curacc2, curprec2, currecall2) = l1.strip().split('\t')
                            thresh2list.append(curthresh2)
                            spec2list.append(curspec2)
                            sens2list.append(cursens2)
                            acc2list.append(curacc2)
                            if (curprec2 != "NA"):
                                prec2list.append(curprec2)
                            recall2list.append(currecall2)


                devlist = list(map(float, devlist))
                mselist = list(map(float, mselist))
                auclist = list(map(float, auclist))
                thresh1list = list(map(float, thresh1list))
                thresh2list = list(map(float, thresh2list))
                spec1list = list(map(float, spec1list))
                sens1list = list(map(float, sens1list))
                spec2list = list(map(float, spec2list))
                sens2list = list(map(float, sens2list))
                acc1list = list(map(float, acc1list))
                prec1list = list(map(float, prec1list))
                recall1list = list(map(float, recall1list))
                acc2list = list(map(float, acc2list))
                prec2list = list(map(float, prec2list))
                recall2list = list(map(float, recall2list))
                #print(auclist)
                meandev = statistics.mean(devlist)
                sddev = statistics.stdev(devlist)
                meanmse = statistics.mean(mselist)
                sdmse = statistics.stdev(mselist)
                meanauc = statistics.mean(auclist)
                sdauc = statistics.stdev(auclist)
                meanthresh1 = statistics.mean(thresh1list)
                sdthresh1 = statistics.stdev(thresh1list)
                meanthresh2 = statistics.mean(thresh2list)
                sdthresh2 = statistics.stdev(thresh2list)
                meanspec1 = statistics.mean(spec1list)
                sdspec1 = statistics.stdev(spec1list)
                meansens1 = statistics.mean(sens1list)
                sdsens1 = statistics.stdev(sens1list)
                meanacc1 = statistics.mean(acc1list)
                sdacc1 = statistics.stdev(acc1list)
                meanprec1 = statistics.mean(prec1list)
                sdprec1 = statistics.stdev(prec1list)
                meanrecall1 = statistics.mean(recall1list)
                sdrecall1 = statistics.stdev(recall1list)

                meanspec2 = statistics.mean(spec2list)
                sdspec2 = statistics.stdev(spec2list)
                meansens2 = statistics.mean(sens2list)
                sdsens2 = statistics.stdev(sens2list)
                meanacc2 = statistics.mean(acc2list)
                sdacc2 = statistics.stdev(acc2list)
                meanprec2 = statistics.mean(prec2list)
                sdprec2 = statistics.stdev(prec2list)
                meanrecall2 = statistics.mean(recall2list)
                sdrecall2 = statistics.stdev(recall2list)

                pfline = f'{paramprefix}\t{meandev}\t{sddev}\t{meanmse}\t{sdmse}\t{meanauc}\t{sdauc}\t{meanthresh1}\t{sdthresh1}\t{meanthresh2}\t{sdthresh2}\t{meanspec1}\t{sdspec1}\t{meansens1}\t{sdsens1}\t{meanspec2}\t{sdspec2}\t{meansens2}\t{sdsens2}\t{meanacc1}\t{meanprec1}\t{meanrecall1}\t{meanacc2}\t{meanprec2}\t{meanrecall2}\n'
                fo.write(pfline)


if __name__ == '__main__':
    summarize_cv(sys.argv[1], sys.argv[2])
