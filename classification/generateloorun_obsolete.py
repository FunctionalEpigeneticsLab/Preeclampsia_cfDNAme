#!/usr/bin/env python

import sys
import os
import numpy as np

def mypowerfun(x):
    return(10^x)

def write_cmd(runfhprefix, sampleinfo, inputdir, flagindexfh, normalization, outmatprefix, lowvarfilter, outprefix):
    #seqlist0 = list(np.arange(-3,0.5,0.05))
    exescript = "callooeval.R"
    runoutfh1 = f'{runfhprefix}.p1.pbs'
    runoutfh2 = f'{runfhprefix}.p2.pbs'
    alphalist1 = [0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.12, 0.15, 0.18, 0.2]
    alphalist2 = [0.035, 0.045, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1]
    lambdalist  = [0.01000, 0.01259, 0.01585, 0.01995, 0.02512, 0.03981, 0.05012, 0.06310, 0.07943, 0.10000, 0.15849, 0.19953, 0.25119, 0.39811, 0.50119, 0.63096, 0.79433, 1.00000]

    with open(runoutfh1, 'w') as fo:
        fo.write("#!/bin/bash -l\n\n#PBS -l walltime=48:00:00\n#PBS -l nodes=1:ppn=2\n#PBS -l pmem=5gb\n#PBS -N split1loop1\n\n")
        for i in alphalist1:
            for j in lambdalist:
                myalpha = i
                mylambda = j
                pfline1 = f'Rscript {exescript} {sampleinfo} {inputdir} {flagindexfh} {normalization} {outmatprefix}.p1.tsv {lowvarfilter} {myalpha} {mylambda} {outprefix}.glm.loo.assess.alpha.{myalpha}.lambda.{mylambda}.pred.out {outprefix}.glm.loo.assess.alpha.{myalpha}.lambda.{mylambda}.coef.out {outprefix}.glm.loo.assess.alpha.{myalpha}.lambda.{mylambda}.perf.pdf\n\n'
                fo.write(pfline1)

    with open(runoutfh2, 'w') as fo:
        fo.write("#!/bin/bash -l\n\n#PBS -l walltime=48:00:00\n#PBS -l nodes=1:ppn=2\n#PBS -l pmem=5gb\n#PBS -N split1loop2\n\n")
        for i in alphalist2:
            for j in lambdalist:
                myalpha = i
                mylambda = j
                pfline1 = f'Rscript {exescript} {sampleinfo} {inputdir} {flagindexfh} {normalization} {outmatprefix}.p2.tsv {lowvarfilter} {myalpha} {mylambda} {outprefix}.glm.loo.assess.alpha.{myalpha}.lambda.{mylambda}.pred.out {outprefix}.glm.loo.assess.alpha.{myalpha}.lambda.{mylambda}.coef.out {outprefix}.glm.loo.assess.alpha.{myalpha}.lambda.{mylambda}.perf.pdf\n\n'
                fo.write(pfline1)

if __name__ == '__main__':
    write_cmd(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8])
