#!/usr/bin/env python

import sys
import os
import re
import subprocess
from collections import defaultdict

def mergeseqmetainfo(codefh, metafh, seqdepfh, bsconversion, outmergefh):
    '''
       merge replicates; fetch sequencing depth, bisulfite conversion information
       1. flag BS conversion failed samples
       2. merge replicates sequencing depth info
       3. fetch meta info and flag failed subjects
    '''
    bscdict = defaultdict()
    seqdict = defaultdict()
    codedict = defaultdict()
    metadict = defaultdict()

    with open(bsconversion, 'r') as fh0:
        for l0 in fh0:
            bscdict[l0.strip()] = 1

    with open(seqdepfh, 'r') as fh1:
        for l1 in fh1:
            if not l1.startswith('sample'):
                (seqsamid, totalbp, meandep, qt3dep, meddep, qt1dep, base15) = l1.strip().split("\t")
                seqdict[seqsamid] = f'{totalbp}:{meandep}:{qt3dep}:{meddep}:{qt1dep}:{base15}'

    with open(codefh, 'r') as fh2:
        for l2 in fh2:
            (subjid, sampleid) = l2.strip().split('\t')
            if subjid in codedict:
                if sampleid in bscdict:
                    stval = f'{codedict[subjid]}:{sampleid}_FBS'
                else:
                    stval = f'{codedict[subjid]}:{sampleid}'
                codedict[subjid] = stval
            else:
                if sampleid in bscdict:
                    stval = f'{sampleid}_FBS'
                else:
                    stval = f'{sampleid}'
                codedict[subjid] = stval

    with open(metafh, 'r') as fh3:
        headerline = fh3.readline().rstrip()
        headerl = headerline.split('\t')
        studyid_pos = int(headerl.index('StudyID'))

        if re.search(r'Study.Code.NIPT2', headerline):
            msubid_pos = int(headerl.index('Study.Code.NIPT2'))
        elif re.search(r'Study.Code.NIPT', headerline):
            msubid_pos = int(headerl.index('Study.Code.NIPT'))
        elif re.search(r'Study.Code.Tissue', headerline):
            msubid_pos = int(headerl.index('Study.Code.Tissue'))
        elif re.search(r'Study.Code.B2', headerline):
            msubid_pos = int(headerl.index('Study.Code.B2'))
        else:
            raise ValueError("sample id field not found")

        pheno0_pos = int(headerl.index('Case.Ctrl'))
        pheno1_pos = int(headerl.index('Case.Control'))

        if re.search(r'GA.NIPT.Days', headerline):
            gad_pos = int(headerl.index('GA.NIPT.Days'))
        elif re.search(r'GA.NIPT2.Days', headerline):
            gad_pos = int(headerl.index('GA.NIPT2.Days'))
        elif re.search(r'GAD.days', headerline):
            gad_pos = int(headerl.index('GAD.days'))
        else:
            raise ValueError("gestational age field not found")

        for l3 in fh3:
            if not l3.startswith('StudyID'):
                l3info = l3.strip().split('\t')
                msubid = l3info[msubid_pos]
                msval = f'{l3info[pheno0_pos]}:{l3info[pheno1_pos]}:{l3info[gad_pos]}'
                metadict[msubid] = msval


    with open(outmergefh, 'w') as fo:
        fo.write("SubjectID\tPhenotype\tDescription\tGA\tSampleID\tMMeanDep\tMMedianDep\tBSFlag\n")
        for ssid, stval in sorted(codedict.items(), key=lambda item: item[0]):
            samlist = stval.split(":")
            (mertotalbp, mermeandep, merqt3dep, mermeddep, merqt1dep, merbase15) = (0, 0, 0, 0, 0, 0)

            (pheno0, pheno1, gad) = metadict[ssid].split(':')
            if len(samlist) > 1:
                BSflag = 'BSPass'
                for k in range(0,len(samlist)):
                    findsamid = samlist[k]
                    if '_FBS' in findsamid:
                        #skip merging seq depth if the sample failed BS conversion
                        BSflag = 'BSFailed'
                    else:
                        (totalbp, meandep, qt3dep, meddep, qt1dep, base15) = seqdict[findsamid].split(':')
                        mermeandep = mermeandep + float(meandep)
                        mermeddep = mermeddep + float(meddep)

            else:
                findsamid = stval.split('_FBS')[0]
                (mertotalbp, mermeandep, merqt3dep, mermeddep, merqt1dep, merbase15) = seqdict[findsamid].split(':')
                if '_FBS' in stval:
                    BSflag = "BSFailed"
                else:
                    BSflag = "BSPass"

            pfline = f'{ssid}\t{pheno0}\t{pheno1}\t{gad}\t{stval}\t{mermeandep}\t{mermeddep}\t{BSflag}\n'
            fo.write(pfline)

def mergerepcount(outmergefh, countdir, outdir, depthresh):
    depthresh = float(depthresh)

    with open(outmergefh, 'r') as fh0:
        for l0 in fh0:
            if not l0.startswith('SubjectID'):
                l0info = l0.strip().split('\t')
                (subjid, samids, meandep, bsflag) = (l0info[0], l0info[4], l0info[5], l0info[7])
                subjoutfh = f'{outdir}/{subjid}.mavg.count.merge.tsv'
                if float(meandep) >= depthresh:
                    samlist = samids.split(':')
                    if len(samlist) == 1:
                        samid = samids
                        if bsflag != "BSFailed":
                            saminfh = f'{countdir}/{samid}.mavg.count.tsv'
                            shcmd0 = f'cp {saminfh} {subjoutfh}'
                            subprocess.call(shcmd0, shell=True)

                        else:
                            pass
                    else:
                        mergedict = defaultdict()
                        (totm, totum) = (0, 0)
                        for i in range(0,len(samlist)):
                            samid = samlist[i]
                            if 'FBS' in samid:
                                pass
                            else:
                                saminfh = f'{countdir}/{samid}.mavg.count.tsv'
                                with open(saminfh, 'r') as subfh:
                                    for subl in subfh:
                                        if not subl.startswith('Chromosome'):
                                            (Chromosome,Start,End,subIndex,Methylated,Unmethylated,Probe) = subl.strip().split('\t')

                                            if subIndex in mergedict:
                                                ssubval = mergedict[subIndex].split(':')
                                                (totm, totum) = (ssubval[0], ssubval[1])
                                                totm = int(totm) + int(Methylated)
                                                totum = int(totum) + int(Unmethylated)
                                                subval = f'{totm}:{totum}:{Chromosome}:{Start}:{End}:{Probe}'
                                                mergedict[subIndex] = subval
                                            else:
                                                subval = f'{Methylated}:{Unmethylated}:{Chromosome}:{Start}:{End}:{Probe}'
                                                mergedict[subIndex] = subval

                        with open(subjoutfh, 'w') as fo:
                            fo.write("Chromosome\tStart\tEnd\tIndex\tMethylated\tUnmethylated\tProbe\n")
                            for subIndex, subval in sorted(mergedict.items(), key=lambda item: int(item[0])):
                                (totm,totum,Chromosome,Start,End,Probe) = subval.split(':')
                                pfline = f'{Chromosome}\t{Start}\t{End}\t{subIndex}\t{totm}\t{totum}\t{Probe}\n'
                                fo.write(pfline)

def mergeproc(codefh, metafh, seqdepfh, bsconversion, outmergefh, countdir, outdir, depthresh):
    mergeseqmetainfo(codefh, metafh, seqdepfh, bsconversion, outmergefh)
    mergerepcount(outmergefh, countdir, outdir, depthresh)

if __name__ == '__main__':
    mergeproc(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8])
