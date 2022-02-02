#!/usr/bin/env python

import sys
import os
import re

def mergeseqmetainfo(codefh, metafh, seqdepfh, bsconversion, outmergefh):
    '''
       merge replicates; fetch sequencing depth, bisulfite conversion information
       1. flag BS conversion failed samples
       2. merge replicates sequencing depth info
       3. fetch meta info and flag failed subjects
    '''
    bscdict = dict()
    seqdict = dict()
    codedict = dict()
    metadict = dict()

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
        else:
            raise ValueError("sample id field not found")

        pheno0_pos = int(headerl.index('Case.Ctrl'))
        pheno1_pos = int(headerl.index('Case.Control'))

        if re.search(r'GA.NIPT1.Days', headerline):
            gad_pos = int(headerl.index('GA.NIPT1.Days'))
        elif re.search(r'GA.NIPT2.Days', headerline):
            gad_pos = int(headerl.index('GA.NIPT2.Days'))
        elif re.search(r'GA.last.pregnancy.days', headerline):
            gad_pos = int(headerl.index('GA.last.pregnancy.days'))
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

if __name__ == '__main__':
    mergeseqmetainfo(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
