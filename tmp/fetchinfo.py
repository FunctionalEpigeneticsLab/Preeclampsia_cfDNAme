#!/usr/bin/env python

import sys
import os


def mergecovinfo(removefh, seqcovfh, outmergefh):
    pedict = dict()
    rmdict = dict()

    with open(removefh, 'r') as fh0:
        for l0 in fh0:
            (studyid, gccode, rmcause) = l0.strip().split('\t')
            rmdict[gccode] = studyid

    with open(outmergefh, 'w') as fo:
        with open(seqcovfh, 'r') as fh:
            for line in fh:
                mergccode = 'Null'
                (mertotalbp, mermeandep, merqt3dep, mermeddep, merqt1dep, merbase15) = (0, 0, 0, 0, 0, 0)

                (peid, gccode, totalbp, meandep, qt3dep, meddep, qt1dep, base15) = line.strip().split('\t')
                if gccode in rmdict:
                    continue
                else:
                    if peid in pedict:
                        (expeid, exgccode, extotalbp, exmeandep, exqt3dep, exmeddep, exqt1dep, exbase15) = pedict[peid].split('\t')
                        mergccode = f'{exgccode};{gccode}'
                        mertotalbp = int(extotalbp) + int(totalbp)
                        mermeandep = '{:.2f}'.format(float(exmeandep) + float(meandep))
                        merqt3dep = float(exqt3dep) + float(qt3dep)
                        mermeddep = float(exmeddep) + float(meddep)
                        merqt1dep = float(exqt1dep) + float(qt1dep)
                        merbase15 = '{:.1f}'.format(float(exbase15) + float(base15))
                        storeinfo = f'{peid}\t{mergccode}\t{mertotalbp}\t{mermeandep}\t{merqt3dep}\t{mermeddep}\t{merqt1dep}\t{merbase15}\n'
                        pedict[peid] = storeinfo
                    else:
                        pedict[peid] = line

        for item in pedict:
            fo.write(pedict[item])


def fetchmetainfo(metadata, covmergefh, outfh):
    namedict = dict()

    with open(metadata, 'r') as fh0:
        for l0 in fh0:
            if not l0.startswith('Study'):
                studyid = l0.strip().split('\t')[0]
                studyn1id = f'{studyid}N1'
                namedict[studyn1id] = l0.strip()

    with open(outfh, 'w') as fo:
        fo.write('MergeGCcode\tMeanDep\tMedianDep\tStudyID\tCase.Ctrl\tCase.Control\tMAC\tGAD.weeks\tGA.NIPT.Weeks\tFF\tNIPT.result\tStudy.Code.NIPT\tbisulfite.conv.NIPT1\t%regions.with.less.than.50.NIPT1\tBMI\tProteinuria\tThrombocytopenia\tLiver.transaminase\tSevere.PET\tGestational.DM\tSex.newborn\tSex.newborn2\n')
        with open(covmergefh, 'r') as fh1:
            for l1 in fh1:
                (matchid, mergccode, mertotalbp, mermeandep, merqt3dep, mermeddep, merqt1dep, merbase15) = l1.strip().split('\t')

                if matchid in namedict:
                    pfline = f'{mergccode}\t{mermeandep}\t{mermeddep}\t{namedict[matchid]}\n'
                    fo.write(pfline)


if __name__ == '__main__':
    #mergecovinfo(sys.argv[1], sys.argv[2], sys.argv[3])
    fetchmetainfo(sys.argv[1], sys.argv[2], sys.argv[3])
