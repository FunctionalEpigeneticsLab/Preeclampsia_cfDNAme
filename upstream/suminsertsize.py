#!/usr/bin/env python

import sys
import os
import statistics

def get_insert_metrics(workdir,inmetafh,outfh):
    with open(outfh,'w') as fo:
        fo.write('SubjectID\tPhenotype\tMeanSize\tMedianSize\tModeSize\tSeqFF\tcfDNAme\tLibConcentration\n')
        with open(inmetafh,'r') as fh:
            for line in fh:
                if not line.startswith("SubjectID"):
                    #subfragtotsize = 0
                    #subfragtotcnt = 0
                    subsizelist = []
                    linfo = line.strip().split('\t')
                    (subjid,pheno,sampleids,seqff,cfdname,libcon) = (linfo[0],linfo[1],linfo[4],linfo[11],linfo[15],linfo[16])
                    if ":" in sampleids:
                        samlist = sampleids.split(':')
                        for samid in samlist:
                            curfh = f'{workdir}/{samid}/{samid}.merge.deduplicated.filtered.sorted.bam.finsert.tsv'
                            with open(curfh,'r') as subfh:
                                for subline in subfh:
                                    (fragsize,fragcount) = subline.strip().split('\t')
                                    subsizelist = subsizelist + [int(fragsize)]*int(fragcount)
                    else:
                        curfh = f'{workdir}/{samid}/{samid}.merge.deduplicated.filtered.sorted.bam.finsert.tsv'
                        with open(curfh,'r') as subfh:
                            for subline in subfh:
                                (fragsize,fragcount) = subline.strip().split('\t')
                                subsizelist = subsizelist + [int(fragsize)]*int(fragcount)

                    meansize = statistics.mean(subsizelist)
                    mediansize = statistics.median(subsizelist)
                    modesize = max([p[0] for p in statistics._counts(subsizelist)])
                    pfline = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (subjid,pheno,meansize,mediansize,modesize,seqff,cfdname,libcon)
                    fo.write(pfline)


if __name__ == '__main__':
    get_insert_metrics(sys.argv[1],sys.argv[2],sys.argv[3])
