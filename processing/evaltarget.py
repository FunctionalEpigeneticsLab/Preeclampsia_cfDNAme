#!/usr/bin/env python

import sys
import os
import time
import re

'''
   examine problematic capture regions from a set of samples
   - input file: GATK -T DiagnoseTargets output (vcf file)
   - input parameter: a float indicating proportion of samples to define Pass/Fail of a region
   - output file: tab-separated file with Pass/Fail QCFlag
'''

def collectdiaginfo(invcf, outfh, prop):
    prop = float(prop)
    with open(outfh, 'w') as fo:
        fo.write('Chromosome\tStart\tEnd\tQCFlag\tBAD_MATE\tCOVERAGE_GAPS\tEXCESSIVE_COVERAGE\tLOW_COVERAGE\tNO_READS\tPASS\tPOOR_QUALITY\n')
        with open(invcf, 'r') as fh:
            for line in fh:
                (nbadmate, ngap, nexcess, nlow, nnoread, npass, npoorq) = tuple([0]*7)
                QCFlag = "Null"
                if line.startswith("##"):
                    pass
                elif line.startswith("#CHROM"):
                    lineinfo = line.split('\t')
                    try:
                        info_pos = int(lineinfo.index("INFO"))
                        format_pos = int(lineinfo.index("FORMAT"))
                    except ValueError:
                        print("check vcf field infomartion.")
                    nsam = len(lineinfo)-format_pos-1
                    print(f'Total number of samples: {nsam}')
                else:
                    lineinfo = line.split('\t')
                    (chrom, tstart) = (lineinfo[0],int(lineinfo[1])-1)
                    tend = re.split('=|;',lineinfo[info_pos])[1]
                    for i in range(format_pos,len(lineinfo)):
                        #samFT = lineinfo[i].split(':')[0]
                        #count the first feature
                        samFT = re.split(':|;',lineinfo[i])[0]

                        if samFT == "BAD_MATE":
                            nbadmate += 1
                        elif samFT == "COVERAGE_GAPS":
                            ngap += 1
                        elif samFT == "EXCESSIVE_COVERAGE":
                            nexcess += 1
                        elif samFT == "LOW_COVERAGE":
                            nlow += 1
                        elif samFT == "NO_READS":
                            nnoread += 1
                        elif samFT == "PASS":
                            npass += 1
                        elif samFT == "POOR_QUALITY":
                            npoorq += 1
                    #if bad targets present in more than 70% of samples
                    #if (nbadmate+ngap+nlow+nnoread+npoorq) > nsam*0.3:
                    if (nbadmate+ngap+nlow+nnoread+npoorq) > nsam*prop:
                        QCFlag = "Fail"
                    else:
                        QCFlag = "Pass"
                    pfline = f'{chrom}\t{tstart}\t{tend}\t{QCFlag}\t{nbadmate}\t{ngap}\t{nexcess}\t{nlow}\t{nnoread}\t{npass}\t{npoorq}\n'
                    fo.write(pfline)

if __name__ == '__main__':
    collectdiaginfo(sys.argv[1], sys.argv[2], sys.argv[3])
