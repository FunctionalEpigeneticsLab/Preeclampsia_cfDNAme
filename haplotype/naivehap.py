#!/usr/bin/env python

import sys
import os
import re
import pysam
import gzip
import subprocess
from collections import defaultdict

'''
Naive try: without breaking long target and predefine regions
trial format: seqID\tchr\tstart\tend\tXMstring

level1 iterator
per read/fragment: index pos and mC
 - For R1 R2 overlap, processing concensus call

level2 iterator
per target: pile pos

level3 iterator
all target
'''

def check_fexist(f):
    if os.path.isfile(f) == True:
        return(1)
    else:
        return(0)

#def get_target(targetfh):
#    tardict = dict()
#    tartype = None
#    with open(targetfh, 'r') as fh:
#        for line in fh:
#            (chrom, start, end, tarname) = line.strip().split('\t')
#            tarkey = f'{chrom}:{start}:{end}'
#            if ('cg' in tarname):
#                tartype = "Blood"
#            elif ('PET' in tarname):
#                tartype = "Placenta"
#            tardict[tarkey] = tartype

#    return(tardict)

#def target_generator(targetfh,infragchr):
#    with open(targetfh, 'r') as fh:
#        for ita in fh:
#            chrom = ita.strip().split('\t')[0]
#            if infragchr == chrom:
#                yield ita

def find_read_pair(inbam, mapregion=None):
    #initialize two-level dict
    matedict = defaultdict(lambda: defaultdict(None))
    for read in inbam.fetch(region=mapregion):
        if not read.is_proper_pair or read.is_secondary or read.is_duplicate or read.is_supplementary or read.is_unmapped or read.mapping_quality < 10:
            continue
        matename = read.query_name
        if matename not in matedict:
            if read.is_read1:
                matedict[matename][0] = read
            else:
                matedict[matename][1] = read
        else:
            if read.is_read1:
                yield read, matedict[matename][1]
            else:
                yield matedict[matename][0], read
            del matedict[matename]

def resolvedeletion(read1, read2):
    ###todo
    ##in case of deletion: deleted bases are not coded in tag
    ##1. overlapped deletion; 2. deletion on one read
    read1cigar = read1.cigarstring
    read2cigar = read2.cigarstring

def resolveinsertion(read1, read2):
    ###todo
    ##
    read1cigar = read1.cigarstring
    read2cigar = read2.cigarstring

def index_frag_mc(inbam, outfh):
    buffer_lines = 100000
    #regions = get_target(targetfh)

    with gzip.open(outfh, 'wt') as fo:
        with pysam.AlignmentFile(inbam, mode='rb') as fh:
            buffchunk = ""
            buffl = 0

            for read1, read2 in find_read_pair(fh):
                #initialize
                (readname, read1ref, fxmstring) = (None, None, None)
                (fleftmost, frightmost, fracC) = (-1, -1, -1)
                read1ref = read1.reference_name
                read2ref = read2.reference_name
                #skip reads with deletion temporarily
                read1cigar = read1.cigarstring
                read2cigar = read2.cigarstring
                if ('D' in read1cigar) or ('D' in read2cigar) or ('I' in read1cigar) or ('I' in read2cigar):
                    continue
                else:
                    if (read1ref == read2ref):
                        read1start = read1.reference_start + 1
                        read2start = read2.reference_start + 1
                        read1end = read1.reference_end
                        read2end = read2.reference_end
                        ftlen = abs(read1.template_length)
                        #pysam reference_start: 0-based; match bam file to 1-based
                        fleftmost = min(read1start, read2start)
                        frightmost = fleftmost + ftlen - 1

                        readname = read1.query_name
                        read1xm = read1.get_tag('XM')
                        read2xm = read2.get_tag('XM')
                        frag_dict = defaultdict(int)

                        #evalute top/bottom strand
                        fstrand = None
                        read1xg = read1.get_tag('XG')
                        read2xg = read2.get_tag('XG')
                        if (read1xg == read2xg == "GA"):
                            fstrand = "BOT"
                        elif (read1xg == read2xg == "CT"):
                            fstrand = "TOP"
                        else:
                            print(f'error in strand {readname}')

                        #evalute letfmost read
                        if fleftmost == read1start:
                            lreadstart = read1start
                            matestart = read2start
                            fragstartXM = read1xm
                            fragmergeXM = read2xm

                        else:
                            lreadstart = read2start
                            matestart = read1start
                            fragstartXM = read2xm
                            fragmergeXM = read1xm

                        #evalute righmost read
                        if frightmost == read1end:
                            lreadend = read2end
                            mateend = read1end
                        else:
                            lreadend = read1end
                            mateend = read2end

                        (mC, umC, totalC) = (0, 0, 0)
                        Zmarkpos = []
                        zmarkpos = []

                        #evaluate overlap
                        if (read1end < read2start):
                            ##use N to fill gap between R1 and R2
                            ##to do check softclip
                            fillgapN = "N"*(read2start - read1end - 1)
                            fxmstring = f'{fragstartXM}{fillgapN}{fragmergeXM}'
                            for fragpos, xm_base in enumerate(fxmstring):
                                if fstrand == "BOT":
                                    fragpos = fragpos - 1

                                if xm_base in 'Z':
                                    mC += 1
                                    Zmarkpos.append(str(lreadstart + fragpos))
                                elif xm_base in 'z':
                                    umC += 1
                                    zmarkpos.append(str(lreadstart + fragpos))
                                else:
                                    pass

                        elif (read1start <= read2start and read1end >= read2end) or (read2start <= read1start and read2end >= read1end):
                            #full overlap; one read covers the mate entirely
                            #code unmatched mC call from R1 and R2 as 'D'
                            overlapxm = ""
                            overlapSpos = matestart - lreadstart
                            overlapEpos = lreadend - matestart + overlapSpos + 1
                            leftxmstring = f'{fragstartXM[0:(overlapSpos)]}'
                            rightxmstring = f'{fragstartXM[overlapEpos:]}'
                            for i in range(overlapSpos,overlapEpos):
                                if (fragstartXM[i] == fragmergeXM[i-overlapSpos]):
                                    overlapxm += fragstartXM[i]
                                else:
                                    overlapxm += 'D'
                            fxmstring = f'{leftxmstring}{overlapxm}{rightxmstring}'

                            for fragpos, xm_base in enumerate(fxmstring):
                                if fstrand == "BOT":
                                    fragpos = fragpos -	1

                                if xm_base in 'Z':
                                    mC += 1
                                    Zmarkpos.append(str(lreadstart + fragpos))
                                elif xm_base in 'z':
                                    umC += 1
                                    zmarkpos.append(str(lreadstart + fragpos))
                                else:
                                    pass

                        else:
                            #partial overlap; merge overlap here
                            overlapxm = ""
                            overlapSpos = matestart - lreadstart
                            overlapEpos = lreadend - matestart + overlapSpos + 1
                            leftxmstring = f'{fragstartXM[0:(overlapSpos)]}'
                            rightxmstring = f'{fragmergeXM[(overlapEpos-overlapSpos):]}'
                            for i in range(overlapSpos, overlapEpos):
                                if (fragstartXM[i] == fragmergeXM[i-overlapSpos]):
                                    overlapxm += fragstartXM[i]
                                else:
                                    overlapxm += 'D'
                            fxmstring = f'{leftxmstring}{overlapxm}{rightxmstring}'
                            for fragpos, xm_base in enumerate(fxmstring):
                                if fstrand == "BOT":
                                    fragpos = fragpos -	1

                                if xm_base in 'Z':
                                    mC += 1
                                    Zmarkpos.append(str(lreadstart + fragpos))
                                elif xm_base in 'z':
                                    umC += 1
                                    zmarkpos.append(str(lreadstart + fragpos))
                                else:
                                    pass

                        totalC = mC + umC
                        if (totalC == 0):
                            fracC = 0
                        else:
                            fracC = '{0:.4}'.format(mC/totalC)

                        if not Zmarkpos:
                            Zmarkposl = 'NULL'
                        else:
                            Zmarkposl = ','.join(Zmarkpos)

                        if not zmarkpos:
                            zmarkposl = 'NULL'
                        else:
                            zmarkposl = ','.join(zmarkpos)

                        ##fetch target information
                        #(ftarstart, ftarend) = (0, 0)
                        #ftartype = '-'
                        #for ita in regions:
                        #    (chrom, tarstart, tarend) = ita.split(":")

                        #    if read1ref == chrom:
                        #        if (frightmost > int(tarstart) and frightmost < int(tarend)) or (fleftmost > int(tarstart) and frightmost < int(tarend)) or (fleftmost > int(tarstart) and fleftmost < int(tarend)):
                        #            ftarstart = tarstart
                        #            ftarend = tarend
                        #            ftartype = regions[ita]
                        #            break

                        #pfline = f'{read1ref}\t{ftarstart}\t{ftarend}\t{ftartype}\t{fleftmost}\t{frightmost}\t{len(fxmstring)}\t{fracC}\t{Zmarkposl}\t{zmarkposl}\t{readname}\t{fxmstring}\n'
                        pfline = f'{read1ref}\t{fleftmost}\t{frightmost}\t{len(fxmstring)}\t{fracC}\t{Zmarkposl}\t{zmarkposl}\t{readname}\t{fxmstring}\n'

                        buffl += 1
                        if buffl >= buffer_lines:
                            fo.write(buffchunk)
                            buffchunk = ""
                            buffl = 0

                        buffchunk += pfline
            if (buffl > 0):
                fo.write(buffchunk)
                buffchunk = ""

def fetch_target(inbam, indexfragfh, targetfh):
    #index_frag_mc(inbam, indexfragfh)
    fhprefix = indexfragfh.rsplit('.',2)[0]
    tmpfh0 = f'{fhprefix}.tmp'
    #shcmd0 = f'gunzip -c {indexfragfh} > {tmpfh0}'
    #subprocess.call(shcmd0, shell=True)
    tmpfh1 = f'{fhprefix}.target.tmp'
    #shcmd1 = f'bedtools intersect -a {tmpfh0} -b {targetfh} -wao > {tmpfh1}'
    #subprocess.call(shcmd1, shell=True)
    ##shcmd2 = f'rm {tmpfh0}'
    ##subprocess.call(shcmd2, shell=True)

    ###A fragment may overlap adjacent targets, assign the fragment to the target with greater overlap
    ###prerequsite: duplicated entries are in order
    pname = "NULL"
    lcounter = 0
    pdist = 0
    flushline = ""

    outfh = f'{fhprefix}.target.bed'
    try:
        with open(tmpfh1, 'r') as fh, open(outfh, 'w') as fo:
            for line in fh:
                curline = line.strip()
                (readname, interdist) = (line.strip().split("\t")[7], line.strip().split("\t")[13])
                interdist = int(interdist)
                if lcounter == 0:
                    flushline = curline
                    pname = readname
                    pdist = interdist
                else:
                    if readname == pname:
                        if (interdist > pdist):
                            flushline = curline
                            pname = readname
                            pdist = interdist
                        else:
                            pass

                    else:
                        pfline = flushline
                        fo.write(f'{pfline}\n')
                        flushline = curline
                        pname = readname
                        pdist = interdist

                lcounter += 1


    except OSError:
        print('cannot open', tmpfh1)



if __name__ == '__main__':
    fetch_target(sys.argv[1], sys.argv[2], sys.argv[3])
