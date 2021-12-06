#!/usr/bin/env python

import sys
import os
import re
import pysam
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

def index_frag_mc(inbam):
    with pysam.AlignmentFile(inbam, mode='rb') as fh:
        for read1, read2 in find_read_pair(fh):
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

                    #evaluate overlap
                    if (read1end < read2start):
                        ##use N to fill gap between R1 and R2
                        ##to do check softclip
                        fillgapN = "N"*(read2start - read1end - 1)
                        fxmstring = f'{fragstartXM}{fillgapN}{fragmergeXM}'
                        for fragpos, xm_base in enumerate(fxmstring):
                            if xm_base in 'Z':
                                mC += 1
                            elif xm_base in 'z':
                                umC += 1
                            else:
                                pass
                        totalC = mC + umC
                        if (totalC == 0):
                            fracC = 0
                        else:
                            fracC = '{0:.4}'.format(mC/totalC)
                        pfline = f'{readname}\t{read1ref}\t{fleftmost}\t{frightmost}\t{len(fxmstring)}\t{fracC}\t{fxmstring}'
                        print(pfline)
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

                        for fragpos, xm_base in	enumerate(fxmstring):
                            if xm_base in 'Z':
                                mC += 1
                            elif xm_base in 'z':
                                umC += 1
                            else:
                                pass
                        totalC = mC + umC
                        if (totalC == 0):
                            fracC = 0
                        else:
                            fracC = '{0:.4}'.format(mC/totalC)

                        pfline = f'{readname}\t{read1ref}\t{fleftmost}\t{frightmost}\t{len(fxmstring)}\t{fracC}\t{fxmstring}'
                        print(pfline)

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
                        for fragpos, xm_base in	enumerate(fxmstring):
                            if xm_base in 'Z':
                                mC += 1
                            elif xm_base in 'z':
                                umC += 1
                            else:
                                pass
                        totalC = mC + umC
                        if (totalC == 0):
                            fracC = 0
                        else:
                            fracC = '{0:.4}'.format(mC/totalC)
                        pfline = f'{readname}\t{read1ref}\t{fleftmost}\t{frightmost}\t{len(fxmstring)}\t{fracC}\t{fxmstring}'
                        print(pfline)


if __name__ == '__main__':
    index_frag_mc(sys.argv[1])
