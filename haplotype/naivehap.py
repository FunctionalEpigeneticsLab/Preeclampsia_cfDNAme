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

def index_frag_mc(inbam):
    with pysam.AlignmentFile(inbam, mode='rb') as fh:
        for read1, read2 in find_read_pair(fh):
            read1ref = read1.reference_name
            read2ref = read2.reference_name
            if (read1ref == read2ref):
                read1start = read1.reference_start + 1
                read2start = read2.reference_start + 1
                ftlen = abs(read1.template_length)
                #pysam reference_start: 0-based; match bam file to 1-based
                fleftmost = min(read1start, read2start)
                frightmost = fleftmost + ftlen - 1

                readname = read1.query_name
                read1xm = read1.get_tag('XM')
                read2xm = read2.get_tag('XM')
                frag_dict = defaultdict(int)

                if fleftmost == read1start:
                    fragstartXM = read1xm
                    fragmergeXM = read2xm

                else:
                    fragstartXM = read2xm
                    fragmergeXM = read1xm

                #evaluate overlap
                if (read1.reference_end < read2start):
                    ##use N to fill gap between R1 and R2
                    ##to do check softclip
                    read1end = read1.reference_end
                    print(read1end)
                    fillgapN = "N"*(read2start - read1end - 1)
                    fxmstring = f'{fragstartXM}{fillgapN}{fragmergeXM}'
                    pfline = f'{readname}\t{read1ref}\t{fleftmost}\t{frightmost}\t{fxmstring}'
                    print(pfline)
                else:
                    continue
                    #merge overlap here

if __name__ == '__main__':
    index_frag_mc(sys.argv[1])
