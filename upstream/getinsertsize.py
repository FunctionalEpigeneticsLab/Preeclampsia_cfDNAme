#!/usr/bin/env python

import sys
import os
import time
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as st
from collections import defaultdict

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


def get_insert_size(inbam):
    fsizedict = defaultdict()
    outfh = f'{inbam}.finsert.tsv'
    outfig = f'{inbam}.finsert.pdf'
    farr = []

    with open(outfh, 'w') as fo:
        with pysam.AlignmentFile(inbam, mode='rb') as fh:
            buffchunk = ""
            buffl = 0
            for read1, read2 in find_read_pair(fh):
                finsert = 0
                readname = read1.query_name
                read1ref = read1.reference_name
                read1start = read1.reference_start
                read2start = read1.next_reference_start

                if read1.is_reverse:
                    read1cigarA = read1.get_cigar_stats()[0]
                    read1mbp = read1cigarA[0]+read1cigarA[2]
                    finsert = read1start + read1mbp - read2start
                else:
                    read2cigarA = read2.get_cigar_stats()[0]
                    read2mbp = read2cigarA[0]+read2cigarA[2]
                    finsert = read2start + read2mbp - read1start
                #pfline = f'{readname}\t{read1ref}\t{read1start}\t{read2start}\t{finsert}\n'
                #fo.write(pfline)
                farr.append(finsert)
                if finsert in fsizedict:
                    fsizedict[finsert] += 1
                else:
                    fsizedict[finsert] = 1

        for record, count in sorted(fsizedict.items(), key=lambda item: item[0]):
            pfline = f'{record}\t{count}\n'
            fo.write(pfline)

        figtitle = f'{os.path.basename(inbam)}.Histrogram'
        #plt.hist(farr, density=True, bins=500, label="InsertSize")
        #mn, mx = plt.xlim()
        #kde_xs = np.linspace(mn,mx,300)
        #kde = st.gaussian_kde(farr)
        #plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

        num_bins = min(max(farr),500)
        fig, ax = plt.subplots()
        sns.histplot(farr, binwidth=1, ax=ax, color="lightblue")
        ax2 = ax.twinx()
        bin_width = (max(farr) - min(farr)) / num_bins
        hist_area = len(farr) * bin_width
        ax2.set_ylim(ymax=ax.get_ylim()[1] / hist_area)
        sns.kdeplot(farr,ax=ax2)
        plt.xlabel("InsertSize")
        plt.tight_layout()
        plt.setp(ax.patches, linewidth=0)
        plt.title(figtitle)
        plt.savefig(outfig)
        plt.close(outfig)


if __name__ == '__main__':
    get_insert_size(sys.argv[1])
