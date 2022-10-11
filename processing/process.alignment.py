#!/usr/bin/env python

import sys
import os
import time
import helper

(fastqdir, workdir, runscriptdir) = helper.get_dir_config()
(samtools, samtoolsdir, fastqc, trimtool, bismark, bismarkdedup, bismarkextractor, bismarkcytosine, bowtiedir, java, gatk, picard, addsourcetool, python, bismarkreference, reference, target, targetflank50, targetflank100) = helper.get_run_config()
DNASampleIDs = helper.get_DNA_sample()
(DNAtrimoption) = helper.get_option()
(noticeaccount, creditaccount) = helper.get_userinfo()

def write_sub_script_qsub(taskscriptdir, taskname, taskworkdir, stderrdir):
    '''
    gather sample task together for qsub
    taskname added in case of multiple qsub in one scriptdir
    stderrdir added in case output to terminal or stdout
    '''
    qsub = 'qsub'

    task_qsub = taskscriptdir + '/' + taskname + '.qsub.sh'
    with open(task_qsub, 'w+') as fo:
        fo.write('#!/bin/bash -l\n\n')

        for script in os.listdir(taskscriptdir):
            tasksearch = taskname + '.task.pbs'
            if script.endswith(tasksearch):
                #myqsub = '%s -d %s -e %s -o %s %s\n' % (qsub, taskworkdir, stderrdir, logdir, script)
                myqsub = '%s %s\n' % (qsub, script)
                fo.write(myqsub)
    return(task_qsub)


def write_sub_script_fastqc(SampleIDs, fastqdir, scriptdir, workdir, logdir):
    fastqc_scriptdir = scriptdir + '/1_fastqc'
    fastqc_workdir = workdir + '/1_qc'
    helper.mkdir(fastqc_scriptdir)
    helper.mkdir(fastqc_workdir)

    for sampleid in SampleIDs:

        samplescript = fastqc_scriptdir + '/' + sampleid + '.fastqc.task.pbs'

        with open(samplescript, 'w+') as fo:
            pbsjobtime = '#PBS -l walltime=%s:00:00\n' % (12)
            pbsnode = '#PBS -l nodes=1:ppn=%s\n' % (1)
            pbsuser = '#PBS -M %s\n' % (noticeaccount)
            pbsjobname = '#PBS -N fastqc_%s\n' % (sampleid)
            pbscredit = '#PBS -A %s\n\n' % (creditaccount)
            pbsaddsource = 'source %s\n\n' % (addsourcetool)
            fo.write('#!/bin/bash -l\n\n')
            fo.write(pbsjobtime)
            fo.write(pbsnode)
            fo.write('#PBS -l pmem=5gb\n')
            fo.write(pbsuser)
            fo.write('#PBS -m abe\n')
            fo.write(pbsjobname)
            fo.write(pbscredit)
            fo.write(pbsaddsource)

            for fastqfile in os.listdir(fastqdir):
                if fastqfile.startswith(sampleid) and (fastqfile.endswith('R1.fastq.gz') or fastqfile.endswith('R2.fastq.gz') or fastqfile.endswith('R1_001.fastq.gz') or fastqfile.endswith('R2_001.fastq.gz')):
                    qc_task = '%s -t 2 %s/%s -o %s\n' % (fastqc, fastqdir, fastqfile, fastqc_workdir)
                    fo.write(qc_task)

    write_sub_script_qsub(fastqc_scriptdir, 'fastqc', fastqc_workdir, logdir)


def write_sub_script_DNAtrimming(DNASampleIDs, fastqdir, scriptdir, workdir, logdir, DNAtrimoption):
    trim_scriptdir = scriptdir + '/2_DNAtrim'
    trim_workdir = workdir + '/2_DNAtrim'
    helper.mkdir(trim_scriptdir)
    helper.mkdir(trim_workdir)

    for sampleid in DNASampleIDs:
        samplescript = trim_scriptdir + '/' + sampleid + '.DNAtrim.task.pbs'
        with open(samplescript, 'w') as fo:
            pbsjobtime = '#PBS -l walltime=%s:00:00\n' % (12)
            pbsnode = '#PBS -l nodes=1:ppn=%s\n' % (2)
            pbsuser = '#PBS -M %s\n' % (noticeaccount)
            pbsjobname = '#PBS -N trim_%s\n' % (sampleid)
            pbscredit = '#PBS -A %s\n\n' % (creditaccount)
            pbsaddsource = 'source %s\n\n' % (addsourcetool)
            fo.write('#!/bin/bash -l\n\n')
            fo.write(pbsjobtime)
            fo.write(pbsnode)
            fo.write('#PBS -l pmem=5gb\n')
            fo.write(pbsuser)
            fo.write('#PBS -m abe\n')
            fo.write(pbsjobname)
            fo.write(pbscredit)
            fo.write(pbsaddsource)

            for fastqfile in os.listdir(fastqdir):
                if fastqfile.startswith(sampleid) and (fastqfile.endswith('R1.fastq.gz') or fastqfile.endswith('R1_001.fastq.gz')):
                    curreadR2 = fastqfile.replace("R1","R2")
                    readpairs = '%s/%s %s/%s' % (fastqdir, fastqfile, fastqdir, curreadR2)

                    if (DNAtrimoption == "predefined"):
                        addparam = '--clip_R1 10 --clip_R2 15 --three_prime_clip_R1 10 --three_prime_clip_R2 10'
                    else:
                        addparam = DNAtrimoption

                    trim_task = '%s %s -o %s --paired %s --fastqc\n\n' % (trimtool, addparam, trim_workdir, readpairs)
                    fo.write(trim_task)

    write_sub_script_qsub(trim_scriptdir, 'DNAtrim', trim_workdir, logdir)


def write_sub_script_DNAalign(DNASampleIDs, fastqdir, scriptdir, workdir, logdir):
    trim_workdir = workdir + '/2_DNAtrim'
    DNAalign_scriptdir = scriptdir + '/3_Bismarkalign'
    DNAalign_workdir = workdir + '/3_Bismark'
    helper.mkdir(DNAalign_scriptdir)
    helper.mkdir(DNAalign_workdir)

    for sampleid in DNASampleIDs:
        samplescript = DNAalign_scriptdir + '/' + sampleid + '.bismark.DNAalign.task.pbs'
        sample_sub_workdir = DNAalign_workdir + '/' + sampleid
        helper.mkdir(sample_sub_workdir)
        lanefastqR1 = []
        lanefastqR2 = []

        with open(samplescript, 'w') as fo:
            pbsjobtime = '#PBS -l walltime=%s:00:00\n' % (24)
            pbsnode = '#PBS -l nodes=1:ppn=%s\n' % (12)
            pbsuser = '#PBS -M %s\n' % (noticeaccount)
            pbsjobname = '#PBS -N DNAalign_%s\n' % (sampleid)
            pbscredit = '#PBS -A %s\n\n' % (creditaccount)
            pbsaddsource = 'source %s\n\n' % (addsourcetool)
            fo.write('#!/bin/bash -l\n\n')
            fo.write(pbsjobtime)
            fo.write(pbsnode)
            fo.write('#PBS -l pmem=5gb\n')
            fo.write(pbsuser)
            fo.write('#PBS -m abe\n')
            fo.write(pbsjobname)
            fo.write(pbscredit)
            fo.write(pbsaddsource)

            for fastqfile in os.listdir(trim_workdir):
                if fastqfile.startswith(sampleid) and (fastqfile.endswith('R1_val_1.fq.gz') or fastqfile.endswith('R1_001_val_1.fq.gz')):
                    curfastqR1 = '%s/%s' % (trim_workdir, fastqfile)
                    lanefastqR1.append(curfastqR1)

            trimFQR1 = ','.join(lanefastqR1)
            trimFQR2 = trimFQR1.replace("R1_","R2_").replace("val_1","val_2")

            #directional library
            DNAalign_task = '%s -N 1 -L 20 --multicore 6 --dovetail --gzip --path_to_bowtie %s --samtools_path %s --genome %s -1 %s -2 %s --output_dir %s --temp_dir %s\n\n' % (bismark, bowtiedir, samtoolsdir, bismarkreference, trimFQR1, trimFQR2, sample_sub_workdir, sample_sub_workdir)
            fo.write(DNAalign_task)

            alignbams = trimFQR1.replace(trim_workdir, sample_sub_workdir).replace(".fq.gz","_bismark_bt2_pe.bam").replace(",", " ")
            #check number of lanes
            if (',' in trimFQR1):
                dedup_task = '%s -p -o %s --output_dir %s --multiple %s\n\n' % (bismarkdedup, sampleid, sample_sub_workdir, alignbams)
                dedupbam = '%s/%s.multiple.deduplicated.bam' % (sample_sub_workdir, sampleid)
            else:
                dedup_task = '%s -p -o %s --output_dir %s %s\n\n' % (bismarkdedup, sampleid, sample_sub_workdir, alignbams)
                dedupbam = '%s/%s.deduplicated.bam' % (sample_sub_workdir, sampleid)
            fo.write(dedup_task)

            #filterbam = '%s/%s.merge.deduplicated.filtered.bam' % (sample_sub_workdir, sampleid)
            #filter_task = "%s view -bS -@ 4 -F 4 -F 256 -F 2048 -q 40 %s > %s\n\n" % (samtools, dedupbam, filterbam)
            #fo.write(filter_task)

            extract_task = '%s -p --gzip --parallel 4 --samtools_path %s -o %s --bedGraph %s\n\n' % (bismarkextractor, samtoolsdir, sample_sub_workdir, filterbam)
            fo.write(extract_task)

    write_sub_script_qsub(DNAalign_scriptdir, 'DNAalign', DNAalign_workdir, logdir)


def write_sub_script_coverage(DNASampleIDs, scriptdir, workdir, logdir):
    cov_scriptdir = scriptdir + '/4_cov'
    cov_workdir = workdir + '/4_cov'
    helper.mkdir(cov_scriptdir)
    helper.mkdir(cov_workdir)

    for sampleid in DNASampleIDs:

        samplescript = cov_scriptdir + '/' + sampleid + '.cov.task.pbs'
        samplescript2 = cov_scriptdir + '/' + sampleid + '.readcount.task.pbs'
        cov_sub_workdir = cov_workdir + '/' + sampleid
        helper.mkdir(cov_sub_workdir)

        with open(samplescript, 'w+') as fo:
            pbsjobtime = '#PBS -l walltime=%s:00:00\n' % (12)
            pbsnode = '#PBS -l nodes=1:ppn=%s\n' % (3)
            pbsuser = '#PBS -M %s\n' % (noticeaccount)
            pbsjobname = '#PBS -N cov_%s\n' % (sampleid)
            pbscredit = '#PBS -A %s\n\n' % (creditaccount)
            pbsaddsource = 'source %s\n\n' % (addsourcetool)
            fo.write('#!/bin/bash -l\n\n')
            fo.write(pbsjobtime)
            fo.write(pbsnode)
            fo.write('#PBS -l pmem=5gb\n')
            fo.write(pbsuser)
            fo.write('#PBS -m abe\n')
            fo.write(pbsjobname)
            fo.write(pbscredit)
            fo.write(pbsaddsource)

            #each step based on previous step; maybe wrap previous step into this step
            DNAalign_workdir = workdir + '/3_Bismark'
            DNAalign_sub_workdir = DNAalign_workdir + '/' + sampleid
            filteredbam = '%s/%s.merge.deduplicated.filtered.bam' % (DNAalign_sub_workdir, sampleid)
            ###bismark bam not in standard GATK bam format, need to add readgroup and sort before feeding into GATK programs
            infilteredbam = '%s/%s.merge.deduplicated.filtered.sorted.bam' % (DNAalign_sub_workdir, sampleid)
            addrgsort_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=BSseq RGPL=Illumina RGSM=%s RGPU=GC SORT_ORDER=coordinate CREATE_INDEX=true\n\n' % (java, logdir, picard, filteredbam, infilteredbam, sampleid)
            fo.write(addrgsort_task)
            indexbam_task = '%s index %s\n\n' % (samtools, infilteredbam)
            fo.write(indexbam_task)

            calfinsert_task = '%s %s/upstream/getinsertsize.py %s\n\n' % (python, runscriptdir, infilteredbam)
            fo.write(calfinsert_task)

            #count coverage in terms of fragments (overlapped pair-end reads)
            fragcov = '%s.merge.deduplicated.filtered.sorted.fragcov' % (sampleid)
            fragcov_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T DepthOfCoverage -R %s -o %s/%s -I %s -L %s --countType COUNT_FRAGMENTS -drf DuplicateRead\n\n' % (java, logdir, gatk, reference, cov_sub_workdir, fragcov, infilteredbam, target)
            fo.write(fragcov_task)
            ontarbam = '%s/%s.merge.deduplicated.filtered.sorted.ontar.bam' % (DNAalign_sub_workdir, sampleid)
            targetbam_task = '%s view -hbS %s -L %s > %s\n\n' % (samtools, infilteredbam, target, ontarbam)
            fo.write(targetbam_task)
            allinsert = '%s.merge.deduplicated.filtered.sorted.insert_size_metrics' % (sampleid)
            allinsertS_task = '%s -Xmx4G -jar %s CollectInsertSizeMetrics I=%s O=%s/%s.txt H=%s/%s.pdf DEVIATIONS=50 INCLUDE_DUPLICATES=false\n\n' % (java, picard, infilteredbam, cov_sub_workdir, allinsert, cov_sub_workdir, allinsert)
            fo.write(allinsertS_task)
            ontarinsert = '%s.merge.deduplicated.filtered.sorted.ontar.insert_size_metrics' % (sampleid)
            ontarinsertS_task = '%s -Xmx4G -jar %s CollectInsertSizeMetrics I=%s O=%s/%s.txt H=%s/%s.pdf DEVIATIONS=50 INCLUDE_DUPLICATES=false\n\n' % (java, picard, ontarbam, cov_sub_workdir, ontarinsert, cov_sub_workdir, ontarinsert)
            fo.write(ontarinsertS_task)

        ##output of CountReads is problematic, needs to be solved
        with open(samplescript2, 'w+') as fo2:
            pbsjobtime = '#PBS -l walltime=%s:00:00\n' % (12)
            pbsnode = '#PBS -l nodes=1:ppn=%s\n' % (3)
            pbsuser = '#PBS -M %s\n' % (noticeaccount)
            pbsjobname = '#PBS -N readcount_%s\n' % (sampleid)
            pbscredit = '#PBS -A %s\n\n' % (creditaccount)
            pbsaddsource = 'source %s\n\n' % (addsourcetool)
            fo2.write('#!/bin/bash -l\n\n')
            fo2.write(pbsjobtime)
            fo2.write(pbsnode)
            fo2.write('#PBS -l pmem=5gb\n')
            fo2.write(pbsuser)
            fo2.write('#PBS -m abe\n')
            fo2.write(pbsjobname)
            fo2.write(pbscredit)
            fo2.write(pbsaddsource)

            #count total filtered reads
            #CountReads write to stdout
            #readcountall = '%s.sorted.merged.md.filtered.readcountall.out' % (sampleid)
            readcountall_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T CountReads -R %s -I %s -drf DuplicateRead\n' % (java, logdir, gatk, reference, infilteredbam)
            fo2.write(readcountall_task)
            #count on target reads
            #readcountontar = '%s.sorted.merged.md.filtered.readcountOnTar.out' % (sampleid)
            readcountontar_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T CountReads -R %s -I %s -L %s -drf DuplicateRead\n' % (java, logdir, gatk, reference, infilteredbam, target)
            fo2.write(readcountontar_task)
            #count flanked on target reads
            #readcountontar50 = '%s.sorted.merged.md.filtered.readcountOnTar.50bpflank.out' % (sampleid)
            readcountontar50_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T CountReads -R %s -I %s -L %s -drf DuplicateRead\n' % (java, logdir, gatk, reference, infilteredbam, targetflank50)
            fo2.write(readcountontar50_task)

            #readcountontar100 = '%s.sorted.merged.md.filtered.readcountOnTar.100bpflank.out' % (sampleid)
            readcountontar100_task = '%s -Xmx4G -XX:-UsePerfData -XX:-UseParallelGC -Djava.io.tmpdir=%s -jar %s -T CountReads -R %s -I %s -L %s -drf DuplicateRead\n' % (java, logdir, gatk, reference, infilteredbam, targetflank100)
            fo2.write(readcountontar100_task)

    write_sub_script_qsub(cov_scriptdir, 'cov', cov_workdir, logdir)
    write_sub_script_qsub(cov_scriptdir, 'readcount', cov_workdir, cov_workdir)


def write_scripts():
    mytime = time.strftime('%y_%m_%d_%H_%M')
    scriptdir = workdir + '/scripts_' + mytime
    helper.mkdir(scriptdir)
    logdir = scriptdir + '/log'
    helper.mkdir(logdir)

    write_sub_script_fastqc(DNASampleIDs, fastqdir, scriptdir, workdir, logdir)
    write_sub_script_DNAtrimming(DNASampleIDs, fastqdir, scriptdir, workdir, logdir, DNAtrimoption)
    write_sub_script_DNAalign(DNASampleIDs, fastqdir, scriptdir, workdir, logdir)
    write_sub_script_coverage(DNASampleIDs, scriptdir, workdir, logdir)


if __name__ == '__main__':
    write_scripts()
