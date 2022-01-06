#!/usr/bin/env python

import sys
import os
import configparser

curdir = os.getcwd()

def mkdir(sdir):
    if not os.path.isdir(sdir):
        os.makedirs(sdir)

def check_fexist(f):
    if os.path.isfile(f) == True:
        return(1)
    else:
        return(0)

def check_dirslash(d):
    if d.strip()[-1] == '/':
        return(d)
    else:
        d = d + '/'
        return(d)
    print(d)

def file_remove(fh):
    try:
        os.remove(fh)
    except OSError:
        raise

def get_dir_config():
    parser = configparser.ConfigParser()
    ConfigFile = curdir + '/ConfigFile'
    parser.read('ConfigFile')

    fastqdir = parser.get('configDir', 'fastqdir')
    workdir = parser.get('configDir', 'workdir')
    runscriptdir = parser.get('configDir', 'runscriptdir')

    return(fastqdir, workdir, runscriptdir)

def get_run_config():
    parser = configparser.ConfigParser()
    ConfigFile = curdir + '/ConfigFile'
    parser.read('ConfigFile')

    samtools = parser.get('configTool', 'samtools')
    samtoolsdir = parser.get('configTool', 'samtoolsdir')
    fastqc = parser.get('configTool', 'fastqc')
    trimtool = parser.get('configTool', 'trimtool')
    bismark = parser.get('configTool', 'bismark')
    bismarkdedup = parser.get('configTool', 'bismarkdedup')
    bismarkextractor = parser.get('configTool', 'bismarkextractor')
    bismarkcytosine = parser.get('configTool', 'bismarkcytosine')
    java = parser.get('configTool', 'java')
    gatk = parser.get('configTool', 'gatk')
    picard = parser.get('configTool', 'picard')
    addsourcetool = parser.get('configTool', 'addsourcetool')
    python = parser.get('configTool', 'python')

    bismarkreference = parser.get('configFile', 'bismarkreference')
    reference = parser.get('configFile', 'reference')
    target = parser.get('configFile', 'target')
    targetflank50 = parser.get('configFile', 'targetflank50')
    targetflank100 = parser.get('configFile', 'targetflank100')

    return(samtools, samtoolsdir, fastqc, trimtool, bismark, bismarkdedup, bismarkextractor, bismarkcytosine, java, gatk, picard, addsourcetool, python, bismarkreference, reference, target, targetflank50, targetflank100)

def get_DNA_sample():
    parser = configparser.ConfigParser()
    ConfigFile = curdir + '/ConfigFile'
    parser.read('ConfigFile')

    DNAsamplelist = parser.get('configFile', 'samplelist')

    DNAsampleID = []

    with open(DNAsamplelist, 'r') as fh0:
        for l0 in fh0:
            DNAsampleID.append(l0.rstrip('\n'))
    fh0.close()

    return(DNAsampleID)

def get_DNAproc_config():
    parser = configparser.ConfigParser()
    ConfigFile = curdir + '/ConfigFile'
    parser.read('ConfigFile')

    runscriptdir = parser.get('configDir', 'runscriptdir')
    Rscript = parser.get('configTool', 'Rscript')
    return(runscriptdir, Rscript)

def get_option():
    parser = configparser.RawConfigParser()
    ConfigFile = curdir + '/ConfigFile'
    parser.read('ConfigFile')

    DNAtrimoption = parser.get('configParam','DNAtrimoption')

    return(DNAtrimoption)

def get_userinfo():
    parser = configparser.RawConfigParser()
    ConfigFile = curdir + '/ConfigFile'
    parser.read('ConfigFile')

    noticeaccount = parser.get('configUser','noticeaccount')
    creditaccount = parser.get('configUser','creditaccount')

    return(noticeaccount, creditaccount)
