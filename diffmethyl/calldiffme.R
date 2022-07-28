script.dir <- "./"
limmatest <- file.path(script.dir,"limmatest.R")
source(limmatest)

args <- commandArgs(TRUE)
sampleinfo <- args[1]
inputdir <- args[2]
flagindexfh <- args[3] ## ./data/probe.index.n1flag.calibrate.tsv
cntoption <- args[4] ## betaval || mval
normalization <- args[5] ## meannorm || mediannorm || unnorm
material <- args[6] ## Freshplacenta || Blood || cfDNAatDiagnosis || OxcfDNAatDiagnosis || cfDNAfirstT || cfDNAfirstT_valid
outprefix <- args[7]

GetLimmaMatrix(sampleinfo, inputdir, flagindexfh, cntoption, normalization, material, outprefix)