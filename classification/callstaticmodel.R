script.dir <- "/staging/leuven/stg_00064/Huiwen/project/3_cfDNA/11_NMrevision/7_behindthemoon/script"
glmeval <- file.path(script.dir,"glmeval.R")
source(glmeval)

args <- commandArgs(TRUE)
sampleinfotrain <- args[1]
sampleinfovalid <- args[2]
inputdir <- args[3]
flagindexfh <- args[4]
cntoption <- args[5]
batchcorrection <- args[6]
normalization <- args[7]
lowvarfilter <- args[8]
outmattrain <- args[9]
outmatvalid <- args[10]
alpha <- args[11]
selected.feat <- args[12]
modeldir <- args[13]
modelname <- args[14]
mylambda <- args[15]
outfh <- args[16]
outcoef <- args[17]
outfig <- args[18]

#cntoption = "avgme"
autosomeonly = TRUE

buildmat = GetBuildMatrixFeat(sampleinfotrain, inputdir, flagindexfh, cntoption, autosomeonly, batchcorrection, normalization, lowvarfilter, outmattrain)
predmat = GetPredMatrixFeat(sampleinfotrain, sampleinfovalid, inputdir, flagindexfh, cntoption, autosomeonly, batchcorrection, normalization, lowvarfilter, outmatvalid)
GenerateGLMmodel(buildmat, flagindex, alpha, selected.feat=NA, modeldir, modelname)
ApplyGLMmodel(predmat, selected.feat=NA, mylambda, modeldir, modelname, outfh, outcoef, outfig)
