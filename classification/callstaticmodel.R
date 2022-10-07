script.dir <- "./"
glmeval <- file.path(script.dir,"glmeval.R")
source(glmeval)

args <- commandArgs(TRUE)
sampleinfotrain <- args[1]
sampleinfovalid <- args[2]
inputdir <- args[3]
flagindexfh <- args[4]
normalization <- args[5]
lowvarfilter <- args[6]
outmattrain <- args[7]
outmatvalid <- args[8]
alpha <- args[9]
selected.feat <- args[10]
modeldir <- args[11]
modelname <- args[12]
mylambda <- args[13]
outfh <- args[14]
outcoef <- args[15]
outfig <- args[16]

cntoption = "mval"
autosomeonly = TRUE

buildmat = GetBuildMatrixFeat(sampleinfo, inputdir, flagindexfh, cntoption, autosomeonly, normalization=NA, lowvarfilter, outmat)
predmat = GetPredMatrixFeat(buildmat, predsampleinfo, inputdir, flagindexfh, cntoption, autosomeonly, normalization, lowvarfilter, predoutmat)
GenerateGLMmodel(ftmat, flagindex, alpha, selected.feat=NA, modeldir, modelname)
ApplyGLMmodel(predmat, selected.feat=NA, mylambda, modeldir, modelname, outfh, outcoef, outfig)
