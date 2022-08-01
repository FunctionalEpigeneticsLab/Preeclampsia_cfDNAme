script.dir <- "./"
glmeval <- file.path(script.dir,"glmeval.R")
source(glmeval)

set.seed(2022)
args <- commandArgs(TRUE)
sampleinfo <- args[1]
inputdir <- args[2]
flagindexfh <- args[3]
normalization <- args[4]
outmat <- args[5]
lowvarfilter <- args[6]
alpha <- args[7]
mylambda <- args[8]
nfold <- args[9]
numrep <- args[10]
coefout <- args[11]
coefsumout <- args[12]
predresout <- args[13]
devmseout <- args[14]
perfoutfig <- args[15]
selected.feat <- args[16]

cntoption = "mval"

ftmat = FilterCountMatrixFeat(sampleinfo, inputdir, flagindexfh, cntoption, normalization, outmat, lowvarfilter)
RunGLMAssessReplicates(ftmat, flagindex, alpha, mylambda, nfold, numrep, coefout, coefsumout, predresout, devmseout, perfoutfig, selected.feat=NA)