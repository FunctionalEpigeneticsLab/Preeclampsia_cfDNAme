script.dir <- "/staging/leuven/stg_00064/Huiwen/project/3_cfDNA/11_NMrevision/7_behindthemoon/script"
glmeval <- file.path(script.dir,"glmeval.R")
source(glmeval)

set.seed(2022)
args <- commandArgs(TRUE)
sampleinfo <- args[1]
inputdir <- args[2]
flagindexfh <- args[3]
cntoption <- args[4]
normalization <- args[5]
batchcorrection <- args[6] #znorm/lmfit/nocorrection
lowvarfilter <- args[7]
alpha <- args[8]
mylambda <- args[9]
nfold <- args[10]
numrep <- args[11]
coefout <- args[12]
coefsumout <- args[13]
predresout <- args[14]
devmseout <- args[15]
perfoutfig <- args[16]

autosomeonly = TRUE

ftmatnorm = GetTrainMatrixFeat(sampleinfo, inputdir, flagindexfh, cntoption, autosomeonly, normalization)
RunGLMAssessReplicates(ftmatnorm, batchcorrection, lowvarfilter, flagindex, alpha, mylambda, nfold, numrep, coefout, coefsumout, predresout, devmseout, perfoutfig, selected.feat=NA)
