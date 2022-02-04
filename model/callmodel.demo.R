script.dir <- dirname(sys.frame(1)$ofile)
glmeval <- file.path(script.dir,"glmeval.R")
source(glmeval)

#set.seed(2022)
args <- commandArgs(TRUE)
sampleinfo <- args[1]
inputdir <- args[2]
indexfh <- args[3]
outmat <- args[4]
flagindexfh <- args[5]
alpha <- args[6]
lambda <- args[7]
outpred <- args[8]
outfig <- args[9]

cntoption <- "mval"

ftmat <- FilterCountMatrixFeat(sampleinfo, inputdir, indexfh, cntoption, outmat, flagindexfh)
AssessGLMLoo(ftmat, alpha, lambda, outpred, outfig)
