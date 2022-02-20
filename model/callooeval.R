source(glmeval)

#set.seed(2022)
args <- commandArgs(TRUE)
sampleinfo <- args[1]
inputdir <- args[2]
flagindexfh <- args[3]
normalization <- args[4]
outmat <- args[5]
lowvarfilter <- args[6]
alpha <- args[7]
lambda <- args[8]
outpred <- args[9]
outcoef <- args[10]
outfig <- args[11]

cntoption = "mval"

ftmat = FilterCountMatrixFeat(sampleinfo, inputdir, flagindexfh, cntoption, normalization, outmat, lowvarfilter)
FeatureGLMLoo(ftmat, flagindexfh, alpha, lambda, outpred, outcoef, outfig, selected.feat=NA)
