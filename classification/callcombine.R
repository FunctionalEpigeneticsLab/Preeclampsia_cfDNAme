script.dir <- "./"
combinescore <- file.path(script.dir,"combinescore.R")
source(combinescore)

args <- commandArgs(TRUE)
infh <- args[1]
FMFscore <- args[2] #FMFSCORE34 or FMFSCORE37
##output arguments
trainrocfig <- args[3]
trainoutfh <- args[4]
traincombrocfig <- args[5]
validoutfh <- args[6]
validrocfig <- args[7]
valid1rocfig <- args[8]
valid2rocfig <- args[9]
classoutfig <- args[10]

TrainingCombineFMF(infh,FMFscore,trainrocfig, trainoutfh, traincombrocfig)
ValidCombineFMF(infh, FMFscore, validoutfh, validrocfig, valid1rocfig, valid2rocfig, classoutfig)
