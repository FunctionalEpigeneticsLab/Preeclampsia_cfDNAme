library(data.table)
library(reshape2)
library(caret)
library(e1071)
library(pROC)
library(glmnet)
library(limma)
# model training and prediction

LoadProbeIndex <- function(indexfh) {
    idx <- fread(indexfh, header=TRUE, sep="\t", data.table=FALSE)
    idx <- idx[order(idx$Index),]
    return(idx)
}

GetRawCountMatrix <- function(sampleinfo, inputdir, flagindexfh, cntoption="avgme") {
    #' @title Read multiple samples into a data matrix
    #' @param sampleinfo A tab-separated file contains structured information 'SubjectID\tPhenotype\tDescription\tGA\tSampleID\tMMeanDep\tMMedianDep\tBSFlag\tCenter\tSeqBatch\tYear\n'
    #' @param inputdir A directory path to methylation count of all samples
    #' @param flagindexfh A tab-separated file contains capture information 'Chromosome\tStart\tEnd\tIndex\tProbe\tFlagIndex'
    
    #saminfo <- fread(sampleinfo, header=TRUE, sep="\t", colClasses=c("character","character","character","numeric","character","numeric","numeric","character","character","character","numeric"), data.table=FALSE)
    saminfo <- fread(sampleinfo, header=TRUE, sep="\t", data.table=FALSE)
    flagidx <- LoadProbeIndex(flagindexfh)
    #datalist <- lapply()
    inmat <- matrix()
    
    for (i in 1:nrow(saminfo)) {
        print(paste("reading in subject ", saminfo[i,"SubjectID"]))
        samfile <- paste0(saminfo[i,"SubjectID"], ".ind.mavg.count.merge.tsv")
        samfh <- file.path(inputdir, samfile)
        curfh <- fread(samfh, header=TRUE, sep="\t", data.table=FALSE)
        curfh <- merge(flagidx,curfh,by=c("Chromosome","Start","End","Index","Probe"),all.x=TRUE)
        curfh <- curfh[order(curfh$Index),]
		
        if (cntoption == "avgme") {
            curfh$calval <- curfh$MePer/100
        } else if (cntoption == "logavgme") {
            curfh$calval <- log2(curfh$MePer/100+0.0001)
        } else {
            print("specify value to be calculated")
        }
            
        if (i == 1) {
            inmat <- matrix(curfh$calval,nrow=1)
        } else {
            inmat <- rbind(inmat, matrix(curfh$calval,nrow=1))
        }
    }
    #colmedian <- apply(inmat,2,median,na.rm=TRUE)
    #medianmat<-t(matrix(colmedian,ncol=nrow(inmat),nrow=ncol(inmat)))
    #inmat[is.na(inmat)] <- medianmat[is.na(inmat)]
    if (cntoption == "avgme") {
        inmat[is.na(inmat)] <- 0
    } else if (cntoption == "logavgme") {
        inmat[is.na(inmat)] <- log2(0.0001)
    }
    
    colnames(inmat) <- flagidx$Index
    rownames(inmat) <- paste0(saminfo$SubjectID, ":", saminfo$Phenotype, ":", saminfo$Batch)
    #write.table(inmat,outmat,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    return(inmat)
}

FlagFailedFilterFeat <- function(inmat, autosomeonly=FALSE) {
    flagidx <- LoadProbeIndex(flagindexfh)

    if (autosomeonly) {
        keepindex <- flagidx$Index[flagidx$FlagIndex==1 & flagidx$Chromosome != "X"]
    } else {
        keepindex <- flagidx$Index[flagidx$FlagIndex==1]
    }
    
    ftmat <- inmat[,colnames(inmat) %in% keepindex]
    return(ftmat)
}

NormalizeCountMatrix <- function(ftmat,cntoption,normalization=NA) {
    if (normalization=="meannorm") {
        if (cntoption == "logavgme") {
	    print("Normalize per sample methylation count by abs log MEAN")
	    #normfactor <- apply(ftmat, 1, function(x) abs(log2(mean(2^x))))
            normfactor <- apply(ftmat, 1, function(x) abs(mean(x)))
	} else if (cntoption == "avgme") {
	    print("Normalize per sample methylation count by MEAN")
            normfactor <- apply(ftmat, 1, function(x) abs(mean(x)))
	}
	print(normfactor)
        ftmatnorm <- ftmat/normfactor
    } else if (normalization=="mediannorm") {
        if (cntoption == "logavgme") {
            print("Normalize per sample methylation count by abs log MEDIAN")
            #normfactor <- apply(ftmat, 1, function(x) abs(log2(median(2^x)))
            normfactor <- apply(ftmat, 1, function(x) abs(median(x)))
        } else if (cntoption == "avgme") {
            print("Normalize per sample methylation count by MEDIAN")
            normfactor <- apply(ftmat, 1, function(x) abs(median(x)))
        }
	print(normfactor)
        ftmatnorm <- ftmat/normfactor
    } else if (normalization=="ffnorm") {
        #normalize on fetal fraction
        #consider to to
        print("to be implemented")
    } else {
	ftmatnorm <- ftmat
    }
    return(ftmatnorm)
}

FlagHighVarianceCtrl <- function(ftmat) {
    ctrlftmat <- ftmat[sapply(strsplit(row.names(ftmat), split=':', fixed=TRUE), function(x) (x[2]))=="Ctrl",]
    print(paste0("Control samples identified: ", dim(ctrlftmat)[1]))
    ctrlftvar <- apply(ctrlftmat, 2, var)
    lowvarindex <- colnames(ctrlftmat[,ctrlftvar < quantile(ctrlftvar, 0.75)])
    #elimftmat <- ftmat[,colnames(ftmat) %in% lowvarindex]
    #return(elimftmat)
    return(lowvarindex)
}

GetZCorrectionFromTrainCtrl <- function(batchftmat,batch) {
    ctrlftmat <- batchftmat[sapply(strsplit(row.names(batchftmat), split=':', fixed=TRUE), function(x) (x[2]))=="Ctrl",]
    print(paste0("Control samples identified to get reference mean and sd from batch ",batch,": ", dim(ctrlftmat)[1]))
    ctrlmean <- apply(ctrlftmat, 2, mean)
    ctrlsd <- apply(ctrlftmat, 2, sd)
    rescaledt <- data.frame(ctrlmean,ctrlsd)
    colnames(rescaledt) <- c("RefMean","RefSD")
    return(rescaledt)
}

ApplyZCorrection <- function(ftmat) {
    batch <- sapply(strsplit(row.names(ftmat), split=':', fixed=TRUE), function(x) (x[3]))
    zftmat <- data.frame()
    for (hb in unique(batch)) {
        batchftmat <- ftmat[sapply(strsplit(row.names(ftmat), split=':', fixed=TRUE), function(x) (x[3]))==hb,]
        refval <- GetZCorrectionFromTrainCtrl(batchftmat,hb)
        bftmat <- t((t(batchftmat)-refval$RefMean)/refval$RefSD)
        zftmat <- rbind(zftmat,bftmat)
    }
    return(zftmat)
}

GetCorrectionCoef <- function(ftmat, covariates=NULL) {
    ctrlftmat <- ftmat[sapply(strsplit(row.names(ftmat), split=':', fixed=TRUE), function(x) (x[2]))=="Ctrl",]
    batch <- sapply(strsplit(row.names(ctrlftmat), split=':', fixed=TRUE), function(x) (x[3]))
    print(paste0("Control samples identified to perform linear regression from batch ",unique(batch),": ", dim(ctrlftmat)[1]))
    design <- matrix(1,ncol(t(ctrlftmat)),1)

    if(is.null(batch) && is.null(covariates)) {
        return(as.matrix(t(ftmat)))
    }

    if(!is.null(batch)) {
	batch <- as.factor(batch)
	contrasts(batch) <- contr.sum(levels(batch))
	batch <- model.matrix(~batch)[,-1,drop=FALSE]
	}

    if(!is.null(covariates)) covariates <- as.matrix(covariates)
    X.batch <- cbind(batch,covariates)
    fit <- lmFit(t(ctrlftmat),cbind(design, X.batch))
    beta <- fit$coefficients[,-(1:ncol(design)),drop=FALSE]
    beta[is.na(beta)] <- 0
    return(beta)
}

ApplyBatchCorrection <- function(transftmat, beta, covariates=NULL, design=matrix(1,ncol(transftmat),1)) {
    batch <- sapply(strsplit(colnames(transftmat), split=':', fixed=TRUE), function(x) (x[3]))

    if(is.null(batch) && is.null(covariates)) {
        return(as.matrix(transftmat))
    }

    if(!is.null(batch)) {
	batch <- as.factor(batch)
	contrasts(batch) <- contr.sum(levels(batch))
	batch <- model.matrix(~batch)[,-1,drop=FALSE]
    }

    if(!is.null(covariates)) covariates <- as.matrix(covariates)
    X.batch <- cbind(batch,covariates)
    corftmat <- as.matrix(transftmat) - beta %*% t(X.batch)
    return(corftmat)
}


GetTrainMatrixFeat <- function(sampleinfo, inputdir, flagindexfh, cntoption, autosomeonly, normalization) {
    inmat <- GetRawCountMatrix(sampleinfo, inputdir, flagindexfh, cntoption)
    
    ftmat <- FlagFailedFilterFeat(inmat, autosomeonly)
    ftmatnorm <- NormalizeCountMatrix(ftmat, cntoption, normalization)
    print(paste0("Keep ", dim(ftmatnorm)[2], " features for ", dim(ftmatnorm)[1], " subjects for training assessment"))
    #write.table(ftmatnorm,outmat,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    return(ftmatnorm)
}

GetBuildMatrixFeat <- function(sampleinfo, inputdir, flagindexfh, cntoption, autosomeonly, batchcorrection, normalization=NA, lowvarfilter, outmat) {
    ftmatnorm <- GetTrainMatrixFeat(sampleinfo, inputdir, flagindexfh, cntoption, autosomeonly, normalization)

    if (batchcorrection == "znorm") {
        print("Applying z-normalization batch correction to the entire training set")
        ftmatnorm <- ApplyZCorrection(ftmatnorm)
        ftmatnorm <- as.matrix(ftmatnorm)
        ftmatnorm[!is.finite(ftmatnorm)] <- 0
    } else if (batchcorrection == "lmfit") {
        print("Applying linear regression batch correction to the entire training set")
        betacoef <- GetCorrectionCoef(ftmatnorm)
        ftmatnorm <- t(ApplyBatchCorrection(t(ftmatnorm),betacoef))
    } else {
        print("Applying no batch correction to the entire training set")
    }

    if (lowvarfilter) {
        print("Discard high varaince features identified in controls")
        lowvarindex <- FlagHighVarianceCtrl(ftmatnorm)
        ftmatnorm <- ftmatnorm[,colnames(ftmatnorm) %in% lowvarindex]
    }

    print(paste0("Keep ", dim(ftmatnorm)[2], " features for ", dim(ftmatnorm)[1], " subjects"))
    write.table(ftmatnorm,outmat,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    return(ftmatnorm)
}

GetPredMatrixFeat <- function(trainsampleinfo, predsampleinfo, inputdir, flagindexfh, cntoption, autosomeonly, batchcorrection, normalization, lowvarfilter, predoutmat) {
    trainftmatnorm <- GetTrainMatrixFeat(trainsampleinfo, inputdir, flagindexfh, cntoption, autosomeonly,normalization)

    predinmat <- GetRawCountMatrix(predsampleinfo, inputdir, flagindexfh, cntoption)
    predftmat <- FlagFailedFilterFeat(predinmat, autosomeonly)
    predftmatnorm <- NormalizeCountMatrix(predftmat, cntoption, normalization)

    if (batchcorrection == "znorm") {
        print("Applying z-normalization batch correction to the validation set")
        batch <- sapply(strsplit(row.names(trainftmatnorm), split=':', fixed=TRUE), function(x) (x[3]))
        zftmat <- data.frame()
        for (hb in unique(batch)) {
            batchtrain <- trainftmatnorm[sapply(strsplit(row.names(trainftmatnorm), split=':', fixed=TRUE), function(x) (x[3]))==hb,]
            refval <- GetZCorrectionFromTrainCtrl(batchtrain,hb)
            batchtest <- predftmatnorm[sapply(strsplit(row.names(predftmatnorm), split=':', fixed=TRUE), function(x) (x[3]))==hb,]
            bftmat <- t((t(batchtest)-refval$RefMean)/refval$RefSD)
            zftmat <- rbind(zftmat,bftmat)
        }
        predftmatnorm <- as.matrix(zftmat)
        predftmatnorm[!is.finite(predftmatnorm)] <- 0
        trainftmatnorm <- ApplyZCorrection(trainftmatnorm)
        trainftmatnorm <- as.matrix(trainftmatnorm)
        trainftmatnorm[!is.finite(trainftmatnorm)] <- 0
    } else if (batchcorrection == "lmfit") {
        print("Applying linear regression batch correction to the validation set")
        betacoef <- GetCorrectionCoef(trainftmatnorm)
        predftmatnorm <- t(ApplyBatchCorrection(t(predftmatnorm),betacoef))
        trainftmatnorm <- t(ApplyBatchCorrection(t(trainftmatnorm),betacoef))
    } else {
        print("Applying no batch correction to the validation set")
    }

    if (lowvarfilter) {
        lowvarindex <- FlagHighVarianceCtrl(trainftmatnorm)
        trainftmatnorm <- trainftmatnorm[,colnames(trainftmatnorm) %in% lowvarindex]
    }
    traindex <- colnames(trainftmatnorm)

    predftmatnorm <- predftmatnorm[,colnames(predftmatnorm) %in% traindex]
    print(paste0("*********Apply final model on ",dim(predftmatnorm)[1], " subjects*********"))
    print(paste0("Keep ", dim(predftmatnorm)[2]," features that were used in training input for prediction"))
    write.table(predftmatnorm,predoutmat,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    return(predftmatnorm)
}

FeatureGLMCV <- function(ftmatnorm, batchcorrection, normalization=NA, lowvarfilter, sidgroup, phenogroup, alpha, mylambda, nfold, curcycle, predresout, devmseout) {
    folds <- createFolds(phenogroup, k=nfold)
    predsmsel <- c()
    predsdevl <- c()
    alltruegroup <- phenogroup[unlist(folds)]
    allpredgroup <- c()
    predres <- data.frame()

    for (i in 1:nfold) {
        x <- folds[[i]]
        print(paste("Fold ",i))
        traindt <- ftmatnorm[-x,]
        traingroup <- phenogroup[-x]
        wts <- as.vector(1/(table(traingroup)[traingroup]/length(traingroup)))
        testdt <- ftmatnorm[x,]
        testgroup <- phenogroup[x]

        if (batchcorrection == "znorm") {
            batch <- sapply(strsplit(row.names(traindt), split=':', fixed=TRUE), function(x) (x[3]))
            zftmat <- data.frame()
            for (hb in unique(batch)) {
                batchtrain <- traindt[sapply(strsplit(row.names(traindt), split=':', fixed=TRUE), function(x) (x[3]))==hb,]
                refval <- GetZCorrectionFromTrainCtrl(batchtrain,hb)
                batchtest <- testdt[sapply(strsplit(row.names(testdt), split=':', fixed=TRUE), function(x) (x[3]))==hb,]
                bftmat <- t((t(batchtest)-refval$RefMean)/refval$RefSD)
                zftmat <- rbind(zftmat,bftmat)
            }
            testdt <- as.matrix(zftmat)
            testdt[!is.finite(testdt)] <- 0
            traindt <- ApplyZCorrection(traindt)
            traindt <- as.matrix(traindt)
            traindt[!is.finite(traindt)] <- 0
        } else if (batchcorrection == "lmfit") {
            betacoef <- GetCorrectionCoef(traindt)
            traindt <- t(ApplyBatchCorrection(t(traindt),betacoef))
            testdt <- t(ApplyBatchCorrection(t(testdt),betacoef))
        }

        if (lowvarfilter) {
            print("Discard high varaince features identified in controls")
            lowvarindex <- FlagHighVarianceCtrl(traindt)
            traindt <- traindt[,colnames(traindt) %in% lowvarindex]
            testdt <- testdt[,colnames(testdt) %in% lowvarindex]
        }
        
        cvmodel <- glmnet(traindt, traingroup, family="binomial", weights=wts, alpha=alpha, lambda=round(10^(seq(-3.5, 1.5, 0.05)),digit=5))
        predsdev <- assess.glmnet(cvmodel, newx=testdt, newy=testgroup, family="binomial", s=mylambda)$deviance
        predsdevl <- c(predsdevl, predsdev[1])
        predsmse <- assess.glmnet(cvmodel, newx=testdt, newy=testgroup, family="binomail", s=mylambda)$mse
        predsmsel <- c(predsmsel, predsmse[1])

        cv.preds.type <- predict(cvmodel, testdt, type="class", s=mylambda)
        cv.preds.res <- predict(cvmodel, testdt, type="response", s=mylambda)
        cv.preds.coef <- predict(cvmodel, testdt, type="coefficient", s=mylambda)

        coefdt <- data.frame(cv.preds.coef@Dimnames[[1]][cv.preds.coef@i + 1], cv.preds.coef@x)
        colnames(coefdt) <- c("feature", paste0("I",i,"R",curcycle))
        if (i == 1) {
                popcoef <- coefdt
                
            } else {
                popcoef <- merge(popcoef, coefdt, by=c("feature"),all=TRUE)
            }
        allpredgroup <- c(allpredgroup, cv.preds.res)
        predres <- cbind(sidgroup[x], cv.preds.type, cv.preds.res)
        write.table(predres,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
    }
    predsdevdt <- t(data.frame(c(paste0("DEV_R",curcycle),predsdevl)))
    write.table(predsdevdt,devmseout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    predsmsedt <- t(data.frame(c(paste0("MSE_R",curcycle),predsmsel)))
    write.table(predsmsedt,devmseout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)

    glmperfm <- roc(response=alltruegroup, predictor=allpredgroup, ci=TRUE,levels=c("Ctrl","Case"), direction="<")
    repauc <- paste0("AUC: ", glmperfm$auc)
    repci <- paste("CI: ", glmperfm$ci)
    youdenbest <- coords(glmperfm, "best", ret="threshold", transpose = FALSE, best.method="youden")
    topleftbest <- coords(glmperfm, "best", ret="threshold", transpose = FALSE, best.method="closest.topleft")
    print(paste0("Youden best cutoff: ", youdenbest))
    print(paste0("Topleft best cutoff: ", topleftbest))
    ssupper <- 1.0
    sslower <- 0.55
    spupper <- 1.0
    splower <- 0.7
    glmcoords <- coords(roc=glmperfm, x = "all", transpose = FALSE)
    coordacc <- glmcoords[(glmcoords$specificity >= splower & glmcoords$specificity <= spupper & glmcoords$sensitivity >= sslower & glmcoords$sensitivity <= ssupper),]
    print(coordacc)
    youdenaim <- coords(roc=glmperfm, "best", ret=c("threshold","specificity","sensitivity","accuracy","precision","recall"), transpose=FALSE, best.method="youden")
    topaim <- coords(roc=glmperfm, "best", ret=c("threshold","specificity","sensitivity","accuracy","precision","recall"), transpose=FALSE, best.method="closest.topleft")
    write.table(repauc,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
    write.table(repci,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
    write.table(c("Youden: ", youdenaim), predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
    write.table(c("Topleft: ", topaim), predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
    if (curcycle==1) {
        plot(glmperfm,print.auc=TRUE,print.auc.y=0.6-0.05*(curcycle-1),print.auc.x=0.3)
    } else {
        plot(glmperfm,add=TRUE,print.auc=TRUE,print.auc.y=0.6-0.05*(curcycle-1),print.auc.x=0.3)
    }
    return(popcoef)
}


RunGLMAssessReplicates <- function(ftmatnorm, batchcorrection, lowvarfilter, flagindex, alpha, mylambda, nfold, numrep, coefout, coefsumout, predresout, devmseout, perfoutfig, selected.feat=NA) {
    nfold <- as.numeric(nfold)
    numrep <- as.numeric(numrep)
    
    alpha <- as.numeric(alpha)
    mylambda <- as.numeric(mylambda)
    glmparam <- paste("#GLM-alpha: ",alpha,", lambda: ",mylambda,"; using selected features: ",selected.feat)
    print(glmparam)
    write.table(glmparam,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)

    sidgroup <- row.names(ftmatnorm)
    phenogroup <- sapply(strsplit(sidgroup, split=':', fixed=TRUE), function(x) (x[2]))
    phenogroup <- as.factor(phenogroup)
    if (length(levels(phenogroup)) == 2) {
        levels(phenogroup) <- list("Ctrl"=0, "Case"=1)
    } else {
        stop("number of group greater than 2")
    }
	
    if (!is.na(selected.feat)) {
        featidx <- fread(selected.feat, header=TRUE, sep="\t",data.table=FALSE)
        fidx <- featidx[order(featidx$Index),]
        sfeat <- fidx$Index
        ftmat <- ftmat[,colnames(ftmat) %in% sfeat]
        coefheader <- t(data.frame(sfeat))
        write.table(coefheader,coefout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    } else {
        flagidx <- fread(flagindexfh, header=TRUE, sep="\t", data.table=FALSE)
        fidx <- flagidx[flagidx$FlagIndex==1,]
        fidx <- fidx[order(fidx$Index),]
        coefheader <- t(data.frame(fidx$Index))
        write.table(coefheader,coefout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    }
	
    repcoef <- data.frame()
    pdf(perfoutfig,width=9,height=9)
    for (j in 1:numrep) {
        #set.seed(j)
        pfline <- paste("Replicate ",j)
        print(pfline)
        predresheader <- paste("##CV Rep",j)
        write.table(predresheader,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
        #write.table(predresheader,devmseout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
        popcoef <- FeatureGLMCV(ftmatnorm, batchcorrection, normalization, lowvarfilter, sidgroup, phenogroup, alpha, mylambda, nfold, curcycle=j, predresout, devmseout)

        if (j == 1) {
            repcoef <- popcoef
        } else {
            repcoef <- merge(repcoef, popcoef, by=c("feature"),all=TRUE)
        }
    }
    dev.off()
    write.table(repcoef,coefout,row.names=FALSE,quote=FALSE,sep="\t")
    freqthresh <- 0.75*(nfold*numrep)
    freqindex <- repcoef[rowSums(is.na(repcoef[,-1]))<=freqthresh,]
    freqindex$coefmean <- apply(freqindex[,-1], 1, mean, na.rm=TRUE)
    freqindex$coefsd <- apply(freqindex[,-1], 1, sd, na.rm=TRUE)
    outdf <- freqindex[,c("feature","coefmean","coefsd")]
    write.table(outdf,coefsumout,row.names=FALSE,quote=FALSE,sep="\t")
}


RunGLMcoefReplicates <- function(ftmat, alpha, nfold, numrep, coefout, coefsumout, mseaucout) {
    nfold <- as.numeric(nfold)
    numrep <- as.numeric(numrep)
    
    alpha <- as.numeric(alpha)
    glmparam <- paste("Pass 1 variable selection\nGLM-alpha: ",alpha)
    print(glmparam)

    sidgroup <- row.names(ftmat)
    phenogroup <- sapply(strsplit(sidgroup, split=':', fixed=TRUE), function(x) (x[2]))
    phenogroup <- as.factor(phenogroup)

    if (length(levels(phenogroup)) == 2) {
        levels(phenogroup) <- list("Ctrl"=0, "Case"=1)
    } else {
        stop("number of group greater than 2")
    }
    
    repcoef <- data.frame()
    for (j in 1:numrep) {
        #set.seed(j)
        pfline <- paste("Replicate ",j)
        print(pfline)
        popcoef <- GLMvariableCV(ftmat, phenogroup, alpha, nfold, curcycle=j, mseaucout)
        if (j == 1) {
            repcoef <- popcoef
        } else {
            repcoef <- merge(repcoef, popcoef, by=c("feature"),all=TRUE)
        }
    }
    write.table(repcoef,coefout,row.names=FALSE,quote=FALSE,sep="\t")
    freqthresh <- 0.5*(nfold*numrep)
    freqindex <- repcoef[rowSums(is.na(repcoef[,-1]))<=freqthresh,]
    freqindex$coefmean <- apply(freqindex[,-1], 1, mean, na.rm=TRUE)
    freqindex$coefsd <- apply(freqindex[,-1], 1, sd, na.rm=TRUE)
    outdf <- freqindex[,c("feature","coefmean","coefsd")]
    write.table(outdf,coefsumout,row.names=FALSE,quote=FALSE,sep="\t")
}

GenerateGLMmodel <- function(ftmatnorm, flagindex, alpha, selected.feat=NA, modeldir, modelname) {
    alpha <- as.numeric(alpha)
    glmparam <- paste("Generating GLM model - alpha: ",alpha,"; using specified features: ",selected.feat)
    print(glmparam)
    
    sidgroup <- row.names(ftmatnorm)
    phenogroup <- sapply(strsplit(sidgroup, split=':', fixed=TRUE), function(x) (x[2]))
    phenogroup <- as.factor(phenogroup)

    if (length(levels(phenogroup)) == 2) {
        levels(phenogroup) <- list("Ctrl"=0, "Case"=1)
    } else {
        stop("number of group greater than 2")
    }

    if (!is.na(selected.feat)) {
        featidx <- fread(selected.feat, header=TRUE, sep="\t",data.table=FALSE)
        fidx <- featidx[order(featidx$Index),]
        sfeat <- fidx$Index
        ftmat <- ftmat[,colnames(ftmat) %in% sfeat]
    } else {
        flagidx <- fread(flagindexfh, header=TRUE, sep="\t", data.table=FALSE)
        fidx <- flagidx[flagidx$FlagIndex==1,]
        fidx <- fidx[order(fidx$Index),]
    }
    
    wts <- as.vector(1/(table(phenogroup)[phenogroup]/length(phenogroup)))
    finalrgml <- glmnet(ftmatnorm, phenogroup, family="binomial", weights=wts, alpha=alpha, lambda=round(10^(seq(-3.5, 1.5, 0.05)),digit=5))
    print("writing GLM model")
    savedmodel <- paste0(modeldir, "/", modelname, ".rds")
    saveRDS(finalrgml, savedmodel)
}

ApplyGLMmodel <- function(predmatnorm, selected.feat=NA, mylambda, modeldir, modelname, outfh, outcoef, outfig) {
    mylambda <- as.numeric(mylambda)
    sidgroup <- row.names(predmatnorm)
    annoid <- sapply(strsplit(sidgroup, split=':', fixed=TRUE), function(x) (x[2]))
    annogroup <- sapply(strsplit(sidgroup, split=':', fixed=TRUE), function(x) (x[2]))

    if (!is.na(selected.feat)) {
        featidx <- fread(selected.feat, header=TRUE, sep="\t",data.table=FALSE)
        fidx <- featidx[order(featidx$Index),]
        sfeat <- fidx$Index
        predmatnorm <- predmatnorm[,colnames(predmatnorm) %in% sfeat]
        coefheader <- t(data.frame(sfeat))
        write.table(coefheader,outcoef,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    } else {
        flagidx <- fread(flagindexfh, header=TRUE, sep="\t", data.table=FALSE)
        fidx <- flagidx[flagidx$FlagIndex==1,]
        fidx <- fidx[order(fidx$Index),]
        coefheader <- t(data.frame(fidx$Index))
        write.table(coefheader,outcoef,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    }
    print(paste0("Classifying ",dim(predmatnorm)[1]," subjects using ", modelname, ".rds"))
    write.table(paste0("Using ", modeldir, "/", modelname, ".rds\n"),outfh,row.names=FALSE,col.names=FALSE,quote=FALSE)
    usemodel <- readRDS(paste0(modeldir, "/", modelname, ".rds"))
    nsam <- nrow(predmatnorm)
    predlist <- sapply(1:nsam, function(x) {
        predsample <- predmatnorm[x,,drop=FALSE]
        rglm.pred.type <- predict(usemodel, predsample, type="class", s=mylambda)
        rglm.pred.res <- predict(usemodel, predsample, type="response", s=mylambda)
        rglm.pred.coef <- predict(usemodel, predsample, type="coefficient", s=mylambda)
        coefdt <- data.frame(rglm.pred.coef@Dimnames[[1]][rglm.pred.coef@i + 1], rglm.pred.coef@x)
        colnames(coefdt) <- c("Index", paste0("coef.",sidgroup[x]))
        #merge to flagindex, discard intercept
        coefdt <- merge(fidx, coefdt, by=c("Index"), all.x=TRUE)
        coefdt <- coefdt[order(coefdt$Index),]
        loocoef <- matrix(coefdt[,c(paste0("coef.",sidgroup[x]))],nrow=1)
	write.table(loocoef,outcoef,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
        write.table(paste("predict", sidgroup[x], "as", rglm.pred.type, "-", rglm.pred.res),outfh,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
        return(rglm.pred.res)
    })
    glmperfm <- roc(response=annogroup, predictor=predlist, ci=TRUE,levels=c("Ctrl","Case"), direction="<")
    pdf(outfig, width=8, height=8)
    plot(glmperfm,col="#b2182b",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.5,ci=TRUE,lwd=4)
    legend("bottomright", legend=c("Validation"),col=c("#b2182b"), lwd=2)
    dev.off()
    print(glmperfm$auc)
    print(glmperfm$ci)
    youdenbest <- coords(glmperfm, "best", ret="threshold", transpose = FALSE, best.method="youden")
    topleftbest <- coords(glmperfm, "best", ret="threshold", transpose = FALSE, best.method="closest.topleft")
    print(paste0("Youden best cutoff: ", youdenbest))
    print(paste0("Topleft best cutoff: ", topleftbest))
    ssupper <- 1.0
    sslower <- 0.55
    spupper <- 1.0
    splower <- 0.65
    glmcoords <- coords(roc=glmperfm, x = "all", transpose = FALSE)
    coordacc <- glmcoords[(glmcoords$specificity >= splower & glmcoords$specificity <= spupper & glmcoords$sensitivity >= sslower & glmcoords$sensitivity <= ssupper),]
    print(coordacc)
}
