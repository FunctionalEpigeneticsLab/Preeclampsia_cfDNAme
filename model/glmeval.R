library(data.table)
library(reshape2)
library(caret)
library(e1071)
library(pROC)
library(glmnet)

LoadProbeIndex <- function(indexfh) {
    idx <- fread(indexfh, header=TRUE, sep="\t", data.table=FALSE)
    idx <- idx[order(idx$Index),]
    return(idx)
}

GetRawCountMatrix <- function(sampleinfo, inputdir, flagindexfh, cntoption="mval") {
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
	samfile <- paste0(saminfo[i,"SubjectID"], ".mavg.count.merge.tsv")
	samfh <- file.path(inputdir, samfile)
	curfh <- fread(samfh, header=TRUE, sep="\t", data.table=FALSE)
	curfh <- merge(flagidx,curfh,by=c("Chromosome","Start","End","Index","Probe"),all.x=TRUE)
	curfh <- curfh[order(curfh$Index),]
		
	if (cntoption == "mval") {
	    curfh$calval <- log2((curfh$Methylated+1)/(curfh$Unmethylated+1))
	} else if (cntoption == "betaval") {
	    curfh$calval <- curfh$Methylated/(curfh$Methylated+curfh$Unmethylated)
	} else {
	    print("specify value to be calculated")
	}
		
	if (i == 1) {
	    inmat <- matrix(curfh$calval,nrow=1)
	} else {
	    inmat <- rbind(inmat, matrix(curfh$calval,nrow=1))
	}
    }
    colnames(inmat) <- flagidx$Index
    rownames(inmat) <- paste0(saminfo$SubjectID,":",saminfo$Phenotype)
    #write.table(inmat,outmat,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    return(inmat)
}

FlagFailedFeat <- function(inmat) {
    flagidx <- LoadProbeIndex(flagindexfh)
    keepindex <- flagidx$Index[flagidx$FlagIndex==1]
    ftmat <- inmat[,colnames(inmat) %in% keepindex]
    return(ftmat)
}

NormalizeCountMatrix <- function(ftmat, normalization=NA) {
    if (normalization=="meannorm") {
	print("Normalize per sample methylation count by absMEAN")
	normfactor <- apply(ftmat, 1, function(x) abs(mean(x)))
	ftmatnorm <- ftmat/normfactor
    } else if (normalization=="mediannorm") {
	print("Normalize per sample methylation count by absMEDIAN")
	normfactor <- apply(ftmat, 1, function(x) abs(median(x)))
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

GetTrainMatrixFeat <- function(sampleinfo, inputdir, flagindexfh, cntoption, normalization=NA, outmat) {
    flagidx <- LoadProbeIndex(flagindexfh)
    inmat <- GetRawCountMatrix(sampleinfo, inputdir, flagindexfh, cntoption)
    
    ftmat <- FlagFailedFeat(inmat)
    ftmatnorm <- NormalizeCountMatrix(ftmat, normalization)
    print(paste0("Keep ", dim(ftmatnorm)[2], " features for ", dim(ftmatnorm)[1], " subjects for training assessment"))
    write.table(ftmatnorm,outmat,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    return(ftmatnorm)
}

GetBuildMatrixFeat <- function(sampleinfo, inputdir, flagindexfh, cntoption, normalization=NA, lowvarfilter=FALSE, outmat) {
    flagidx <- LoadProbeIndex(flagindexfh)
    inmat <- GetRawCountMatrix(sampleinfo, inputdir, flagindexfh, cntoption)
    
    ftmat <- FlagFailedFeat(inmat)
    ftmatnorm <- NormalizeCountMatrix(ftmat, normalization)
	
    if (lowvarfilter) {
        print("Discard high varaince features identified in controls")
        lowvarindex <- FlagHighVarianceCtrl(ftmatnorm)
        ftmatnorm <- ftmatnorm[,colnames(ftmatnorm) %in% lowvarindex]
    } 
    print(paste0("Keep ", dim(ftmatnorm)[2], " features for ", dim(ftmatnorm)[1], " subjects"))
    write.table(ftmatnorm,outmat,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    return(ftmatnorm)
}

GetPredMatrixFeat <- function(trainftmat, predsampleinfo, inputdir, flagindexfh, cntoption, normalization, lowvarfilter, predoutmat) {
    flagidx <- LoadProbeIndex(flagindexfh)
    traindex <- colnames(trainftmat)
    predinmat <- GetRawCountMatrix(predsampleinfo, inputdir, flagindexfh, cntoption)
    predftmat <- FlagFailedFeat(predinmat)
    predftmatnorm <- NormalizeCountMatrix(predftmat, normalization)
    predftmatnorm <- predftmatnorm[,colnames(predftmatnorm) %in% traindex]
    print(paste0("*********Apply final model on ",dim(predftmatnorm)[1], " subjects*********"))
    print(paste0("Keep ", dim(predftmatnorm)[2]," features that were used in training input for prediction"))
    write.table(predftmatnorm,predoutmat,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    return(predftmatnorm)
}

obsolete_GetCountMatrix <- function(sampleinfo, inputdir, flagindexfh, cntoption="mval", normalization=NA, outmat) {
    #' @title Read multiple samples into a data matrix
    #' @param sampleinfo A tab-separated file contains structured information 'SubjectID\tPhenotype\tDescription\tGA\tSampleID\tMMeanDep\tMMedianDep\tBSFlag\tCenter\tSeqBatch\tYear\n'
    #' @param inputdir A directory path to methylation count of all samples
    #' @param flagindexfh A tab-separated file contains capture information 'Chromosome\tStart\tEnd\tIndex\tProbe\tFlagIndex'
    #' @param outmat An output file for combined matrix
    
    saminfo <- fread(sampleinfo, header=TRUE, sep="\t", colClasses=c("character","character","character","numeric","character","numeric","numeric","character","character","character","numeric"), data.table=FALSE)
    flagidx <- LoadProbeIndex(flagindexfh)
    #datalist <- lapply()
    inmat <- matrix()
    
    for (i in 1:nrow(saminfo)) {
        print(paste("reading in subject ", saminfo[i,"SubjectID"]))
	samfile <- paste0(saminfo[i,"SubjectID"], ".mavg.count.merge.tsv")
	samfh <- file.path(inputdir, samfile)
	curfh <- fread(samfh, header=TRUE, sep="\t", data.table=FALSE)
	curfh <- merge(flagidx,curfh,by=c("Chromosome","Start","End","Index","Probe"),all.x=TRUE)
	curfh <- curfh[order(curfh$Index),]
		
	if (cntoption == "mval") {
	    curfh$calval <- log2((curfh$Methylated+1)/(curfh$Unmethylated+1))
	} else if (cntoption == "betaval") {
	    curfh$calval <- curfh$Methylated/(curfh$Methylated+curfh$Unmethylated)
	} else {
	    print("specify value to be calculated")
	}
		
	if (normalization=="meannorm") {
	    normfactor <- abs(mean(curfh$calval))
            print(paste0("Normalize the sample methylation count by absMEAN ", normfactor))
	    curfh$calval <- curfh$calval/normfactor
	} else if (normalization=="mediannorm") {
	    normfactor <- abs(median(curfh$calval))
	    print(paste0("Normalize the sample methylation count by absMEDIAN ", normfactor))
	    curfh$calval <- curfh$calval/normfactor
	} else if (normalization=="ffnorm") {
	    #normalize on fetal fraction
	    #consider to to
	    print("to be implemented")
	}
		
	if (i == 1) {
	    inmat <- matrix(curfh$calval,nrow=1)
	} else {
	    inmat <- rbind(inmat, matrix(curfh$calval,nrow=1))
	}
    }
    colnames(inmat) <- flagidx$Index
    rownames(inmat) <- paste0(saminfo$SubjectID,":",saminfo$Phenotype)
    write.table(inmat,outmat,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    return(inmat)
}

obsolete_AssessGLMLoo <- function(ftmat, flagindexfh, alpha, lambda, outpred, outcoef, outfig) {
    alpha <- as.numeric(alpha)
    lambda <- as.numeric(lambda)
    glmparam <- paste("GLM-alpha: ",alpha,", lambda: ",lambda)
    print(glmparam)
    flagidx <- fread(flagindexfh, header=TRUE, sep="\t", data.table=FALSE)
    fidx <- flagidx[flagidx$FlagIndex==1,]
    fidx <- fidx[order(fidx$Index),]
    coefheader <- t(data.frame(fidx$Index))
    write.table(coefheader,outcoef,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    
    sidgroup <- row.names(ftmat)
    phenogroup <- sapply(strsplit(sidgroup, split=':', fixed=TRUE), function(x) (x[2]))
    phenogroup <- as.factor(phenogroup)

    if (length(levels(phenogroup)) == 2) {
        levels(phenogroup) <- list("Ctrl"=0, "Case"=1)
    } else {
        stop("number of group greater than 2")
    }
    
    n_train <- nrow(ftmat)
    outheader <- paste(glmparam, "\nClassified samples:\n")
    write.table(outheader,outpred,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
    alltruegroup <- c()
    allpredgroup <- c()
    allpredval <- c()

    loopreds <- sapply(1:n_train, function(x) {
        traindt <- ftmat[-x,]
        traingroup <- phenogroup[-x]
        wts <- as.vector(1/(table(traingroup)[traingroup]/length(traingroup)))
        testdt <- t(as.data.frame(ftmat[x,]))
        testgroup <- phenogroup[x]
        rglm.model <- glmnet(traindt, traingroup, family="binomial", alpha=alpha, weights=wts, lambda=round(10^(seq(-3.5, 1.5, 0.05)),digit=5))
        rglm.preds.type <- predict(rglm.model, testdt, type="class", s=lambda)
        rglm.preds.res <- predict(rglm.model, testdt, type="response", s=lambda)
        rglm.preds.coef <- predict(rglm.model, testdt, type="coefficient", s=lambda)
        alltruegroup <- c(alltruegroup, testgroup)
        allpredval <- c(allpredval, rglm.preds.res)
        write.table(paste0("predict ", sidgroup[x], "\t", testgroup, " as ", rglm.preds.type, "\t", rglm.preds.res),outpred,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
	coefdt <- data.frame(rglm.preds.coef@Dimnames[[1]][rglm.preds.coef@i + 1], rglm.preds.coef@x)
	colnames(coefdt) <- c("Index", paste0("coef.",sidgroup[x]))
	coefdt <- merge(fidx, coefdt, by=c("Index"), all.x=TRUE)
	coefdt <- coefdt[order(coefdt$Index),]
        loocoef <- matrix(coefdt[,c(paste0("coef.",sidgroup[x]))],nrow=1)
	write.table(loocoef,outcoef,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
        return(allpredval)
    })
    
    glmperfm <- roc(response=phenogroup, predictor=loopreds, ci=TRUE, levels=c("Ctrl","Case"), direction="<")
    pdf(outfig, width=8, height=8)
    plot(glmperfm,print.auc=TRUE,print.auc.y=0.6,print.auc.x=0.3)
    dev.off()
    print(glmperfm$auc)
    print(glmperfm$ci)
    runacc <- table(phenogroup, ifelse(loopreds>0.5, "Case", "Ctrl"))
    print(runacc)
}

CoefRegular <- function(coeffh, flagindexfh, freqthresh, outfeat) {
    freqthresh <- as.numeric(freqthresh)
    cfdt <- fread(coeffh, header=TRUE, sep="\t", data.table=FALSE)
    cfdf <- data.frame(t(cfdt))
    cfdf$Index <- row.names(cfdf)
    flagidx <- fread(flagindexfh, header=TRUE, sep="\t", data.table=FALSE)
    dt <- merge(flagidx, cfdf, by=c("Index"),all.x=TRUE)

    featthresh <- (1-freqthresh)*dim(dt[,!(colnames(dt) %in% c("Index","Chromosome","Start","End","Probe","FlagIndex"))])[2]
    featindex <- dt[rowSums(is.na(dt[,!(colnames(dt) %in% c("Index","Chromosome","Start","End","Probe","FlagIndex"))]))<=featthresh,]
    featindex$coefmean <- apply(featindex[,!(colnames(featindex) %in% c("Index","Chromosome","Start","End","Probe","FlagIndex"))], 1, mean, na.rm=TRUE)
    featindex$coefsd <- apply(featindex[,!(colnames(featindex) %in% c("Index","Chromosome","Start","End","Probe","FlagIndex"))], 1, sd, na.rm=TRUE)
    outdf <- featindex[,c("Index","Chromosome","Start","End","Probe","FlagIndex","coefmean","coefsd")]
    write.table(outdf, outfeat, row.names=FALSE, quote=FALSE, sep="\t")
}

FeatureGLMLoo <- function(ftmat, lowvarfilter, flagindexfh, alpha, lambda, outpred, outcoef, outfig, selected.feat=NA) {
    alpha <- as.numeric(alpha)
    lambda <- as.numeric(lambda)
    glmparam <- paste("GLM-alpha: ",alpha,", lambda: ",lambda,"; using selected features: ",selected.feat)
    print(glmparam)
    
    sidgroup <- row.names(ftmat)
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
	write.table(coefheader,outcoef,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    } else {
        flagidx <- fread(flagindexfh, header=TRUE, sep="\t", data.table=FALSE)
	fidx <- flagidx[flagidx$FlagIndex==1,]
	fidx <- fidx[order(fidx$Index),]
	coefheader <- t(data.frame(fidx$Index))
	write.table(coefheader,outcoef,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    }
    
    n_train <- nrow(ftmat)
    outheader <- paste(glmparam, "\nClassified samples:\n")
    write.table(outheader,outpred,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
    alltruegroup <- c()
    allpredgroup <- c()
    allpredval <- c()

    loopreds <- sapply(1:n_train, function(x) {
        traindt <- ftmat[-x,]
	if (lowvarfilter) {
	    print("Discard high varaince features identified in controls")
            lowvarindex <- FlagHighVarianceCtrl(traindt)
	    traindt <- traindt[,colnames(traindt) %in% lowvarindex]
	    testdt <- t(as.data.frame(ftmat[x,colnames(ftmat) %in% lowvarindex]))
	} else {
	    testdt <- t(as.data.frame(ftmat[x,]))
	}
        traingroup <- phenogroup[-x]
        wts <- as.vector(1/(table(traingroup)[traingroup]/length(traingroup)))
	testgroup <- phenogroup[x]
	
        rglm.model <- glmnet(traindt, traingroup, family="binomial", alpha=alpha, weights=wts, lambda=round(10^(seq(-3.5, 1.5, 0.05)),digit=5))
        rglm.preds.type <- predict(rglm.model, testdt, type="class", s=lambda)
        rglm.preds.res <- predict(rglm.model, testdt, type="response", s=lambda)
        rglm.preds.coef <- predict(rglm.model, testdt, type="coefficient", s=lambda)
        alltruegroup <- c(alltruegroup, testgroup)
        allpredval <- c(allpredval, rglm.preds.res)
        write.table(paste0("predict ", sidgroup[x], "\t", testgroup, " as ", rglm.preds.type, "\t", rglm.preds.res),outpred,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
	coefdt <- data.frame(rglm.preds.coef@Dimnames[[1]][rglm.preds.coef@i + 1], rglm.preds.coef@x)
	colnames(coefdt) <- c("Index", paste0("coef.",sidgroup[x]))
	#merge to flagindex, discard intercept
	coefdt <- merge(fidx, coefdt, by=c("Index"), all.x=TRUE)
	coefdt <- coefdt[order(coefdt$Index),]
        loocoef <- matrix(coefdt[,c(paste0("coef.",sidgroup[x]))],nrow=1)
	write.table(loocoef,outcoef,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
        return(allpredval)
    })
    
    glmperfm <- roc(response=phenogroup, predictor=loopreds, ci=TRUE,levels=c("Ctrl","Case"), direction="<")
    pdf(outfig, width=8, height=8)
    plot(glmperfm,print.auc=TRUE,print.auc.y=0.6,print.auc.x=0.3)
    dev.off()
    print(glmperfm$auc)
    print(glmperfm$ci)
    youdenbest <- coords(glmperfm, "best", ret=c("threshold","specificity","sensitivity","accuracy","precision","recall"), transpose = FALSE, best.method="youden")
    topleftbest <- coords(glmperfm, "best", ret=c("threshold","specificity","sensitivity","accuracy","precision","recall"), transpose = FALSE, best.method="closest.topleft")

    print("Youden best cutoff: ")
    print(youdenbest)
    print("Topleft best cutoff: ")
    print(topleftbest)

    ssupper <- 1.0
    sslower <- 0.55
    spupper <- 1.0
    splower <- 0.7
    glmcoords <- coords(roc=glmperfm, x = "all", transpose = FALSE)
    coordacc <- glmcoords[(glmcoords$specificity >= splower & glmcoords$specificity <= spupper & glmcoords$sensitivity >= sslower & glmcoords$sensitivity <= ssupper),]
    print(coordacc)
    runacc <- table(phenogroup, ifelse(loopreds>0.5, "Case", "Ctrl"))
    print(runacc)
}


obsolete_GLMvariableCV <- function(ftmat, phenogroup, alpha, nfold, curcycle, mseaucout) {
    folds <- createFolds(phenogroup, k=nfold)
    predsmsel <- c()
    predsaucl <- c()
    predsdevl <- c()

    for (i in 1:nfold) {
        x <- folds[[i]]
        print(paste("Fold ",i))
        traindt <- ftmat[-x,]
        traingroup <- phenogroup[-x]
        wts <- as.vector(1/(table(traingroup)[traingroup]/length(traingroup)))
        testdt <- ftmat[x,]
        testgroup <- phenogroup[x]
        cvfit <- cv.glmnet(traindt, traingroup, family="binomial", alpha=alpha, weights=wts, lambda=round(10^(seq(-3.5, 1.5, 0.05)),digit=5))
        print(cvfit$lambda.min)
        autocoef <- coef(cvfit,s = "lambda.min")
        coefdt <- data.frame(autocoef@Dimnames[[1]][autocoef@i + 1], autocoef@x)
        colnames(coefdt) <- c("feature", paste0("I",i,"R",curcycle))
	if (i == 1) {
            popcoef <- coefdt
        } else {
            popcoef <- merge(popcoef, coefdt, by=c("feature"),all=TRUE)
        }
	predsdev <- assess.glmnet(cvfit, newx=testdt, newy=testgroup)$deviance
	predsdevl <- c(predsdevl, predsdev[1])
	predsmse <- assess.glmnet(cvfit, newx=testdt, newy=testgroup)$mse
	predsmsel <- c(predsmsel, predsmse[1])
	predsauc <- assess.glmnet(cvfit, newx=testdt, newy=testgroup)$auc
        predsaucl <- c(predsaucl, predsauc[1])
    }
    print(summary(predsdevl))
    print(summary(predsmsel))
    print(summary(predsaucl))
    predsdevdt <- t(data.frame(c(paste0("DEV_R",curcycle),predsdevl)))
    write.table(predsdevdt,mseaucout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    predsmsedt <- t(data.frame(c(paste0("MSE_R",curcycle),predsmsel)))
    write.table(predsmsedt,mseaucout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    predsaucdt <- t(data.frame(c(paste0("AUC_R",curcycle),predsaucl)))
    write.table(predsaucdt,mseaucout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    return(popcoef)
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


FeatureGLMCV <- function(ftmat, lowvarfilter, sidgroup, phenogroup, alpha, mylambda, nfold, curcycle, predresout, devmseout) {
    folds <- createFolds(phenogroup, k=nfold)
    predsmsel <- c()
    predsdevl <- c()
    alltruegroup <- phenogroup[unlist(folds)]
    allpredgroup <- c()
    predres <- data.frame()

    for (i in 1:nfold) {
        x <- folds[[i]]
        print(paste("Fold ",i))
        traindt <- ftmat[-x,]
        traingroup <- phenogroup[-x]
        wts <- as.vector(1/(table(traingroup)[traingroup]/length(traingroup)))
        testdt <- ftmat[x,]
        testgroup <- phenogroup[x]
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

RunGLMAssessReplicates <- function(ftmat, lowvarfilter, flagindex, alpha, mylambda, nfold, numrep, coefout, coefsumout, predresout, devmseout, perfoutfig, selected.feat=NA) {
    nfold <- as.numeric(nfold)
    numrep <- as.numeric(numrep)
    
    alpha <- as.numeric(alpha)
    mylambda <- as.numeric(mylambda)
    glmparam <- paste("#GLM-alpha: ",alpha,", lambda: ",mylambda,"; using selected features: ",selected.feat)
    print(glmparam)
    write.table(glmparam,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)

    sidgroup <- row.names(ftmat)
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
	popcoef <- FeatureGLMCV(ftmat, lowvarfilter, sidgroup, phenogroup, alpha, mylambda, nfold, curcycle=j, predresout, devmseout)

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

GenerateGLMmodel <- function(ftmat, flagindex, alpha, selected.feat=NA, modeldir, modelname) {
    alpha <- as.numeric(alpha)
    glmparam <- paste("Generating GLM model - alpha: ",alpha,"; using specified features: ",selected.feat)
    print(glmparam)
    
    sidgroup <- row.names(ftmat)
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
    finalrgml <- glmnet(ftmat, phenogroup, family="binomial", weights=wts, alpha=alpha, lambda=round(10^(seq(-3.5, 1.5, 0.05)),digit=5))
    print("writing GLM model")
    savedmodel <- paste0(modeldir, "/", modelname, ".rds")
    saveRDS(finalrgml, savedmodel)
}

ApplyGLMmodel <- function(predmat, selected.feat=NA, mylambda, modeldir, modelname, outfh, outcoef, outfig) {
    mylambda <- as.numeric(mylambda)
    sidgroup <- row.names(predmat)
    annoid <- sapply(strsplit(sidgroup, split=':', fixed=TRUE), function(x) (x[2]))
    annogroup <- sapply(strsplit(sidgroup, split=':', fixed=TRUE), function(x) (x[2]))

    if (!is.na(selected.feat)) {
        featidx <- fread(selected.feat, header=TRUE, sep="\t",data.table=FALSE)
        fidx <- featidx[order(featidx$Index),]
        sfeat <- fidx$Index
        predmat <- predmat[,colnames(predmat) %in% sfeat]
	coefheader <- t(data.frame(sfeat))
	write.table(coefheader,outcoef,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    } else {
        flagidx <- fread(flagindexfh, header=TRUE, sep="\t", data.table=FALSE)
        fidx <- flagidx[flagidx$FlagIndex==1,]
        fidx <- fidx[order(fidx$Index),]
	coefheader <- t(data.frame(fidx$Index))
	write.table(coefheader,outcoef,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
    }
    print(paste0("Classifying ",dim(predmat)[1]," subjects using ", modelname, ".rds"))
    write.table(paste0("Using ", modeldir, "/", modelname, ".rds\n"),outfh,row.names=FALSE,col.names=FALSE,quote=FALSE)
    usemodel <- readRDS(paste0(modeldir, "/", modelname, ".rds"))
    nsam <- nrow(predmat)
    predlist <- sapply(1:nsam, function(x) {
        predsample <- predmat[x,,drop=FALSE]
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



#############for input from seqmonk
CombineForwardReverseMatrix <- function(frlist, indexfh, outfh) {
    
    samlist <- fread(frlist, header=FALSE, sep="\t", colClasses=c("character", "character"), data.table=FALSE)
    idx <- LoadProbeIndex(indexfh)
    dtmat <- data.frame()
    for (i in 1:nrow(samlist)) {
         print(paste("reading in", samlist[i,1]))
         fwdfh <- fread(samlist[i,1], header=TRUE, sep="\t", data.table=FALSE)
	 colnames(fwdfh) <- gsub(".cov.gz","",colnames(fwdfh))
	 idxfwdfh <- merge(idx, fwdfh, by=c("Probe","Chromosome","Start","End"), all.x=TRUE)
	 idxfwdfh <- idxfwdfh[order(idxfwdfh$Index),]
	 idxfwdfh <- idxfwdfh[,-c(1:4,6:14)]
	 idxfwdfh_long <- melt(idxfwdfh, id.vars=c("Index"))
	 colnames(idxfwdfh_long) <- c("Index","Sample","Fwdval")

	 print(paste("reading in", samlist[i,2]))
	 revfh <- fread(samlist[i,2], header=TRUE, sep="\t", data.table=FALSE)
	 colnames(revfh) <- gsub(".cov.gz","",colnames(revfh))
	 idxrevfh <- merge(idx, revfh, by=c("Probe","Chromosome","Start","End"), all.x=TRUE)
	 idxrevfh <- idxrevfh[order(idxrevfh$Index),]
	 idxrevfh <- idxrevfh[,-c(1:4,6:14)]
	 idxrevfh_long <- melt(idxrevfh, id.vars=c("Index"))
	 colnames(idxrevfh_long) <- c("Index","Sample","Revval")

	 dt_long <- merge(idxfwdfh_long, idxrevfh_long, by=c("Index","Sample"))
	 #dt_long$Total <- dt_long$Fwdval+dt_long$Revval
	 #dt_long$Avgmc <- dt_long$Fwdval/dt_long$Total
         dt_long$Mval <- log2((dt_long$Fwdval+1)/(dt_long$Revval+1))

	 print(head(dt_long))
	 #dtrt <- dcast(melt(dt_long, id.vars=c("Index","Sample")), Index~Sample+variable)
	 #print(head(dtrt))

         dt <- dcast(dt_long, Index ~ Sample, value.var="Mval")
	 dt <- dt[order(dt$Index),]
	 dx <- t(dt[,-1])
	 colnames(dx) <- paste0("I",dt$Index)
	 dtmat <- rbind(dtmat, dx)
    }

    write.table(dtmat, outfh, sep="\t", row.names=TRUE, quote=FALSE)
}
