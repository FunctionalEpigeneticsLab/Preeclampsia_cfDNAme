library(data.table)
library(reshape2)
library(caret)
library(e1071)
library(pROC)
library(glmnet)

LoadLockModelIndex <- function(indexfh) {
    idx <- fread(indexfh, header=TRUE, sep="\t", data.table=FALSE)
    idx <- idx[order(idx$Index),]
    return(idx)
}

GetValidCountMatrix <- function(validmat, modelindex) {
    indt <- fread(validmat, header=TRUE, sep="\t", data.table=FALSE)
    modelidx <- LoadLockModelIndex(modelindex)
    stidx <- modelidx$Index
    ftdt <- indt[,colnames(indt) %in% stidx]
    ftdt <- as.matrix(ftdt)
    colnames(ftdt) <- stidx
    row.names(ftdt) <- indt$SubjectIDPheno
    return(ftdt)
}

SecondRoundGLMCV <- function(ftmat, sidgroup, phenogroup, alpha, mylambda, nfold, curcycle, predresout, devmseout) {
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
        plot(glmperfm,print.auc=FALSE,print.auc.y=0.6-0.03*(curcycle-1),print.auc.x=0.3)
    } else {
        plot(glmperfm,add=TRUE,print.auc=FALSE,print.auc.y=0.6-0.03*(curcycle-1),print.auc.x=0.3)
    }
    return(popcoef)
}

RunGLM2RoundReplicates <- function(ftmat, modelindex, alpha, mylambda, nfold, numrep, outprefix) {
    nfold <- as.numeric(nfold)
    numrep <- as.numeric(numrep)

    modelidx <- LoadLockModelIndex(modelindex)
    
    alpha <- as.numeric(alpha)
    mylambda <- as.numeric(mylambda)
    glmparam <- paste0("#GLM: alpha.",alpha,".lambda.",mylambda)
    print(glmparam)

    sidgroup <- row.names(ftmat)
    phenogroup <- sapply(strsplit(sidgroup, split=':', fixed=TRUE), function(x) (x[2]))
    phenogroup <- as.factor(phenogroup)
	
    if (length(levels(phenogroup)) == 2) {
        levels(phenogroup) <- list("Ctrl"=0, "Case"=1)
    } else {
        stop("number of group greater than 2")
    }
    
    if (alpha==1) {
	coefout <- paste0(outprefix,".auto.lasso.coef.out")
	predresout <- paste0(outprefix,".auto.lasso.predres.out")
	devmseout <- paste0(outprefix,".auto.lasso.devmse.out")
	coefsumout <- paste0(outprefix,".auto.lasso.coef.freq75per.sum")
	perfoutfig <- paste0(outprefix,".auto.lasso.perf.pdf")
	write.table(glmparam,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
	
	coefheader <- t(data.frame(colnames(ftmat)))
	write.table(coefheader,coefout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)	
	repcoef <- data.frame()
	pdf(perfoutfig,width=9,height=9)
	for (j in 1:numrep) {
	    #set.seed(j)
	    pfline <- paste("Replicate ",j)
	    print(pfline)
	    predresheader <- paste("##CV Rep",j)
	    write.table(predresheader,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
	    #write.table(predresheader,devmseout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
	    popcoef <- SecondRoundGLMCV(ftmat, sidgroup, phenogroup, alpha, mylambda, nfold, curcycle=j, predresout, devmseout)

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
    } else if (alpha==0) {
	for (aa in c(0.997, 0.9968, 0.9965, 0.996, 0.9955, 0.995, 0.9945, 0.994)) {
	    coefout <- paste0(outprefix,".",aa,"quantile.ridge.coef.out")
	    predresout <- paste0(outprefix,".",aa,"quantile.ridge.predres.out")
	    devmseout <- paste0(outprefix,".",aa,"quantile.ridge.devmse.out")
	    coefsumout <- paste0(outprefix,".",aa,"quantile.ridge.coef.freq75per.out")
	    perfoutfig <- paste0(outprefix,".",aa,"quantile.ridge.perf.pdf")
	    write.table(glmparam,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
		
	    submodelidx <- modelidx[abs(modelidx$ENcoef) > quantile(abs(modelidx$ENcoef),aa),]$Index
	    print(submodelidx)
	    subftmat <- ftmat[,colnames(ftmat) %in% submodelidx]
	    coefheader <- t(data.frame(colnames(subftmat)))
	    write.table(coefheader,coefout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
	    repcoef <- data.frame()
	    pdf(perfoutfig,width=9,height=9)
	    for (j in 1:numrep) {
		#set.seed(j)
		pfline <- paste("Replicate ",j)
		print(pfline)
		predresheader <- paste("##CV Rep",j)
		write.table(predresheader,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
		#write.table(predresheader,devmseout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
		popcoef <- SecondRoundGLMCV(subftmat, sidgroup, phenogroup, alpha, mylambda, nfold, curcycle=j, predresout, devmseout)

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
    }
}

RunGLM2RoundRandomSampReplicates <- function(ftmat, modelindex, alpha, mylambda, nfold, numrep, outprefix) {
    nfold <- as.numeric(nfold)
    numrep <- as.numeric(numrep)

    modelidx <- LoadLockModelIndex(modelindex)
    
    alpha <- as.numeric(alpha)
    mylambda <- as.numeric(mylambda)
    glmparam <- paste0("#GLM: alpha.",alpha,".lambda.",mylambda)
    print(glmparam)

    sidgroup <- row.names(ftmat)
    phenogroup <- sapply(strsplit(sidgroup, split=':', fixed=TRUE), function(x) (x[2]))
    phenogroup <- as.factor(phenogroup)
	
    if (length(levels(phenogroup)) == 2) {
        levels(phenogroup) <- list("Ctrl"=0, "Case"=1)
    } else {
        stop("number of group greater than 2")
    }
    
    if (alpha==1) {
	coefout <- paste0(outprefix,".auto.lasso.coef.out")
	predresout <- paste0(outprefix,".auto.lasso.predres.out")
	devmseout <- paste0(outprefix,".auto.lasso.devmse.out")
	coefsumout <- paste0(outprefix,".auto.lasso.coef.freq75per.sum")
	perfoutfig <- paste0(outprefix,".auto.lasso.perf.pdf")
	write.table(glmparam,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
	
	coefheader <- t(data.frame(colnames(ftmat)))
	write.table(coefheader,coefout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)	
	repcoef <- data.frame()
	pdf(perfoutfig,width=9,height=9)
	for (j in 1:numrep) {
	    #set.seed(j)
	    pfline <- paste("Replicate ",j)
	    print(pfline)
	    predresheader <- paste("##CV Rep",j)
	    write.table(predresheader,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
	    #write.table(predresheader,devmseout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
	    popcoef <- SecondRoundGLMCV(ftmat, sidgroup, phenogroup, alpha, mylambda, nfold, curcycle=j, predresout, devmseout)

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
    } else if (alpha==0) {
	for (aa in c(11)) {
	    coefout <- paste0(outprefix,".",aa,"prop.random.ridge.coef.out")
	    predresout <- paste0(outprefix,".",aa,"prop.random.ridge.predres.out")
	    devmseout <- paste0(outprefix,".",aa,"prop.random.ridge.devmse.out")
	    coefsumout <- paste0(outprefix,".",aa,"prop.random.ridge.coef.freq75per.out")
	    perfoutfig <- paste0(outprefix,".",aa,"prop.random.ridge.perf.pdf")
	    write.table(glmparam,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
		
	    submodelidx <- sort(sample(modelidx$Index, aa, replace=FALSE))
	    print(submodelidx)
	    subftmat <- ftmat[,colnames(ftmat) %in% submodelidx]
	    coefheader <- t(data.frame(colnames(subftmat)))
	    write.table(coefheader,coefout,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
	    repcoef <- data.frame()
	    pdf(perfoutfig,width=9,height=9)
	    for (j in 1:numrep) {
		#set.seed(j)
		pfline <- paste("Replicate ",j)
		print(pfline)
		predresheader <- paste("##CV Rep",j)
		write.table(predresheader,predresout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
		#write.table(predresheader,devmseout,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
		popcoef <- SecondRoundGLMCV(subftmat, sidgroup, phenogroup, alpha, mylambda, nfold, curcycle=j, predresout, devmseout)

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
    }
}
