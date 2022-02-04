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

GetCountMatrix <- function(sampleinfo, inputdir, indexfh, cntoption="mval", outmat) {
    #' @title Read multiple samples into a data matrix
    #' @param sampleinfo A tab-separated file contains structured information 'SubjectID\tPhenotype\tDescription\tGA\tSampleID\tMMeanDep\tMMedianDep\tBSFlag'
    #' @param inputdir A directory path to methylation count of all samples
    #' @param indexfh A tab-separated file contains capture information 'CHROM\tSTART\tEND\tINDEX\tPROBE'
    #' @param outmat An output file for combined matrix
    
    saminfo <- fread(sampleinfo, header=TRUE, sep="\t", colClasses=c("character","character","character","numeric","character","numeric","numeric","character"), data.table=FALSE)
    idx <- LoadProbeIndex(indexfh)
    #datalist <- lapply()
    inmat <- matrix()
    
    for (i in 1:nrow(saminfo)) {
        print(paste("reading in subject ", saminfo[i,"SubjectID"]))
	samfile <- paste0(saminfo[i,"SubjectID"], ".mavg.count.merge.tsv")
	samfh <- file.path(inputdir, samfile)
	curfh <- fread(samfh, header=TRUE, sep="\t", data.table=FALSE)
	curfh <- merge(idx,curfh,by=c("Chromosome","Start","End","Index","Probe"),all.x=TRUE)
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
    colnames(inmat) <- idx$Index
    rownames(inmat) <- paste0(saminfo$SubjectID,":",saminfo$Phenotype)
    write.table(inmat,outmat,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
    return(inmat)
}

FilterCountMatrixFeat <- function(sampleinfo, inputdir, indexfh, cntoption, outmat, flagindexfh) {
    inmat <- GetCountMatrix(sampleinfo, inputdir, indexfh, cntoption, outmat)
    flagidx <- fread(flagindexfh, header=TRUE, sep="\t", data.table=FALSE)
    keepindex <- flagidx$Index[flagidx$FlagIndex==1]
    ftmat <- inmat[,colnames(inmat) %in% keepindex]
    print(paste0("Keep ", dim(ftmat)[2], " targets for ", dim(ftmat)[1], " subjects"))
    return(ftmat)
}

AssessGLMLoo <- function(ftmat, alpha, lambda, outpred, outfig) {
    alpha <- as.numeric(alpha)
    lambda <- as.numeric(lambda)
    glmparam <- paste("GLM-alpha: ",alpha,", lambda: ",lambda)
    print(glmparam)
    
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
        rglm.model <- glmnet(traindt, traingroup, family="binomial", alpha=alpha, weights=wts, lambda=10^(seq(-3, 0.5, 0.05)))
        rglm.preds.type <- predict(rglm.model, testdt, type="class", s=lambda)
        rglm.preds.res <- predict(rglm.model, testdt, type="response", s=lambda)
        rglm.preds.coef <- predict(rglm.model, testdt, type="coefficient", s=lambda)
        alltruegroup <- c(alltruegroup, testgroup)
        allpredval <- c(allpredval, rglm.preds.res)
        write.table(paste0("predict ", sidgroup[x], "\t", testgroup, " as ", rglm.preds.type, "\t", rglm.preds.res),outpred,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
        return(allpredval)
    })
    
    glmperfm <- roc(response=phenogroup, predictor=loopreds, ci=TRUE)
    pdf(outfig, width=8, height=8)
    plot(glmperfm,print.auc=TRUE,print.auc.y=0.6,print.auc.x=0.3)
    dev.off()
    print(glmperfm$auc)
    print(glmperfm$ci)
    runacc <- table(phenogroup, ifelse(loopreds>0.5, "Case", "Ctrl"))
    print(runacc)
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