library(data.table)
library(reshape2)
library(caret)
library(e1071)
library(pROC)
library(glmnet)

#args <- commandArgs(TRUE)

LoadProbeIndex <- function(indexfh) {
    idx <- fread(indexfh, header=TRUE, sep="\t", data.table=FALSE)
    return(idx)
}

GetCountMatrix <- function(samplelist, indexfh, cntoption="mval") {
    #' @title Read multiple samples into a data matrix
    #' @param samplelist A tab-separated file contains structured sample id and sample filename 'SAMPLE.ID\tSAMPLE.FILENAME'
    #' @param indexfh A tab-separated file contains capture information 'CHROM\tSTART\tEND\tINDEX\tPROBE'
    samlist <- fread(samplelist, header=TRUE, sep="\t", colClasses=c("character", "character"), data.table=FALSE)
    idx <- LoadProbeIndex(indexfh)
    #datalist <- lapply()
    dtmat <- data.frame()
    
    for (i in 1:nrow(samlist)) {
        print(paste("reading in", samlist[i,1]))
	curfh <- fread(samlist[i,1], header=FALSE, sep="\t", data.table=FALSE)
	colnames(curfh) <- c("Chromosome","Start","End","Index","mCount","umCount","Probe")
	curfh <- merge(idx,curfh,by=c("Chromosome","Start","End","Index","Probe"),all.x=TRUE)
	curfh <- curfh[order(curfh$Index),]
	
	if (cntoption == "mval") {
	    curfh$calval <- log2((curfh$mCount+1)/(curfh$umCount+1))
	} else if (cntoption == "betaval") {
	    curfh$calval <- curfh$mCount/(curfh$mCount+curfh$umCount)
	} else {
	    print("specify value to be calculated")
	}

	if (i == 1) {
	    dtmat <- data.frame(t(curfh$calval))
	} else {
	    dtmat <- rbind(dtmat, data.frame(t(curfh$calval)))
	}
    }
    colnames(dtmat) <- samlist$SAMPLE.ID
    return(dtmat)
}

CombineForwardReverseMatrix <- function(frlist, indexfh, outfh) {
    #for input from seqmonk count matrix
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