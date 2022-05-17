library(data.table)
library(limma)
library(ggplot2)
library(EnhancedVolcano)

args <- commandArgs(TRUE)
sampleinfo <- args[1]
inputdir <- args[2]
flagindexfh <- args[3]
cntoption <- args[4]
material <- args[5]
outprefix <- args[6]

LoadProbeIndex <- function(indexfh) {
    idx <- fread(indexfh, header=TRUE, sep="\t", data.table=FALSE)
    idx <- idx[order(idx$Index),]
    return(idx)
}

FlagFailedFeat <- function(inmat) {
    flagidx <- LoadProbeIndex(flagindexfh)
    keepindex <- flagidx$Index[flagidx$FlagIndex==1]
    ftmat <- inmat[,colnames(inmat) %in% keepindex]
    return(ftmat)
}

myvolcano <- function(mytoptable) {
    EnhancedVolcano(mytoptable,
       lab=rownames(mytoptable),
       ylab=bquote(~-Log[10]~ 'Adjusted' ~ italic('P')),
       x="logFC",
       y="adj.P.Val",
       title=NULL,
       subtitle=NULL,
       pCutoff=0.05,
       FCcutoff=1,
       ylim=c(0,max(-log10(mytoptable$adj.P.Val))+0.1),
       pointSize=1.5,
       labSize=0,
       shape=c(16, 16, 16, 16),
       caption=NULL,
       legendPosition='none',
       cutoffLineType='blank',
       hline=c(1.30103),
       vline=c(0),
       hlineWidth=2,
       vlineWidth=2,
       legendLabels=c('NS','Log2FC','Adj. p-value','Adj. p-value & Log2FC'),
       col = c("grey30", "cyan", "forestgreen", "red2"),
       colAlpha = 0.6,
       gridlines.minor = FALSE)
}

GetLimmaMatrix <- function(sampleinfo, inputdir, flagindexfh, cntoption, material, outprefix) {
    saminfo <- fread(sampleinfo, header=TRUE, sep="\t", data.table=FALSE)
    flagidx <- LoadProbeIndex(flagindexfh)
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
	    curfh$calval <- 100*curfh$Methylated/(curfh$Methylated+curfh$Unmethylated)
	} else {
	    print("specify value to be calculated")
	}
		
	if (i == 1) {
	    inmat <- matrix(curfh$calval,nrow=1)
	} else {
	    inmat <- rbind(inmat, matrix(curfh$calval,nrow=1))
	}
    }
    #in case of beta matrix, fill missing value
    colmedian <- apply(inmat,2,median,na.rm=TRUE)
    medianmat<-t(matrix(colmedian,ncol=nrow(inmat),nrow=ncol(inmat)))
    inmat[is.na(inmat)] <- medianmat[is.na(inmat)]

    colnames(inmat) <- flagidx$Index
    ftmat <- FlagFailedFeat(inmat)
    ftmat <- as.data.frame(ftmat)
    names(ftmat) <- paste0("I",colnames(ftmat))
    
    if (material == "Freshplacenta" || material == "FFPEplacenta" || material == "Blood") {
        Phenotype <- factor(saminfo$Phenotype)
	GA <- saminfo$GA

        design0 <- model.matrix(~ 0+Phenotype)
	myobj0 <- lmFit(t(ftmat),design0)
        contrastmat0 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design0))
	myfit0 <- contrasts.fit(myobj0, contrastmat0)
        myfit0 <- eBayes(myfit0)
        limmatable0 <- topTable(myfit0,number = ncol(ftmat))
        outfh0 <- paste0(outprefix,".",material,".",cntoption,".limma.tsv")
        write.table(limmatable0,outfh0,row.names=TRUE,sep="\t",quote=FALSE)
	outfig0 <- paste0(outprefix,".",material,".",cntoption,".limma.FC.adjP.volcano.pdf")
        pdf(outfig0,height=7,width=7)
        p0 <- myvolcano(limmatable0)
        print(p0)
        dev.off()

	design1 <- model.matrix(~ 0+Phenotype+GA)
	myobj1 <- lmFit(t(ftmat),design1)
	contrastmat1 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design1))
	myfit1 <- contrasts.fit(myobj1, contrastmat1)
	myfit1 <- eBayes(myfit1)
	limmatable1 <- topTable(myfit1,number = ncol(ftmat))
	outfh1 <- paste0(outprefix,".",material,".",cntoption,".limma.factor_GA.tsv")
	write.table(limmatable1,outfh1,row.names=TRUE,sep="\t",quote=FALSE)
	outfig1 <- paste0(outprefix,".",material,".",cntoption,".limma.factor_GA.FC.adjP.volcano.pdf")
	pdf(outfig1,height=7,width=7)
	p1 <- myvolcano(limmatable1)
	print(p1)
	dev.off()
    } else if (material == "cfDNAatDiagnosis" || material=="OxcfDNAatDiagnosis" || material == "cfDNAfirstT") {
        Phenotype <- factor(saminfo$Phenotype)
	GA <- saminfo$GA
	absMean <- apply(ftmat, 1, function(x) abs(mean(x)))

        design0 <- model.matrix(~ 0+Phenotype)
	myobj0 <- lmFit(t(ftmat),design0)
        contrastmat0 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design0))
	myfit0 <- contrasts.fit(myobj0, contrastmat0)
        myfit0 <- eBayes(myfit0)
        limmatable0 <- topTable(myfit0,number = ncol(ftmat))
        outfh0 <- paste0(outprefix,".",material,".",cntoption,".limma.tsv")
        write.table(limmatable0,outfh0,row.names=TRUE,sep="\t",quote=FALSE)
	outfig0 <- paste0(outprefix,".",material,".",cntoption,".limma.FC.adjP.volcano.pdf")
        pdf(outfig0,height=7,width=7)
        p0 <- myvolcano(limmatable0)
        print(p0)
        dev.off()

	design1 <- model.matrix(~ 0+Phenotype+GA)
	myobj1 <- lmFit(t(ftmat),design1)
	contrastmat1 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design1))
	myfit1 <- contrasts.fit(myobj1, contrastmat1)
	myfit1 <- eBayes(myfit1)
	limmatable1 <- topTable(myfit1,number = ncol(ftmat))
	outfh1 <- paste0(outprefix,".",material,".",cntoption,".limma.factor_GA.tsv")
	write.table(limmatable1,outfh1,row.names=TRUE,sep="\t",quote=FALSE)
	outfig1 <- paste0(outprefix,".",material,".",cntoption,".limma.factor_GA.FC.adjP.volcano.pdf")
	pdf(outfig1,height=7,width=7)
	p1 <- myvolcano(limmatable1)
	print(p1)
	dev.off()

        design2 <- model.matrix(~ 0+Phenotype+GA+absMean)
        myobj2 <- lmFit(t(ftmat),design2)
        contrastmat2 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design2))
        myfit2 <- contrasts.fit(myobj2, contrastmat2)
        myfit2 <- eBayes(myfit2)
        limmatable2 <- topTable(myfit2,number = ncol(ftmat))
        outfh2 <- paste0(outprefix,".",material,".",cntoption,".limma.factor_GA_absMean.tsv")
        write.table(limmatable2,outfh2,row.names=TRUE,sep="\t",quote=FALSE)
        outfig2 <- paste0(outprefix,".",material,".",cntoption,".limma.factor_GA_absMean.FC.adjP.volcano.pdf")
        pdf(outfig2,height=7,width=7)
        p2 <- myvolcano(limmatable2)
        print(p2)
        dev.off()
    }
}

GetLimmaMatrix(sampleinfo, inputdir, flagindexfh, cntoption, material, outprefix)