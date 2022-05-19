library(data.table)
library(limma)
library(ggplot2)
library(EnhancedVolcano)

args <- commandArgs(TRUE)
sampleinfo <- args[1]
inputdir <- args[2]
flagindexfh <- args[3]
cntoption <- args[4]
normalization <- args[5]
material <- args[6]
outprefix <- args[7]

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

unadjvolcano <- function(mytoptable) {
    EnhancedVolcano(mytoptable,
       lab=rownames(mytoptable),
       xlab=bquote('Fold Change'),
       x="logFC",
       y="P.Value",
       title=NULL,
       subtitle=NULL,
       pCutoff=0.05,
       FCcutoff=1,
       ylim=c(0,max(-log10(mytoptable$P.Value))+0.1),
       pointSize=1.5,
       labSize=0,
       shape=c(16, 16, 16, 16),
       caption=NULL,
       legendPosition='none',
       cutoffLineType='blank',
       hline=0.05,
       vline=0,
       hlineWidth=1,
       vlineWidth=1,
       col = c("grey30", "cyan", "forestgreen", "red2"),
       colAlpha = 0.6,
       gridlines.minor = FALSE)
}

adjvolcano <- function(mytoptable) {
    EnhancedVolcano(mytoptable,
       lab=rownames(mytoptable),
       xlab=bquote('Fold Change'),
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
       hline=0.05,
       vline=0,
       hlineWidth=1,
       vlineWidth=1,
       legendLabels=c('NS','Log2FC','Adj. p-value','Adj. p-value & Log2FC'),
       col = c("grey30", "cyan", "forestgreen", "red2"),
       colAlpha = 0.6,
       gridlines.minor = FALSE)
}

GetLimmaMatrix <- function(sampleinfo, inputdir, flagindexfh, cntoption, normalization, material, outprefix) {
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
	} else if (cntoption == "log2betaval") {
	    curfh$calval <- log2(100*curfh$Methylated/(curfh$Methylated+curfh$Unmethylated)+0.001)
	} else {
	    print("specify value to be calculated")
	}
		
	if (i == 1) {
	    inmat <- matrix(curfh$calval,nrow=1)
	} else {
	    inmat <- rbind(inmat, matrix(curfh$calval,nrow=1))
	}
    }
    #in case of beta matrix, fill missing value, remove 0 variance column
    colmedian <- apply(inmat,2,median,na.rm=TRUE)
    medianmat<-t(matrix(colmedian,ncol=nrow(inmat),nrow=ncol(inmat)))
    inmat[is.na(inmat)] <- medianmat[is.na(inmat)]
    #inmat <- inmat[ , apply(inmat, 2, function(x) !any(is.na(x)))]
    #inmat[, colSums(inmat != 0) > 0]
    #print(dim(inmat))

    colnames(inmat) <- flagidx$Index
    ftmat <- FlagFailedFeat(inmat)
    ftmat <- as.data.frame(ftmat)
    names(ftmat) <- paste0("I",colnames(ftmat))
    ftmat <- ftmat[,apply(ftmat,2,function(x) !any(is.na(x)))]
    ftmat <- ftmat[,colSums(ftmat!=0)>0]
    #print(dim(ftmat))
    ftmat <- NormalizeCountMatrix(ftmat, normalization)
    
    if (material == "Freshplacenta" || material == "FFPEplacenta" || material == "Blood") {
        Phenotype <- factor(saminfo$Phenotype)
	GA <- saminfo$GA

        design0 <- model.matrix(~ 0+Phenotype)
	myobj0 <- lmFit(t(ftmat),design0)
        contrastmat0 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design0))
	myfit0 <- contrasts.fit(myobj0, contrastmat0)
        myfit0 <- eBayes(myfit0)
        limmatable0 <- topTable(myfit0,number = ncol(ftmat))
	#print(max(-log10(limmatable0$P.Value)))
	#print(max(-log10(limmatable0$adj.P.Val)))
        outfh0 <- paste0(outprefix,".",material,".",cntoption,".limma.tsv")
        write.table(limmatable0,outfh0,row.names=TRUE,sep="\t",quote=FALSE)
	
	#transform FC from betaval to log2 FC - not strictly the way limma deal with expression data
	outfig0 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.FC.adjP.volcano.pdf")
        pdf(outfig0,height=7,width=7)
        p00 <- adjvolcano(limmatable0)
        print(p00)
        dev.off()

        outfig01 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.FC.unadjP.volcano.pdf")
        pdf(outfig01,height=7,width=7)
        p01 <- unadjvolcano(limmatable0)
        print(p01)
        dev.off()

	design1 <- model.matrix(~ 0+Phenotype+GA)
	myobj1 <- lmFit(t(ftmat),design1)
	contrastmat1 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design1))
	myfit1 <- contrasts.fit(myobj1, contrastmat1)
	myfit1 <- eBayes(myfit1)
	limmatable1 <- topTable(myfit1,number = ncol(ftmat))
	outfh1 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_GA.tsv")
	write.table(limmatable1,outfh1,row.names=TRUE,sep="\t",quote=FALSE)
	
	outfig10 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_GA.FC.adjP.volcano.pdf")
	pdf(outfig10,height=7,width=7)
	p10 <- adjvolcano(limmatable1)
	print(p10)
	dev.off()

        outfig11 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_GA.FC.unadjP.volcano.pdf")
        pdf(outfig11,height=7,width=7)
        p11 <- unadjvolcano(limmatable1)
        print(p11)
        dev.off()
    } else if (material == "cfDNAatDiagnosis" || material=="OxcfDNAatDiagnosis" || material == "cfDNAfirstT") {
        Phenotype <- factor(saminfo$Phenotype)
	GA <- saminfo$GA
	SeqDep <- saminfo$MMeanDep
	MeanMeth <- apply(ftmat, 1, function(x) mean(x))

        design0 <- model.matrix(~ 0+Phenotype)
	myobj0 <- lmFit(t(ftmat),design0)
        contrastmat0 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design0))
	myfit0 <- contrasts.fit(myobj0, contrastmat0)
        myfit0 <- eBayes(myfit0)
        limmatable0 <- topTable(myfit0,number = ncol(ftmat))
        outfh0 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.tsv")
        write.table(limmatable0,outfh0,row.names=TRUE,sep="\t",quote=FALSE)
	outfig00 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.FC.adjP.volcano.pdf")
        pdf(outfig00,height=7,width=7)
        p00 <- adjvolcano(limmatable0)
        print(p00)
        dev.off()

        outfig01 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.FC.unadjP.volcano.pdf")
        pdf(outfig01,height=7,width=7)
        p01 <- unadjvolcano(limmatable0)
        print(p01)
        dev.off()

	design1 <- model.matrix(~ 0+Phenotype+MeanMeth)
	myobj1 <- lmFit(t(ftmat),design1)
	contrastmat1 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design1))
	myfit1 <- contrasts.fit(myobj1, contrastmat1)
	myfit1 <- eBayes(myfit1)
	limmatable1 <- topTable(myfit1,number = ncol(ftmat))
	outfh1 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_MeanMeth.tsv")
	write.table(limmatable1,outfh1,row.names=TRUE,sep="\t",quote=FALSE)
	outfig10 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_MeanMeth.FC.adjP.volcano.pdf")
	pdf(outfig10,height=7,width=7)
	p10 <- adjvolcano(limmatable1)
	print(p10)
	dev.off()

        outfig11 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_MeanMeth.FC.unadjP.volcano.pdf")
        pdf(outfig11,height=7,width=7)
        p11 <- unadjvolcano(limmatable1)
        print(p11)
        dev.off()

        #design2 <- model.matrix(~ 0+Phenotype+GA+MeanMeth)
        #myobj2 <- lmFit(t(ftmat),design2)
        #contrastmat2 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design2))
        #myfit2 <- contrasts.fit(myobj2, contrastmat2)
        #myfit2 <- eBayes(myfit2)
        #limmatable2 <- topTable(myfit2,number = ncol(ftmat))
        #outfh2 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_GA_MeanMeth.tsv")
        #write.table(limmatable2,outfh2,row.names=TRUE,sep="\t",quote=FALSE)
        #outfig20 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_GA_MeanMeth.FC.adjP.volcano.pdf")
        #pdf(outfig20,height=7,width=7)
        #p20 <- adjvolcano(limmatable2)
        #print(p20)
        #dev.off()

        #outfig21 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_GA_MeanMeth.FC.unadjP.volcano.pdf")
        #pdf(outfig21,height=7,width=7)
        #p21 <- unadjvolcano(limmatable2)
        #print(p21)
        #dev.off()

        design3 <- model.matrix(~ 0+Phenotype+SeqDep)
        myobj3 <- lmFit(t(ftmat),design3)
        contrastmat3 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design3))
        myfit3 <- contrasts.fit(myobj3, contrastmat3)
        myfit3 <- eBayes(myfit3)
        limmatable3 <- topTable(myfit3,number = ncol(ftmat))
        outfh3 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_SeqDep.tsv")
        write.table(limmatable3,outfh3,row.names=TRUE,sep="\t",quote=FALSE)
        outfig30 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_SeqDep.FC.adjP.volcano.pdf")
        pdf(outfig30,height=7,width=7)
        p30 <- adjvolcano(limmatable3)
        print(p30)
        dev.off()

        outfig31 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_SeqDep.FC.unadjP.volcano.pdf")
        pdf(outfig31,height=7,width=7)
        p31 <- unadjvolcano(limmatable3)
        print(p31)
	dev.off()

        design4 <- model.matrix(~ 0+Phenotype+SeqDep+MeanMeth)
        myobj4 <- lmFit(t(ftmat),design4)
        contrastmat4 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design4))
        myfit4 <- contrasts.fit(myobj4, contrastmat4)
        myfit4 <- eBayes(myfit4)
        limmatable4 <- topTable(myfit4,number = ncol(ftmat))
        outfh4 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_SeqDep_MeanMeth.tsv")
        write.table(limmatable4,outfh4,row.names=TRUE,sep="\t",quote=FALSE)
        outfig40 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_SeqDep_MeanMeth.FC.adjP.volcano.pdf")
        pdf(outfig40,height=7,width=7)
        p40 <- adjvolcano(limmatable4)
        print(p40)
        dev.off()

        outfig41 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_SeqDep_MeanMeth.FC.unadjP.volcano.pdf")
        pdf(outfig41,height=7,width=7)
        p41 <- unadjvolcano(limmatable4)
        print(p41)
        dev.off()
    }
}

GetLimmaMatrix(sampleinfo, inputdir, flagindexfh, cntoption, normalization, material, outprefix)
