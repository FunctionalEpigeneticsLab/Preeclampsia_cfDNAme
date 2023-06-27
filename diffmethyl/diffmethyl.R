library(data.table)
library(limma)
library(ggplot2)
library(EnhancedVolcano)
library(VennDiagram)
library(eulerr)
library(plyr)

# Differential methylation analysis

args <- commandArgs(TRUE)
sampleinfo <- args[1]
inputdir <- args[2]
flagindexfh <- args[3]
cntoption <- args[4]
autosomeonly <- args[5]
normalization <- args[6]
material <- args[7]
outprefix <- args[8]

LoadProbeIndex <- function(indexfh) {
    idx <- fread(indexfh, header=TRUE, sep="\t", data.table=FALSE)
    idx <- idx[order(idx$Index),]
    return(idx)
}

#FlagFailedFeat <- function(inmat) {
#    flagidx <- LoadProbeIndex(flagindexfh)
#    keepindex <- flagidx$Index[flagidx$FlagIndex==1]
#    ftmat <- inmat[,colnames(inmat) %in% keepindex]
#    return(ftmat)
#}

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
            normfactor <- apply(ftmat, 1, function(x) abs(mean(x)))
	    } else if (cntoption == "avgme") {
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

make_volcano <- function(result.limma,adj = FALSE,Annotate = TRUE,FDRcutoff = 0.05) {
    maxfc <- max(result.limma$logFC)+0.1
    myRange <- c(-maxfc,maxfc)
    yaxlim <- c(0,max(-log10(result.limma$P.Value))+0.1)
    
    if (adj == FALSE) {
        colors <- densCols(result.limma$logFC, -log10(result.limma$P.Value), nbin = 800, colramp = colorRampPalette(c("grey85", "black")))
        colors2 <- densCols(subset(result.limma, subset =  P.Value < 0.05 & logFC < 0)$logFC, -log10(subset(result.limma, subset =  P.Value < 0.05 & logFC < 0)$P.Value), nbin=800, colramp = colorRampPalette(c("#52BDEC", "#0B1E57")))
        colors3 <- densCols(subset(result.limma, subset =  P.Value < 0.05 & logFC > 0)$logFC, -log10(subset(result.limma, subset =  P.Value < 0.05 & logFC > 0)$P.Value), nbin=800, colramp = colorRampPalette(c("#FF665F", "#f71302")))
        p1 <- ggplot(result.limma, aes(x = logFC, y = -log10(P.Value)))+geom_point(col = colors)+theme_bw()+xlim(myRange)+ylim(yaxlim)+geom_point(data = subset(result.limma, subset= P.Value < 0.05 & logFC < 0), aes(x = logFC, y = -log10(P.Value)), col = colors2)+geom_point(data = subset(result.limma, subset =  P.Value < 0.05 & logFC > 0), aes(x = logFC, y = -log10(P.Value)), col = colors3) +geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") + xlab('Fold change PE vs Control') + ylab(bquote(~-Log[10]~ italic('P'))) +theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text = element_text(size = 15))
      
        if (Annotate == TRUE){
            a = nrow(subset(result.limma, subset = P.Value   < 0.05 & logFC < 0))
            b = nrow(subset(result.limma, subset = P.Value   < 0.05 & logFC > 0))
            c = nrow(subset(result.limma, subset = adj.P.Val < FDRcutoff & logFC < 0))
            d = nrow(subset(result.limma, subset = adj.P.Val < FDRcutoff & logFC > 0))
            p1 <- p1+annotate("text", y = yaxlim[2],x = myRange[2], label = "Unadjusted P-Value",hjust = "right", vjust = "top",size = 7.5, fontface = "bold")+annotate("text", y = yaxlim[2]*(5/6),    x = myRange[2], label = "Adjusted P-Value",hjust = "right", vjust = "top",size = 7.5, fontface = "bold") + annotate("text", y = yaxlim[2]*(5.65/6), x = myRange[2], label = paste0("N=",format(a, big.mark = ",", scientific = FALSE)), hjust = "right", vjust = "top", size = 5.5, color = "#52BDEC", fontface = "bold")+annotate("text", y = yaxlim[2]*(5.4/6),  x = myRange[2], label = paste0("N=",format(b, big.mark = ",", scientific = FALSE)), hjust = "right", vjust = "top", size = 5.5, color = "#FF665F", fontface = "bold")+annotate("text", y = yaxlim[2]*(4.65/6), x = myRange[2], label = paste0("N=",format(c, big.mark = ",", scientific = FALSE)), hjust = "right", vjust = "top", size = 5.5, color = "#52BDEC", fontface = "bold")+annotate("text", y = yaxlim[2]*(4.4/6),  x = myRange[2], label = paste0("N=",format(d, big.mark = ",", scientific = FALSE)), hjust = "right", vjust = "top", size = 5.5, color = "#FF665F", fontface = "bold")
        } else if(Annotate == FALSE) {
           p1 <- p1
        }
    } else if (sum(result.limma$adj.P.Val < FDRcutoff) < 50) {
        colors <- densCols(result.limma$logFC, -log10(result.limma$adj.P.Val), nbin = 800, colramp = colorRampPalette(c("grey85", "black")))
        p1 <- ggplot(result.limma, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(col = colors) + theme_bw() + xlim(myRange) + ylim(yaxlim)+geom_point(data = subset(result.limma, subset =  adj.P.Val < FDRcutoff & logFC < 0), aes(x = logFC, y = -log10(adj.P.Val)), col = "#0759ed") +geom_point(data = subset(result.limma, subset =  adj.P.Val < FDRcutoff & logFC > 0), aes(x = logFC, y = -log10(adj.P.Val)), col = "#f71302")+geom_hline(yintercept = -log10(0.05), col="black", linetype = "dashed") + xlab('Fold change PE vs Control')
    } else {
        colors <- densCols(result.limma$logFC, -log10(result.limma$adj.P.Val), nbin = 800, colramp = colorRampPalette(c("grey85", "black")))
        colors2 <- densCols(subset(result.limma, subset =  adj.P.Val < FDRcutoff & logFC < 0)$logFC, -log10(subset(result.limma,subset =  adj.P.Val < FDRcutoff & logFC < 0)$adj.P.Val), nbin=800, colramp = colorRampPalette(c("#93b3ed", "#0759ed")))
        colors3 <- densCols(subset(result.limma, subset =  adj.P.Val < FDRcutoff & logFC > 0)$logFC, -log10(subset(result.limma,subset =  adj.P.Val < FDRcutoff & logFC > 0)$adj.P.Val), nbin=800, colramp = colorRampPalette(c("#ed9993", "#A53D3D")))
        p1 <- ggplot(result.limma, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(col = colors) + theme_bw() + xlim(myRange) + ylim(yaxlim)+geom_point(data = subset(result.limma, subset =  adj.P.Val < FDRcutoff & logFC < 0), aes(x = logFC, y = -log10(adj.P.Val)), col = colors2) +geom_point(data = subset(result.limma, subset =  adj.P.Val < FDRcutoff & logFC > 0), aes(x = logFC, y = -log10(adj.P.Val)), col = colors3)+geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") + xlab('Fold change PE vs Control')
    }
    return(p1)
}

GetLimmaMatrix <- function(sampleinfo, inputdir, flagindexfh, cntoption, autosomeonly, normalization, material, outprefix) {
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
            curfh$calval <- curfh$MePer
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
    ftmat <- FlagFailedFilterFeat(inmat, autosomeonly)

    ftmat <- as.data.frame(ftmat)
    #names(ftmat) <- paste0("I",colnames(ftmat))
    names(ftmat) <- colnames(ftmat)
    ftmat <- ftmat[,apply(ftmat,2,function(x) !any(is.na(x)))]
    ftmat <- ftmat[,colSums(ftmat!=0)>0]
    #print(dim(ftmat))
    ftmat <- NormalizeCountMatrix(ftmat,cntoption,normalization)
    
    if (material == "Freshplacenta" || material == "FFPEplacenta" || material == "Mixplacenta"|| material == "Blood") {
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
	limmatable0$Index <- row.names(limmatable0)
	fetchlimmatable0 <- merge(limmatable0,flagidx,by=c("Index"),all.x=TRUE)
	fetchlimmatable0 <- fetchlimmatable0[order(fetchlimmatable0$adj.P.Val),]
        outfh0 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.tsv")
        write.table(fetchlimmatable0,outfh0,row.names=FALSE,sep="\t",quote=FALSE)
	
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

        outfig02 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.FC.unadjP.protocol.volcano.pdf")
	pdf(outfig02,height=7,width=7)
        p02 <- make_volcano(limmatable0,adj = FALSE,Annotate = TRUE,FDRcutoff = 0.05)
	print(p02)
	dev.off()

	design1 <- model.matrix(~ 0+Phenotype+GA)
	myobj1 <- lmFit(t(ftmat),design1)
	contrastmat1 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design1))
	myfit1 <- contrasts.fit(myobj1, contrastmat1)
	myfit1 <- eBayes(myfit1)
	limmatable1 <- topTable(myfit1,number = ncol(ftmat))
	limmatable1$Index <- row.names(limmatable1)
        fetchlimmatable1 <- merge(limmatable1,flagidx,by=c("Index"),all.x=TRUE)
        fetchlimmatable1 <- fetchlimmatable1[order(fetchlimmatable1$adj.P.Val),]
	outfh1 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_GA.tsv")
	write.table(fetchlimmatable1,outfh1,row.names=FALSE,sep="\t",quote=FALSE)
	
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

    } else if (material == "cfDNAatDiagnosis" || material=="OxcfDNAatDiagnosis" || material == "cfDNAfirstT" || material == "cfDNAfirstT_valid") {
        Phenotype <- factor(saminfo$Phenotype)
	GA <- saminfo$GA
	SeqDep <- saminfo$MMeanDep
	MeanMeth <- apply(ftmat, 1, function(x) mean(x))
    Batch <- saminfo$Batch

        design0 <- model.matrix(~ 0+Phenotype)
	myobj0 <- lmFit(t(ftmat),design0)
        contrastmat0 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design0))
	myfit0 <- contrasts.fit(myobj0, contrastmat0)
        myfit0 <- eBayes(myfit0)
        limmatable0 <- topTable(myfit0,number = ncol(ftmat))
        outfh0 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.tsv")
	limmatable0$Index <- row.names(limmatable0)
        fetchlimmatable0 <- merge(limmatable0,flagidx,by=c("Index"))
        fetchlimmatable0 <- fetchlimmatable0[order(fetchlimmatable0$adj.P.Val),]
        
        write.table(fetchlimmatable0,outfh0,row.names=FALSE,sep="\t",quote=FALSE)
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

        outfig02 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.FC.unadjP.protocol.volcano.pdf")
        pdf(outfig02,height=7,width=7)
        p02 <- make_volcano(limmatable0,adj = FALSE,Annotate = TRUE,FDRcutoff = 0.05)
        print(p02)
        dev.off()

	design1 <- model.matrix(~ 0+Phenotype+Batch)
	myobj1 <- lmFit(t(ftmat),design1)
	contrastmat1 <- makeContrasts(PhenotypeCase-PhenotypeCtrl, levels=colnames(design1))
	myfit1 <- contrasts.fit(myobj1, contrastmat1)
	myfit1 <- eBayes(myfit1)
	limmatable1 <- topTable(myfit1,number = ncol(ftmat))
	limmatable1$Index <- row.names(limmatable1)
        fetchlimmatable1 <- merge(limmatable1,flagidx,by=c("Index"))
        fetchlimmatable1 <- fetchlimmatable1[order(fetchlimmatable1$adj.P.Val),]
	outfh1 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_Batch.tsv")
        write.table(fetchlimmatable1,outfh1,row.names=FALSE,sep="\t",quote=FALSE)

	outfig10 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_Batch.FC.adjP.volcano.pdf")
	pdf(outfig10,height=7,width=7)
	p10 <- adjvolcano(limmatable1)
	print(p10)
	dev.off()

        outfig11 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_Batch.FC.unadjP.volcano.pdf")
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
	limmatable3$Index <- row.names(limmatable3)
        fetchlimmatable3 <- merge(limmatable3,flagidx,by=c("Index"))
        fetchlimmatable3 <- fetchlimmatable3[order(fetchlimmatable3$adj.P.Val),]
        write.table(fetchlimmatable3,outfh3,row.names=FALSE,sep="\t",quote=FALSE)
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
	limmatable4$Index <- row.names(limmatable4)
        fetchlimmatable4 <- merge(limmatable4,flagidx,by=c("Index"))
        fetchlimmatable4 <- fetchlimmatable4[order(fetchlimmatable4$adj.P.Val),]
        outfh4 <- paste0(outprefix,".",material,".",cntoption,".",normalization,".limma.factor_SeqDep_MeanMeth.tsv")
        write.table(limmatable4,outfh4,row.names=FALSE,sep="\t",quote=FALSE)
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

GetLimmaMatrix(sampleinfo, inputdir, flagindexfh, cntoption, autosomeonly, normalization, material, outprefix)

display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

plotvenn <- function(sumfh,outfig1) {
    fh <- fread(sumfh,header=FALSE,sep="\t",data.table=FALSE)
    colnames(fh) <- c("Index","Material","Type")
    FreshTissue <- fh[fh$Material=="FreshTissue",]$Index
    FFPE <- fh[fh$Material=="FFPE",]$Index
    cfDNA_1stTrimester <- fh[fh$Material=="cfDNA_1stTrimester",]$Index
    cfDNA_atDiagnosis <- fh[fh$Material=="cfDNA_atDiagnosis",]$Index
    vlists <- list(FreshTissue,FFPE,cfDNA_atDiagnosis,cfDNA_1stTrimester)
    pdf(outfig1,height=7,width=7)
    p1 <- display_venn(
	    vlists,
	    category.names = c("Fresh Tissue" , "FFPE " , "cfDNA at diagnosis", "cfDNA 1st trimester"),
	    # Circles
	    lwd = 2,
	    lty = 'blank',
	    fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
	    # Numbers
	    cex = .9,
	    fontface = "italic",
	    # Set names
	    cat.cex = 1,
	    cat.fontface = "bold",
	    cat.default.pos = "outer",
	    cat.dist = c(0.055, 0.055, 0.1, 0.1)
    )
    print(p1)
    dev.off()
}

ploteuler <- function(sumfh,outfig2) {
    fh <- fread(sumfh,header=FALSE,sep="\t",data.table=FALSE)
    colnames(fh) <- c("Index","Material","Type")
    FreshTissue <- rep(TRUE,length(fh[fh$Material=="FreshTissue",]$Index))
    names(FreshTissue) <- fh[fh$Material=="FreshTissue",]$Index
    FFPE <- rep(TRUE,length(fh[fh$Material=="FFPE",]$Index))
    names(FFPE) <- fh[fh$Material=="FFPE",]$Index
    cfDNA_1stTrimester <- rep(TRUE,length(fh[fh$Material=="cfDNA_1stTrimester",]$Index))
    names(cfDNA_1stTrimester) <- fh[fh$Material=="cfDNA_1stTrimester",]$Index
    cfDNA_atDiagnosis <- rep(TRUE,length(fh[fh$Material=="cfDNA_atDiagnosis",]$Index))
    names(cfDNA_atDiagnosis) <- fh[fh$Material=="cfDNA_atDiagnosis",]$Index

    #vmat <- rbind.fill.matrix(t(FreshTissue), t(FFPE), t(cfDNA_atDiagnosis), t(cfDNA_1stTrimester))
    #vmat[is.na(vmat)] <- FALSE
    #eulmat <- t(vmat)
    #colnames(eulmat) <- c("Fresh Tissue" , "FFPE " , "cfDNA at diagnosis", "cfDNA 1st trimester")
    #eufit <- eulerr:euler(eulmat,shape="ellipse")

    vmat <- rbind.fill.matrix(t(FreshTissue), t(cfDNA_atDiagnosis), t(cfDNA_1stTrimester))
    vmat[is.na(vmat)] <- FALSE
    eulmat <- t(vmat)
    colnames(eulmat) <- c("Fresh Tissue" , "cfDNA at diagnosis", "cfDNA 1st trimester")
    eufit <- eulerr::euler(eulmat,shape="ellipse")

    pdf(outfig2,height=7,width=7)
    p1 <- plot(eufit, fills = c("#a6cee3", "#b2df8a", "#fdbf6f"), quantities = TRUE)
    print(p1)
    dev.off()
}
