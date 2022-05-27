library(data.table)
library(reshape2)
library(ggplot2)

args <- commandArgs(TRUE)
sumfh <- args[1]
outfig <- args[2]

plotoverlapP <- function(sumfh,outfig) {
    fh <- fread(sumfh,header=FALSE,sep="\t",data.table=FALSE)
    colnames(fh) <- c("Index","Material","Type")
    FreshTissue <- fh[fh$Material=="FreshTissue",]$Index
    FFPE <- fh[fh$Material=="FFPE",]$Index
    Blood <- fh[fh$Material=="Blood",]$Index
    cfDNA_1stTrimester <- fh[fh$Material=="cfDNA_1stTrimester",]$Index
    cfDNA_atDiagnosis <- fh[fh$Material=="cfDNA_atDiagnosis",]$Index

    total <- 31289
    # Fresh tissue vs Blood
    groupT <- length(FreshTissue)
    groupB <- length(Blood)
    overlapTB <- length(Reduce(intersect,list(FreshTissue,Blood)))
    pTB <- phyper(overlapTB-1, groupB, total-groupB, groupT,lower.tail=FALSE)
    
    # Fresh tissue vs cfDNA diagnosis
    groupT <- length(FreshTissue)
    groupC <- length(cfDNA_atDiagnosis)
    overlapTC <- length(Reduce(intersect,list(FreshTissue,cfDNA_atDiagnosis)))
    pTC <- phyper(overlapTC-1, groupC, total-groupC, groupT,lower.tail=FALSE)
    
    # Blood vs cfDNA diagnosis
    groupB <- length(Blood)
    groupC <- length(cfDNA_atDiagnosis)
    overlapBC <- length(Reduce(intersect,list(Blood,cfDNA_atDiagnosis)))
    pBC <- phyper(overlapBC-1, groupC, total-groupC, groupB,lower.tail=FALSE)

    pmat <- matrix(c(0,pTB,pTC,pTB,0,pBC,pTC,pBC,0),nrow=3)
    colnames(pmat) <- row.names(pmat) <- c("Placenta","Blood","cfDNA at diagnosis")
    plotdt <- melt(pmat)
    plotdt$group <- cut(plotdt$value,breaks=c(-1,0,0.05,1))
    print(head(plotdt))
    pdf(outfig,height=5,width=5.5)
    #p1 <- ggplot(data=plotdt,aes(x=Var1, y=Var2, fill=value))+geom_tile(color="black")+theme_bw()+theme(axis.text.x=element_text(size=14,angle=45,vjust=0.5,color="black"),axis.text.y=element_text(size=14,color="black"),axis.title=element_blank(),legend.position="none")+geom_text(aes(Var1, Var2, label = formatC(value,digits=3,format="e")),color="red",size=4)+scale_fill_gradientn(limits=c(0,1),colors = hcl.colors(200, "Greens"))
    p1 <- ggplot(data=plotdt,aes(x=Var1, y=Var2, fill=group))+geom_tile(color="black")+theme_bw()+theme(axis.text.x=element_text(size=14,angle=45,vjust=0.5,color="black"),axis.text.y=element_text(size=14,color="black"),axis.title=element_blank(),legend.position="none")+geom_text(aes(Var1, Var2, label = formatC(value,digits=3,format="e")),color="red",size=4)+scale_fill_manual(breaks=levels(plotdt$group),values=c("#004616","#93CA8B","#F6FBF4"))
    print(p1)
    dev.off()
}

plotoverlapP(sumfh,outfig)
