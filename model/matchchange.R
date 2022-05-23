source("/data/leuven/323/vsc32370/script/cfDNA-methylation-processing-and-AUC/model/glmeval.R")
library(data.table)
library(reshape2)
library(ggplot2)
library(viridis)
library(cowplot)
library(scales)
library(dplyr)

args <- commandArgs(TRUE)
sampleinfo <- args[1]
inputdir <- args[2]
flagindexfh <- args[3]
cntoption <- args[4]
outfigdir <- args[5]
regionlist <- args[6]

LoadProbeIndex <- function(indexfh) {
    idx <- fread(indexfh,header=TRUE,sep="\t",data.table=FALSE)
    idx <- idx[order(idx$Index),]
    return(idx)
}

PlotMatchSampleChange <- function(sampleinfo,inputdir,flagindexfh,cntoption="betaval",outfigdir) {
    saminfo <- fread(sampleinfo,header=TRUE,sep="\t",data.table=FALSE)
    flagidx <- LoadProbeIndex(flagindexfh)
    inmat <- matrix()
    keepindex <- flagidx$Index[flagidx$FlagIndex==1]
    for (i in 1:nrow(saminfo)) {
        print(paste("reading in subject ", saminfo[i,"Individual"]))
        samfile1 <- paste0(saminfo[i,"N2bsID"],".mavg.count.merge.tsv")
        samfh1 <- file.path(inputdir,samfile1)
        curfh1 <- fread(samfh1,header=TRUE,sep="\t",data.table=FALSE)
        #print(head(curfh1))
        curfh1 <- merge(flagidx,curfh1,by=c("Chromosome","Start","End","Index","Probe"),all.x=TRUE)
        #print(head(curfh1))
	curfh1 <- curfh1[order(curfh1$Index),]
	curfh1 <- curfh1[curfh1$Index %in% keepindex,]

        samfile2 <- paste0(saminfo[i,"N2oxbsID"],".mavg.count.merge.tsv")
       	samfh2 <- file.path(inputdir,samfile2)
       	curfh2 <- fread(samfh2,header=TRUE,sep="\t",data.table=FALSE)
        curfh2 <- merge(flagidx,curfh2,by=c("Chromosome","Start","End","Index","Probe"),all.x=TRUE)
        curfh2 <- curfh2[order(curfh2$Index),]
       	curfh2 <- curfh2[curfh2$Index %in% keepindex,]

        if (cntoption == "mval") {
	    curfh1$mval1 <- log2((curfh1$Methylated+1)/(curfh1$Unmethylated+1))
	    curfh2$mval2 <- log2((curfh2$Methylated+1)/(curfh2$Unmethylated+1))
	    curfh1 <- curfh1[,c("Chromosome","Start","End","Index","Probe","FlagIndex","FlagPla","mval1")]
	    curfh2 <- curfh2[,c("Chromosome","Start","End","Index","Probe","FlagIndex","FlagPla","mval2")]
	    combfh <- merge(curfh1,curfh2,by=c("Chromosome","Start","End","Index","Probe","FlagIndex","FlagPla"))
	    combfh <- combfh[order(combfh$Index),]

	    combfh$mchange <- combfh$mval2-combfh$mval1
	    combfh <- combfh[combfh$FlagIndex==1,]
	    #combchange <- combfh[combfh$mchange > 1.5 | combfh$mchange < -1.5,]

	    outfh <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".all.tsv")
	    write.table(combfh,outfh,quote=FALSE,sep="\t",row.names=FALSE)
	    bloodcombfh <- subset(combfh,FlagPla==0)
	    placombfh <- subset(combfh,FlagPla==1)
	    
	    outfig0 <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".all.pdf")
	    pdf(outfig0,height=7,width=7.8)
	    p0 <- ggplot(combfh,aes(x=mval1,y=mval2))+geom_bin2d(bins=100)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14))+scale_fill_continuous(type="viridis")+labs(x="1st sampling methylation level",y="2nd sampling methylation level")+scale_x_continuous(limits = c(-10, 5))+scale_y_continuous(limits = c(-10, 5))
	    print(p0)
	    dev.off()

            outfig1 <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".blood.pdf")
            pdf(outfig1,height=7,width=7.8)
            p1 <- ggplot(bloodcombfh,aes(x=mval1,y=mval2))+geom_bin2d(bins=100)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14))+scale_fill_continuous(type="viridis")+labs(x="1st sampling methylation level",y="2nd sampling methylation level")+scale_x_continuous(limits = c(-10, 5))+scale_y_continuous(limits = c(-10, 5))
            print(p1)
            dev.off()

            outfig2 <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".placenta.pdf")
            pdf(outfig2,height=7,width=7.8)
            p2 <- ggplot(placombfh,aes(x=mval1,y=mval2))+geom_bin2d(bins=100)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14))+scale_fill_continuous(type="viridis")+labs(x="1st sampling methylation level",y="2nd sampling methylation level")+scale_x_continuous(limits = c(-10, 5))+scale_y_continuous(limits = c(-10, 5))
            print(p2)
            dev.off()
	    
	} else if (cntoption == "betaval") {
	    curfh1$betaval1 <- 100*curfh1$Methylated/(curfh1$Methylated+curfh1$Unmethylated)
            curfh2$betaval2 <- 100*curfh2$Methylated/(curfh2$Methylated+curfh2$Unmethylated)
            curfh1 <- curfh1[,c("Chromosome","Start","End","Index","Probe","FlagIndex","FlagPla","betaval1")]
            curfh2 <- curfh2[,c("Chromosome","Start","End","Index","Probe","FlagIndex","FlagPla","betaval2")]
            combfh <- merge(curfh1,curfh2,by=c("Chromosome","Start","End","Index","Probe","FlagIndex","FlagPla"))
            print(head(combfh))
	    combfh <- combfh[order(combfh$Index),]

	    combfh$mchange <- combfh$betaval2-combfh$betaval1
	    combfh <- combfh[combfh$FlagIndex==1,]
	    
            outfh <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".all.tsv")
            write.table(combfh,outfh,quote=FALSE,sep="\t",row.names=FALSE)
	    bloodcombfh	<- subset(combfh,FlagPla==0)
            placombfh <- subset(combfh,FlagPla==1)
	    
	    outfig0 <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".all.pdf")
            pdf(outfig0,height=7,width=7.8)
            #p0 <- ggplot(combfh,aes(x=betaval1,y=betaval2))+geom_bin2d(bins=100)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14))+scale_fill_continuous(type="viridis")+labs(x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
            p0 <- ggplot(combfh,aes(x=betaval1,y=betaval2))+geom_bin2d(bins=150)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14,color="black"))+scale_fill_gradient(low="darkblue",high="cyan",trans="log10")+labs(fill = "Density (log)",x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
            print(p0)
            dev.off()

            outfig1 <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".blood.pdf")
            pdf(outfig1,height=7,width=7.8)
            p1 <- ggplot(bloodcombfh,aes(x=betaval1,y=betaval2))+geom_bin2d(bins=150)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14,color="black"))+scale_fill_gradient(low="darkblue",high="cyan",trans="log10")+labs(fill = "Density (log)",x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
            print(p1)
            dev.off()

            outfig2 <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".placenta.pdf")
            pdf(outfig2,height=7,width=7.8)
            p2 <- ggplot(placombfh,aes(x=betaval1,y=betaval2))+geom_bin2d(bins=150)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14,color="black"))+scale_fill_gradient(low="darkblue",high="cyan",trans="log10")+labs(fill = "Density (log)",x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
            print(p2)
            dev.off()
	} else {
	    print("specify value to be calculated")
	}
    }  
}

#PlotMatchSampleChange(sampleinfo,inputdir,flagindexfh,cntoption,outfigdir)

PlotSumChange <- function(sampleinfo,outfigdir,cntoption) {
    saminfo <- fread(sampleinfo,header=TRUE,sep="\t",data.table=FALSE)
    sumdiff <- data.frame()
    for (i in 1:nrow(saminfo)) {
        mdifffh <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".all.tsv")
        print(mdifffh)
        curfh <- fread(mdifffh,header=TRUE,sep="\t",data.table=FALSE)
        curfh$mchange[is.na(curfh$mchange)] <- 0
        curbindiff <- curfh %>% group_by(group=cut(abs(mchange),breaks=c(-1,1,5,10,20,50,100))) %>% summarise(n=n())
        curbindiff$Individual <- saminfo[i,"Individual"]
        curbindiff$Phenotype <- saminfo[i,"Phenotype"]
        sumdiff <- rbind(sumdiff,curbindiff)
    }
    sumdiff$Phenotype <- factor(sumdiff$Phenotype,levels=c("Ctrl","Case"))
    outfig <- paste0(outfigdir,"/N2.bs.oxbs.methyldiff.",cntoption,".all.sum.boxplot.new.pdf")
    pdf(outfig,width=7,height=5.8)
    p1 <- ggplot(sumdiff,aes(x=group,y=n,fill=Phenotype))+geom_boxplot(position=position_dodge(1))+geom_dotplot(binwidth=0.1,dotsize=2,binaxis='y',stackdir='center',position=position_dodge(1))+theme_bw()+theme(axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black",angle=45,vjust=0.5),legend.text=element_text(size=14),axis.title=element_text(size=14,color="black"))+scale_fill_manual(values=c("#999999","#E69F00"))+labs(x="Methylation difference between BS and OxBS",y="Region count")+scale_x_discrete(labels=c("< 1%","1% - 5%","5% - 10%","10% - 20%","20% - 50%", "> 50%"))
    print(p1)
    dev.off()
}

#PlotSumChange(sampleinfo,outfigdir,cntoption)

PlotHyperHypoChange <- function(sampleinfo,outfigdir,cntoption,regionlist) {
    saminfo <- fread(sampleinfo,header=TRUE,sep="\t",data.table=FALSE)
    regions <- fread(regionlist,header=FALSE,sep="\t",data.table=FALSE)
    colnames(regions) <- c("Index","Material","Mtype")
    sumdiff <- data.frame()
    hypersigdiff <- data.frame()
    hyposigdiff <- data.frame()
    for (i in 1:nrow(saminfo)) {
        mdifffh <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".all.tsv")
        print(mdifffh)
        curfh <- fread(mdifffh,header=TRUE,sep="\t",data.table=FALSE)
        curfh$mchange[is.na(curfh$mchange)] <- 0
	regionfh <- merge(regions,curfh,by=c("Index"),all.x=TRUE)
	hyperfh <- subset(regionfh,Mtype=="Hypermethylated")
	hypofh <- subset(regionfh,Mtype=="Hypomethylated")
	hyperfh$Individual <- saminfo[i,"Individual"]
	hyperfh$Phenotype <- saminfo[i,"Phenotype"]
	hypersigdiff <- rbind(hypersigdiff,hyperfh)
	
	hypofh$Individual <- saminfo[i,"Individual"]
	hypofh$Phenotype <- saminfo[i,"Phenotype"]
	hyposigdiff <- rbind(hyposigdiff,hypofh)

        outfig0 <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".unadjP.hyper.pdf")
        pdf(outfig0,height=5.6,width=5.6)
        #p0 <- ggplot(combfh,aes(x=betaval1,y=betaval2))+geom_bin2d(bins=150)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14,color="black"))+scale_fill_gradient(low="darkblue",high="cyan",trans="log10")+labs(fill = "Density (log)",x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
	hypercorr <- round(cor(hyperfh$betaval1,hyperfh$betaval2,method=c("pearson")),3)
	p0 <- ggplot(hyperfh,aes(x=betaval1,y=betaval2))+geom_point(size=1,color="darkblue",alpha=0.5)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14,color="black"))+geom_text(x=25,y=90,label=paste0("Pearson Correlation: ",hypercorr),color="darkred")+labs(x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
        print(p0)
	dev.off()
	
	outfig1 <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".unadjP.hypo.pdf")
        pdf(outfig1,height=5.6,width=5.6)
        hypocorr <- round(cor(hypofh$betaval1,hypofh$betaval2,method=c("pearson")),3)
	p1 <- ggplot(hypofh,aes(x=betaval1,y=betaval2))+geom_point(size=1,color="darkblue",alpha=0.5)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14,color="black"))+geom_text(x=25,y=90,label=paste0("Pearson Correlation: ",hypocorr),color="darkred")+labs(x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
        print(p1)
        dev.off()
	    
	hyperbindiff <- hyperfh %>% group_by(group=cut(abs(mchange),breaks=c(-1,1,5,10,20,50,100))) %>% summarise(n=n()) %>% mutate(freq = n / sum(n))
	hypobindiff <- hypofh %>% group_by(group=cut(abs(mchange),breaks=c(-1,1,5,10,20,50,100))) %>% summarise(n=n()) %>% mutate(freq = n / sum(n))
	hypobindiff$methyl <- "Hypomethylated"
	hyperbindiff$methyl <- "Hypermethylated"
	curbindiff <- rbind(hypobindiff,hyperbindiff)
        curbindiff$Individual <- saminfo[i,"Individual"]
        curbindiff$Phenotype <- saminfo[i,"Phenotype"]
        sumdiff <- rbind(sumdiff,curbindiff)
    }
    sumdiff$methyl <- factor(sumdiff$methyl,levels=c("Hypomethylated","Hypermethylated"))
    sumdiff$freq <- 100*sumdiff$freq
    outfig <- paste0(outfigdir,"/N2.bs.oxbs.methyldiff.",cntoption,".all.sum.boxplot.N2BShypohyper.pdf")
    pdf(outfig,width=7,height=5.8)
    pp <- ggplot(sumdiff,aes(x=group,y=freq,fill=methyl))+geom_boxplot(position=position_dodge(1))+geom_dotplot(binwidth=0.1,dotsize=2,binaxis='y',stackdir='center',position=position_dodge(1))+theme_bw()+theme(axis.text.y=element_text(size=14,color="black"),axis.text.x=element_text(size=14,color="black",angle=45,vjust=0.5),legend.text=element_text(size=14),axis.title=element_text(size=14,color="black"))+scale_fill_manual(values=c("#e41a1c","#377eb8"))+labs(x="Methylation difference between BS and OxBS",y="Percentage of regions (%)")+scale_x_discrete(labels=c("< 1%","1% - 5%","5% - 10%","10% - 20%","20% - 50%", "> 50%"))

    sumhypersigdiff <- hypersigdiff %>% group_by(Index) %>% summarize(Mean1=mean(betaval1,na.rm=TRUE),Mean2=mean(betaval2,na.rm=TRUE))
    outfig91 <- paste0(outfigdir,"/N2.bs.oxbs.methyldiff.",cntoption,".unadjP.hyper.allmean.pdf")
    pdf(outfig91,width=5.6,height=5.6)
    allhypercorr <- round(cor(sumhypersigdiff$Mean1,sumhypersigdiff$Mean2,method=c("pearson")),3)
    p91 <- ggplot(sumhypersigdiff,aes(x=Mean1,y=Mean2))+geom_point(size=1,color="darkblue",alpha=0.5)+theme_bw()+theme(axis.text=element_text(size=14,color="black"),axis.title=element_text(size=14,color="black"))+geom_text(x=25,y=90,label=paste0("Pearson Correlation: ",hypercorr),color="darkred")+labs(x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
    print(p91)
    dev.off()

    sumhyposigdiff <- hyposigdiff %>% group_by(Index) %>% summarize(Mean1=mean(betaval1,na.rm=TRUE),Mean2=mean(betaval2,na.rm=TRUE))
    outfig92 <- paste0(outfigdir,"/N2.bs.oxbs.methyldiff.",cntoption,".unadjP.hypo.allmean.pdf")
    pdf(outfig92,width=5.6,height=5.6)
    allhypocorr <- round(cor(sumhyposigdiff$Mean1,sumhyposigdiff$Mean2,method=c("pearson")),3)
    p92 <- ggplot(sumhyposigdiff,aes(x=Mean1,y=Mean2))+geom_point(size=1,color="darkblue",alpha=0.5)+theme_bw()+theme(axis.text=element_text(size=14,color="black"),axis.title=element_text(size=14,color="black"))+geom_text(x=25,y=90,label=paste0("Pearson Correlation: ",hypocorr),color="darkred")+labs(x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
    print(p92)
    dev.off()

    sumhypersigdiff2 <- hypersigdiff %>% group_by(Index,Phenotype) %>% summarize(Mean1=mean(betaval1,na.rm=TRUE),Mean2=mean(betaval2,na.rm=TRUE))
    outfig93 <- paste0(outfigdir,"/N2.bs.oxbs.methyldiff.",cntoption,".unadjP.hyper.phenomean.pdf")
    pdf(outfig93,width=9.6,height=5.6)
    ctlhyper <- subset(sumhypersigdiff2,Phenotype=="Ctrl")
    cashyper <- subset(sumhypersigdiff2,Phenotype=="Case")
    ctlhypercorr <- round(cor(ctlhyper$Mean1,ctlhyper$Mean2,method=c("pearson")),3)
    cashypercorr <- round(cor(cashyper$Mean1,cashyper$Mean2,method=c("pearson")),3)

    corrtext <- data.frame(Phenotype=c("Case","Ctrl"),label=c(paste0("Pearson Correlation: ",cashypercorr),paste0("Pearson Correlation: ",ctlhypercorr)))
    p93 <- ggplot(sumhypersigdiff2,aes(x=Mean1,y=Mean2))+geom_point(size=1,color="darkblue",alpha=0.5)+theme_bw()+theme(axis.text=element_text(size=14,color="black"),axis.title=element_text(size=14,color="black"),strip.text.x=element_text(size=14,color="black"),legend.position="none")+geom_text(data=corrtext,x=25,y=85,aes(label=label))+labs(x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))+facet_grid(.~Phenotype)
    print(p93)
    dev.off()

    sumhyposigdiff2 <- hyposigdiff %>% group_by(Index,Phenotype) %>% summarize(Mean1=mean(betaval1,na.rm=TRUE),Mean2=mean(betaval2,na.rm=TRUE))
    outfig94 <- paste0(outfigdir,"/N2.bs.oxbs.methyldiff.",cntoption,".unadjP.hypo.phenomean.pdf")
    pdf(outfig94,width=9.6,height=5.6)
    ctlhypo <- subset(sumhyposigdiff2,Phenotype=="Ctrl")
    cashypo <- subset(sumhyposigdiff2,Phenotype=="Case")
    ctlhypocorr <- round(cor(ctlhypo$Mean1,ctlhypo$Mean2,method=c("pearson")),3)
    cashypocorr <- round(cor(cashypo$Mean1,cashypo$Mean2,method=c("pearson")),3)
    corrtext <- data.frame(Phenotype=c("Case","Ctrl"),label=c(paste0("Pearson Correlation: ",cashypocorr),paste0("Pearson Correlation: ",ctlhypocorr)))
    p94 <- ggplot(sumhyposigdiff2,aes(x=Mean1,y=Mean2))+geom_point(size=1,color="darkblue",alpha=0.5)+theme_bw()+theme(axis.text=element_text(size=14,color="black"),axis.title=element_text(size=14,color="black"),strip.text.x=element_text(size=14,color="black"),legend.position="none")+geom_text(data=corrtext,x=25,y=85,aes(label=label))+labs(x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))+facet_grid(.~Phenotype)
    print(p94)
    dev.off()
    
}

#PlotHyperHypoChange(sampleinfo,outfigdir,cntoption,regionlist)

PlotSigChange <- function(sampleinfo,outfigdir,cntoption,regionlist) {
    saminfo <- fread(sampleinfo,header=TRUE,sep="\t",data.table=FALSE)
    regions <- fread(regionlist,header=FALSE,sep="\t",data.table=FALSE)
    colnames(regions) <- c("Index","Material","Mtype")
	sigdiff <- data.frame()

    for (i in 1:nrow(saminfo)) {
        mdifffh <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".all.tsv")
        print(mdifffh)
        curfh <- fread(mdifffh,header=TRUE,sep="\t",data.table=FALSE)
        curfh$mchange[is.na(curfh$mchange)] <- 0
	regionfh <- merge(regions,curfh,by=c("Index"),all.x=TRUE)
	regionfh$Individual <- saminfo[i,"Individual"]
	regionfh$Phenotype <- saminfo[i,"Phenotype"]
	sigdiff <- rbind(sigdiff,regionfh)

        outfig0 <- paste0(outfigdir,"/",saminfo[i,"Phenotype"],".",saminfo[i,"Individual"],".N2.bs.oxbs.methyldiff.",cntoption,".sig.adjP.pdf")
        pdf(outfig0,height=5.6,width=5.6)
        
	allcorr <- round(cor(regionfh$betaval1,regionfh$betaval2,method=c("pearson")),3)
	p0 <- ggplot(regionfh,aes(x=betaval1,y=betaval2))+geom_point(size=1,color="darkblue",alpha=0.5)+theme_bw()+theme(legend.text=element_text(size=14),text=element_text(size=14,color="black"))+geom_text(x=25,y=90,label=paste0("Pearson Correlation: ",allcorr),color="darkred")+labs(x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
        print(p0)
	dev.off()
    }

    sumsigdiff <- sigdiff %>% group_by(Index) %>% summarize(Mean1=mean(betaval1,na.rm=TRUE),Mean2=mean(betaval2,na.rm=TRUE))
    outfig91 <- paste0(outfigdir,"/N2.bs.oxbs.methyldiff.",cntoption,".sig.adjP.allmean.pdf")
    pdf(outfig91,width=5.6,height=5.6)
    allcorr <- round(cor(sumsigdiff$Mean1,sumsigdiff$Mean2,method=c("pearson")),3)
    p91 <- ggplot(sumsigdiff,aes(x=Mean1,y=Mean2))+geom_point(size=1,color="darkblue",alpha=0.5)+theme_bw()+theme(axis.text=element_text(size=14,color="black"),axis.title=element_text(size=14,color="black"))+geom_text(x=25,y=90,label=paste0("Pearson Correlation: ",allcorr),color="darkred")+labs(x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))
    print(p91)
    dev.off()

    sumsigdiff2 <- sigdiff %>% group_by(Index,Phenotype) %>% summarize(Mean1=mean(betaval1,na.rm=TRUE),Mean2=mean(betaval2,na.rm=TRUE))
    outfig93 <- paste0(outfigdir,"/N2.bs.oxbs.methyldiff.",cntoption,".sig.adjP.phenomean.pdf")
    pdf(outfig93,width=9.6,height=5.6)
    ctldt <- subset(sumsigdiff2,Phenotype=="Ctrl")
    casdt <- subset(sumsigdiff2,Phenotype=="Case")
    ctlcorr <- round(cor(ctldt$Mean1,ctldt$Mean2,method=c("pearson")),3)
    cascorr <- round(cor(casdt$Mean1,casdt$Mean2,method=c("pearson")),3)
    corrtext <- data.frame(Phenotype=c("Case","Ctrl"),label=c(paste0("Pearson Correlation: ",cascorr),paste0("Pearson Correlation: ",ctlcorr)))
    p93 <- ggplot(sumsigdiff2,aes(x=Mean1,y=Mean2))+geom_point(size=1,color="darkblue",alpha=0.5)+theme_bw()+theme(axis.text=element_text(size=14,color="black"),axis.title=element_text(size=14,color="black"),strip.text.x=element_text(size=14,color="black"),legend.position="none")+geom_text(data=corrtext,x=25,y=85,aes(label=label))+labs(x="BS methylation level",y="OxBS methylation level")+scale_x_continuous(limits = c(-1, 100))+scale_y_continuous(limits = c(-1, 100))+facet_grid(.~Phenotype)
    print(p93)
    dev.off()
}

PlotSigChange(sampleinfo,outfigdir,cntoption,regionlist)