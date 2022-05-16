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

PlotSumChange(sampleinfo,outfigdir,cntoption)
