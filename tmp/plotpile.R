#!/usr/bin/Rscript

list.of.packages <- c("data.table", "ggplot2", "dplyr", "cowplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)

args <- commandArgs(TRUE)
#pilefh <- args[1]
#outprefix <- args[2]

#testfh <- args[1]
#outfig <- args[2]
#outfh1 <- args[3]
#outfh2 <- args[4]

GetPileFrag <- function(infh) {
    fh <- fread(infh, header=FALSE, sep="\t", data.table=FALSE)
    fh <- fh[,c(1:7,13)]
    colnames(fh) <- c("chr", "start", "end", "size", "mcrate", "mcpos", "umcpos","type")
    #remove fragments with no CG site
    #subfh <- fh[!(is.na(fh$mcpos) & is.na(fh$umcpos)),]
    subfh <- fh[!(fh$mcpos=="NULL" & fh$umcpos=="NULL"),]
    dt <- subfh %>% filter(!(type==".")) %>% rowwise() %>% mutate(mcnum=ifelse(mcpos=="NULL", 0, length(unlist(strsplit(mcpos,split=","))))) %>% mutate(umcnum=ifelse(umcpos=="NULL", 0, length(unlist(strsplit(umcpos,split=","))))) %>% mutate(totalmc = mcnum+umcnum)
    outfh <- paste0(outprefix, ".filtered.tsv")
    write.table(dt, outfh, quote=FALSE, row.names=FALSE, sep="\t")
    return(subfh)
}

PlotFragmC <- function(infh, outprefix) {
    fh <- GetPileFrag(infh)
    outfig <- paste0(outprefix, ".frag.mc.pdf")
    fh$type[fh$type=="."] <- "Off"
    fh$type <- factor(fh$type, levels=c("Placenta", "Blood", "Mix", "Off"))
    fh$mcrate <- fh$mcrate*100
    ntype <- length(unique(fh$type))
    pdf(outfig, height=7.6, width=9.8)
    #p1 <- ggplot(fh, aes(x=mcrate)) + geom_histogram(aes(y = ..density..), size=0.3, color="black", binwidth=1, alpha=0.3, position="identity") + geom_density(aes(color=type), alpha=0.9)+theme(axis.text.y=element_text(color="black",face="bold"),axis.text.x=element_text(color="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme_bw()+scale_color_manual(values=rainbow(ntype+1))+facet_grid(cols=vars(type))+xlim(-1, 101)+labs(x="mC (%)", y="Density")
    p1 <- ggplot(fh, aes(x=mcrate, color=type, fill=type))+geom_density(alpha=0.9)+theme_bw()+theme(axis.text.y=element_text(color="black",face="bold"),axis.text.x=element_text(color="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none")+scale_color_manual(values=rainbow(ntype+1))+facet_grid(cols=vars(type))+xlim(-1, 101)+labs(x="mC (%)", y="Density")
    p2 <- ggplot(fh, aes(x=mcrate, color=type, fill=type))+geom_histogram(size=0.3, binwidth=1, alpha=0.9, position="identity")+theme_bw()+theme(axis.text.y=element_text(color="black",face="bold"),axis.text.x=element_text(color="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none")+scale_fill_manual(values=rainbow(ntype+1))+scale_color_manual(values=rainbow(ntype+1))+facet_grid(cols=vars(type))+xlim(-1, 101)+labs(x="mC (%)", y="Count")
    p3 <- plot_grid(p1, p2, align='v', ncol=1)
    print(p3)
    dev.off()
}

PlotFragmCnum0 <- function(testfh, outfig, outfh1, outfh2) {
    fh <- fread(testfh, header=TRUE, sep="\t", data.table=FALSE)
    fh <- fh[,c("size","mcrate","type","totalmc")]
    #subfh <- subset(fh, totalmc < 26)
    #subfh$totalmc <- as.factor(subfh$totalmc)
    #subfh$mcrate <- subfh$mcrate*100
    #subfh <-subset(fh, totalmc < 31)
    #subfh$mcrate <- subfh$mcrate*100
    fh$mcrate <- as.numeric(fh$mcrate)
    tt <- fh %>% group_by(type, mcrange=cut(mcrate,breaks=seq(-0.001,1.001,by=0.01)-.Machine$double.eps)) %>% tally()
    write.table(tt,outfh1,row.names=FALSE,quote=FALSE,sep="\t")
    dtgroup <- fh %>% group_by(type, totalmc)
    dtsum <- dtgroup %>% summarize(mcavg = mean(mcrate), mcsd = sd(mcrate))
    write.table(dtsum,outfh2,row.names=FALSE,quote=FALSE,sep="\t")
    subdtsum <- subset(dtsum, totalmc < 31)
    #pdf(outfig, height=8, width=6)
    #p1 <- ggplot(subfh,aes(x=totalmc, y=mcrate, color=type))+geom_boxplot(outlier.shape = NA)+theme_bw()+facet_grid(rows=vars(type))+scale_color_manual(values=rainbow(3))
    pdf(outfig, height=5, width=11.5)
    p1 <- ggplot(subdtsum, aes(x=totalmc, y=mcavg, fill=type))+geom_bar(stat="identity",position=position_dodge(),color="black")+geom_errorbar(aes(ymin=mcavg-mcsd, ymax=mcavg+mcsd), width=0.2,position=position_dodge(0.9))+theme_bw()+scale_fill_manual(values=rainbow(3))
    print(p1)
    dev.off()
}

PlotFragmCnum <- function(testfh, outfig, outfh1, outfh2) {
    fh <- fread(testfh, header=TRUE, sep="\t", data.table=FALSE)
    fh <- fh[,c("size","mcrate","type","totalmc")]
    #subfh <- subset(fh, totalmc < 26)
    #subfh$totalmc <- as.factor(subfh$totalmc)
    #subfh$mcrate <- subfh$mcrate*100
    #subfh <-subset(fh, totalmc < 31)
    #subfh$mcrate <- subfh$mcrate*100
    fh$mcrate <- as.numeric(fh$mcrate)
    tt <- fh %>% group_by(type, mcrate) %>% tally()
    write.table(tt,outfh1,row.names=FALSE,quote=FALSE,sep="\t")
    dtsum1 <- fh %>% group_by(type, totalmc) %>% summarize(mcavg = mean(mcrate), mcsd = sd(mcrate))
    dtsum2 <- fh %>% group_by(type, totalmc) %>% tally()
    dtsum <- merge(dtsum1, dtsum2, by=c("type","totalmc"))
    write.table(dtsum,outfh2,row.names=FALSE,quote=FALSE,sep="\t")
    subdtsum <- subset(dtsum, totalmc < 31)
    #pdf(outfig, height=8, width=6)
    #p1 <- ggplot(subfh,aes(x=totalmc, y=mcrate, color=type))+geom_boxplot(outlier.shape = NA)+theme_bw()+facet_grid(rows=vars(type))+scale_color_manual(values=rainbow(3))
    pdf(outfig, height=8, width=11.5)
    p1 <- ggplot(subdtsum, aes(x=totalmc, y=mcavg, fill=type))+geom_bar(stat="identity",position=position_dodge(),color="black")+geom_errorbar(aes(ymin=mcavg-mcsd, ymax=mcavg+mcsd), width=0.2,position=position_dodge(0.9))+theme_bw()+scale_fill_manual(values=rainbow(3))+labs(x="totalC")
    p2 <- ggplot(subdtsum, aes(x=totalmc, y=n, fill=type))+geom_bar(stat="identity",position=position_dodge(),color="black")+theme_bw()+scale_fill_manual(values=rainbow(3))+labs(x="totalC")
    p3 <- plot_grid(p1, p2, align='v', ncol=1)
    print(p3)
    dev.off()
}

#PlotFragmCnum(testfh, outfig, outfh1, outfh2)
#PlotFragmC(pilefh, outprefix)

samples.file <- args[1]
outfh <- args[2]
outfig <- args[3]
covfh <- args[4]

PlotSumSam <- function(samples.file, outfh, outfig, covfh) {
    cov <- fread(covfh, header=TRUE, sep="\t", data.table=FALSE)
    samples.list <- fread(samples.file, sep="\t", header=TRUE, colClasses=c("character"), data.table=FALSE)
    samples.df <- samples.list$filename
    datalist <- lapply(samples.df, function(x) {
        df <- fread(x, header=TRUE, sep="\t")
	cal1 <- df %>% group_by(type) %>% summarise(expm=sum(exp(df$mcrate)*df$n/sum(n)))
	cal2 <- df %>% filter((mcrate==0 | mcrate==1)) %>% group_by(type) %>% summarise(biprop=n[mcrate==1]/n[mcrate==0])
	caldt <- merge(cal1, cal2, by=c("type"))
	caldt$sample <- unlist(strsplit(basename(x),"[.]"))[1]
	return(caldt)
    })
    dtmatrix <- Reduce(function(x, y) {rbind(x, y)}, datalist)
    write.table(dtmatrix, outfh, row.names=FALSE, quote=FALSE, sep="\t")
    dtplot1 <- reshape2::melt(dtmatrix, id=c("sample","type"))
    dtplot2 <- reshape2::dcast(dtplot1, sample+variable ~ type, value.var="value")
    #dtplot2 <- subset(dtplot2, sample != "tBT2260_217N1" & sample != "GC114391_715N1")
    dtplot2$time <- substr(dtplot2$sample,(nchar(dtplot2$sample)+1)-2,nchar(dtplot2$sample))
    dtplot2$time[dtplot2$time=="_T"] = "Tissue"
    dtplot2$time[dtplot2$time=="PE"] = "FFPE"
    subdt1 <- subset(dtplot2,variable=="expm")
    subdt1 <- merge(subdt1, cov, by=c("sample"))
    subdt2 <- subset(dtplot2,variable=="biprop")
    subdt2 <- merge(subdt2, cov, by=c("sample"))
    pdf(outfig, height=8.5, width=7)
    p1 <- ggplot(subdt1, aes(x=Blood, y=Placenta, color=time, size=meanCov))+geom_point(alpha=0.5)+theme_bw()+theme(axis.text.y=element_text(color="black",face="bold"),axis.text.x=element_text(color="black"),panel.grid.major=element_blank())+scale_color_manual(values=rainbow(4))+scale_size_continuous(range = c(0.5,4))+labs(title="expm")
    p2 <- ggplot(subdt2, aes(x=Blood, y=Placenta, color=time, size=meanCov))+geom_point(alpha=0.5)+theme_bw()+theme(axis.text.y=element_text(color="black",face="bold"),axis.text.x=element_text(color="black"),panel.grid.major=element_blank())+scale_color_manual(values=rainbow(4))+scale_size_continuous(range = c(0.5,4))+labs(title="biprop")
    p3 <- plot_grid(p1, p2, align='v', ncol=1)
    print(p3)
    dev.off()
}

PlotSumSam(samples.file, outfh, outfig, covfh)