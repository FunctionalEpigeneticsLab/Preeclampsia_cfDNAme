#!/usr/bin/Rscript

list.of.packages <- c("data.table", "ggplot2", "dplyr", "cowplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)

args <- commandArgs(TRUE)
pilefh <- args[1]
outprefix <- args[2]

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

PlotFragmC(pilefh, outprefix)