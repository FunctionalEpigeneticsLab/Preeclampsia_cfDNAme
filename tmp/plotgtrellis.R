list.of.packages <- c("data.table", "zoo", "gtrellis")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(data.table)
library(zoo)
library(gtrellis)
library(ggplot2)

args <- commandArgs(TRUE)
workdir <- args[1]
samplename <- args[2]
covthresh <- args[3]

PlotFeatureProfile <- function(workdir, samplename) {
    fh <- paste0(workdir, '/', samplename, '.R1_val_1_bismark_bt2_pe.sorted_by_name.deduplicated.bismark.cov.sorted.bed.ontar.bed')
    infh <- fread(fh, header=FALSE, sep="\t", data.table=FALSE)
    infh <- infh[,c(5:10)]
    colnames(infh) <- c("chromosome","start","end", "mpercent", "mcnt", "umcnt")
    infh$chromosome <- paste0("chr", infh$chromosome)
    infh$cov <- infh$mcnt + infh$umcnt
    infh <- subset(infh, cov > 10)
    covmax <- max(infh$cov)
    outfig <- paste0(fh,".dist.pdf")
    pdf(outfig, width=16.5, height=5.2)
    gtrellis_layout(n_track=2, track_axis=TRUE, border=FALSE, category="chr1", padding=unit(c(1, 1, 1, 1), "mm"), track_ylim=c(0,100,0,100), track_height = unit.c(unit(2, "null"), unit(2, "null")), track_ylab = c("Cov","MePer"), add_ideogram_track = TRUE, add_name_track = TRUE, axis_label_fontsize = 4, lab_fontsize = 5, name_fontsize = 6, title=samplename, title_fontsize=11)
    #add_points_track(infh, infh$cov, pch=1, size=unit(0.1, "mm"), gp=gpar(col="gray"))
    #add_points_track(infh, infh$mpercent, pch=1, size=unit(0.1, "mm"), gp=gpar(col="steelblue"))
    add_rect_track(infh, h1=infh$cov, h2=0, gp=gpar(col="gray"))
    add_rect_track(infh, h1=infh$mpercent, h2=0, gp=gpar(col="gray"))
    dev.off()
}

PlotPer <- function(workdir, samplename, covthresh) {
    autosomes <- c(1:22)
    covthresh <- as.numeric(covthresh)
    topfh <- paste0(workdir, "/CpG_OT_", samplename, ".R1_val_1_bismark_bt2_pe.sorted_by_name.deduplicated.count.sorted.ontar.bed")
    botfh <- paste0(workdir, "/CpG_OB_", samplename, ".R1_val_1_bismark_bt2_pe.sorted_by_name.deduplicated.count.sorted.ontar.bed")
    topinfh <- fread(topfh, header=FALSE, sep="\t", data.table=FALSE)
    botinfh <- fread(botfh, header=FALSE, sep="\t", data.table=FALSE)
    topinfh <- topinfh[,c(4:9)]
    botinfh <- botinfh[,c(4:9)]
    colnames(topinfh) <- c("target","chrom","start","end","topcov","topmethcnt")
    colnames(botinfh) <- c("target","chrom","start","end","botcov","botmethcnt")
    topinfh <- subset(topinfh, topinfh$topcov > covthresh)
    topinfh$topmethper <- topinfh$topmethcnt/topinfh$topcov
    botinfh <- subset(botinfh, botinfh$botcov > covthresh)
    botinfh$botmethper <- botinfh$botmethcnt/botinfh$botcov
    infh <- merge(topinfh, botinfh, by=c("target","chrom", "start", "end"), all=TRUE)
    infh$target <- ifelse(infh$target=="PET_CpG","Placenta","Blood")
    infh <- infh[infh$chrom %in% autosomes,]
    infh <- infh[complete.cases(infh),]
    outfh <- paste0(samplename, ".covthresh" ,covthresh, ".test.topbot.tsv")
    write.table(infh, outfh, row.names=FALSE, quote=FALSE, sep="\t")
    plotdt <- infh
    #plotdt <- subset(infh, chrom=="1")
    print(unique(plotdt$chrom))
    
    outfig <- paste0(samplename, ".covthresh", covthresh, ".test.topbot.pdf")
    pdf(outfig, width=9, height=9)
    p1 <- ggplot(plotdt, aes(x=topmethper, y=botmethper, color=target)) + geom_point(size=0.5, alpha=0.5)+theme_bw()+ylim(-0.01,0.8)+xlim(-0.01,0.8)+scale_colour_manual(values = rainbow(2)) 
    print(p1)
    dev.off()
}


#PlotFeatureProfile(workdir, samplename)
PlotPer(workdir, samplename, covthresh)