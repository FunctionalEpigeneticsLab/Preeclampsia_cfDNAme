library(data.table)
library(VennDiagram)
library(eulerr)
library(plyr)
library(UpSetR)

args <- commandArgs(TRUE)
sumfh <- args[1]
outfigprefix <- args[2]

display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

plotvenn <- function(sumfh,outfigprefix) {
    fh <- fread(sumfh,header=FALSE,sep="\t",data.table=FALSE)
    colnames(fh) <- c("Index","Material","Type")
    FreshTissue_hypo <- fh[fh$Material=="FreshTissue" & fh$Type=="Hypomethylated",]$Index
    FFPE_hypo <- fh[fh$Material=="FFPE" & fh$Type=="Hypomethylated",]$Index
    Blood_hypo <- fh[fh$Material=="Blood" & fh$Type=="Hypomethylated",]$Index
    cfDNA_1stTrimester_hypo <- fh[fh$Material=="cfDNA_1stTrimester" & fh$Type=="Hypomethylated",]$Index
    cfDNA_atDiagnosis_hypo <- fh[fh$Material=="cfDNA_atDiagnosis" & fh$Type=="Hypomethylated",]$Index

    FreshTissue_hyper <- fh[fh$Material=="FreshTissue" & fh$Type=="Hypermethylated",]$Index
    FFPE_hyper <- fh[fh$Material=="FFPE" & fh$Type=="Hypermethylated",]$Index
    Blood_hyper <- fh[fh$Material=="Blood" & fh$Type=="Hypermethylated",]$Index
    cfDNA_1stTrimester_hyper <- fh[fh$Material=="cfDNA_1stTrimester" & fh$Type=="Hypermethylated",]$Index
    cfDNA_atDiagnosis_hyper <- fh[fh$Material=="cfDNA_atDiagnosis" & fh$Type=="Hypermethylated",]$Index
    #vlists <- list(FFPE,Blood,cfDNA_atDiagnosis,cfDNA_1stTrimester,FreshTissue)
    #pdf(outfig1,height=7,width=7)
    #p1 <- display_venn(
    #	    vlists,
    #	    category.names = c("FFPE","Blood","cfDNA at diagnosis","cfDNA 1st trimester","Fresh Tissue"),
	    # Circles
    #	    lwd = 2,
    #	    lty = 'blank',
    #	    fill = c("#E69F00", "#56B4E9", "#009E73", "#bcbddc","#999999"),
	    # Numbers
    #	    cex = .9,
    #	    fontface = "italic",
	    # Set names
    #	    cat.cex = 1,
    #	    cat.fontface = "bold",
    #	    cat.default.pos = "outer",
    #	    cat.dist = c(0.055, 0.065, 0.055, 0.035, 0.08)
    #)
    #print(p1)
    #dev.off()

    vlists1 <- list(FreshTissue_hypo,Blood_hypo,cfDNA_atDiagnosis_hypo)
    outfig1 <- paste0(outfigprefix,".hypo.3way.overlap.venn.pdf")
    pdf(outfig1,height=7,width=7)
    p1 <- display_venn(
            vlists1,
            category.names = c("Placenta","Blood","cfDNA at diagnosis"),
            # Circles
            lwd = 2,
            lty = 'blank',
            fill = c("#a6cee3","#bdbdbd","#b2df8a"),
            # Numbers
            cex = 1.5,
            fontface = "italic",
            # Set names
            cat.cex = 1.5,
            cat.fontface = "bold",
            cat.default.pos = "outer",
            cat.dist = c(0.055, 0.055, 0.055)
    )
    print(p1)
    dev.off()
 
    vlists2 <- list(FreshTissue_hyper,Blood_hyper,cfDNA_atDiagnosis_hyper)
    outfig2 <- paste0(outfigprefix,".hyper.3way.overlap.venn.pdf")
    pdf(outfig2,height=7,width=7)
    p2 <- display_venn(
            vlists2,
            category.names = c("Placenta","Blood","cfDNA at diagnosis"),
            # Circles
            lwd = 2,
            lty = 'blank',
            fill = c("#a6cee3","#bdbdbd","#b2df8a"),
            # Numbers
            cex = 1.5,
            fontface = "italic",
            # Set names
            cat.cex = 1.5,
            cat.fontface = "bold",
            cat.default.pos = "outer",
            cat.dist = c(0.055, 0.055, 0.055)
    )
    print(p2)
    dev.off()
}

plotvenn(sumfh,outfigprefix)

plotupset <- function(sumfh,outfigprefix) {
    fh <- fread(sumfh,header=FALSE,sep="\t",data.table=FALSE)
    colnames(fh) <- c("Index","Material","Type")
    MAY1_hypo <- fh[fh$Material=="FreshTissue" & fh$Type=="Hypomethylated",]$Index
    MAY1_hyper <- fh[fh$Material=="FreshTissue" & fh$Type=="Hypermethylated",]$Index
    MAY2_hypo <- fh[fh$Material=="FFPE" & fh$Type=="Hypomethylated",]$Index
    MAY2_hyper <- fh[fh$Material=="FFPE" & fh$Type=="Hypermethylated",]$Index
    MAY3_hypo <- fh[fh$Material=="Blood" & fh$Type=="Hypomethylated",]$Index
    MAY3_hyper <- fh[fh$Material=="Blood" & fh$Type=="Hypermethylated",]$Index
    MAY4_hypo <- fh[fh$Material=="cfDNA_atDiagnosis" & fh$Type=="Hypomethylated",]$Index
    MAY4_hyper <- fh[fh$Material=="cfDNA_atDiagnosis" & fh$Type=="Hypermethylated",]$Index
    MAY5_hypo <- fh[fh$Material=="cfDNA_1stTrimester" & fh$Type=="Hypomethylated",]$Index
    MAY5_hyper <- fh[fh$Material=="cfDNA_1stTrimester" & fh$Type=="Hypermethylated",]$Index
    
    vlists1 <- list(`Placenta`=MAY1_hypo,Blood=MAY3_hypo,`cfDNA at diagnosis`=MAY4_hypo,`cfDNA 1st trimester`=MAY5_hypo)
    vlists2 <- list(`Placenta`=MAY1_hyper,Blood=MAY3_hyper,`cfDNA at diagnosis`=MAY4_hyper,`cfDNA 1st trimester`=MAY5_hyper)
    
    outfig1 <- paste0(outfigprefix,".4way.hypo.overlap.upsetplot.pdf")
    pdf(outfig1,height=6.3,width=8.9)
    p1 <- upset(fromList(vlists1),sets=c("Placenta","Blood","cfDNA at diagnosis","cfDNA 1st trimester"),keep.order=TRUE,point.size=2,line.size=0.8,mainbar.y.label="Overlap regions",sets.x.label="Regions per material",text.scale = c(1.5, 1.5, 1.5, 1.3, 1.5, 1),main.bar.color="steelblue",empty.intersections="on")
    print(p1,newpage=FALSE)
    dev.off()

    outfig2 <- paste0(outfigprefix,".4way.hyper.overlap.upsetplot.pdf")
    pdf(outfig2,height=6.3,width=8.9)
    p2 <- upset(fromList(vlists2),sets=c("Placenta","Blood","cfDNA at diagnosis","cfDNA 1st trimester"),keep.order=TRUE,point.size=2,line.size=0.8,mainbar.y.label="Overlap regions",sets.x.label="Regions per material",text.scale = c(1.5, 1.5, 1.5, 1.3, 1.5, 1),main.bar.color="steelblue",empty.intersections="on")
    print(p2,newpage=FALSE)
    dev.off()    

    vlists3 <- list(`Placenta`=MAY1_hypo, Blood=MAY3_hypo, `cfDNA at diagnosis`=MAY4_hypo)
    vlists4 <- list(`Placenta`=MAY1_hyper, Blood=MAY3_hyper, `cfDNA at diagnosis`=MAY4_hyper)

    outfig3 <- paste0(outfigprefix,".3way.hypo.Fresh.blood.cfDNA2.overlap.upsetplot.pdf")
    pdf(outfig3,height=5.6,width=8)
    p3 <- upset(fromList(vlists3),sets=c("Placenta","Blood","cfDNA at diagnosis"),keep.order=TRUE,point.size=2.5,line.size=1,mainbar.y.label="Overlap regions",sets.x.label="Regions per material",text.scale = c(1.8, 1.8, 1.8, 1.5, 1.8, 1.3),main.bar.color="steelblue",empty.intersections="on")
    print(p3,newpage=FALSE)
    dev.off()

    outfig4 <- paste0(outfigprefix,".3way.hyper.Fresh.blood.cfDNA2.overlap.upsetplot.pdf")
    pdf(outfig4,height=5.6,width=8)
    p4 <- upset(fromList(vlists4),sets=c("Placenta","Blood","cfDNA at diagnosis"),keep.order=TRUE,point.size=2.5,line.size=1,mainbar.y.label="Overlap regions",sets.x.label="Regions per material",text.scale = c(1.8, 1.8, 1.8, 1.5, 1.8, 1.3),main.bar.color="steelblue",empty.intersections="on")
    print(p4,newpage=FALSE)
    dev.off()

    vlists5 <- list(`Placenta hypo`=MAY1_hypo, `Blood hypo`=MAY3_hypo, `cfDNA at diagnosis hypo`=MAY4_hypo,`Placenta hyper`=MAY1_hyper, `Blood hyper`=MAY3_hyper, `cfDNA at diagnosis hyper`=MAY4_hyper)
    outfig5 <- paste0(outfigprefix,".3way.hypohyper.Fresh.blood.cfDNA2.overlap.upsetplot.pdf")
    pdf(outfig5,height=6.3,width=11.8)
    p5 <- upset(fromList(vlists5),sets=c("Placenta hypo","Blood hypo","cfDNA at diagnosis hypo","Placenta hyper","Blood hyper","cfDNA at diagnosis hyper"),keep.order=TRUE,point.size=2.5,line.size=1,mainbar.y.label="Overlap regions",sets.x.label="Regions per material",text.scale = c(1.8, 1.8, 1.8, 1.5, 1.8, 1.3),main.bar.color="steelblue",empty.intersections="on")
    print(p5,newpage=FALSE)
    dev.off()
}

plotupset(sumfh,outfigprefix)

ploteuler <- function(sumfh,outfigp3refix) {
    fh <- fread(sumfh,header=FALSE,sep="\t",data.table=FALSE)
    colnames(fh) <- c("Index","Material","Type")
    FreshTissue_hypo <- rep(TRUE,length(fh[fh$Material=="FreshTissue" & fh$Type=="Hypomethylated",]$Index))
    names(FreshTissue_hypo) <- fh[fh$Material=="FreshTissue" & fh$Type=="Hypomethylated",]$Index
    FreshTissue_hyper <- rep(TRUE,length(fh[fh$Material=="FreshTissue" & fh$Type=="Hypermethylated",]$Index))
    names(FreshTissue_hyper) <- fh[fh$Material=="FreshTissue" & fh$Type=="Hypermethylated",]$Index

    FFPE_hypo <- rep(TRUE,length(fh[fh$Material=="FFPE" & fh$Type=="Hypomethylated",]$Index))
    names(FFPE_hypo) <- fh[fh$Material=="FFPE" & fh$Type=="Hypomethylated",]$Index
    FFPE_hyper <- rep(TRUE,length(fh[fh$Material=="FFPE" & fh$Type=="Hypermethylated",]$Index))
    names(FFPE_hyper) <- fh[fh$Material=="FFPE" & fh$Type=="Hypermethylated",]$Index

    Blood_hypo <- rep(TRUE,length(fh[fh$Material=="Blood" & fh$Type=="Hypomethylated",]$Index))
    names(Blood_hypo) <- fh[fh$Material=="Blood" & fh$Type=="Hypomethylated",]$Index
    Blood_hyper <- rep(TRUE,length(fh[fh$Material=="Blood" & fh$Type=="Hypermethylated",]$Index))
    names(Blood_hyper) <- fh[fh$Material=="Blood" & fh$Type=="Hypermethylated",]$Index

    cfDNA_1stTrimester_hypo <- rep(TRUE,length(fh[fh$Material=="cfDNA_1stTrimester" & fh$Type=="Hypomethylated",]$Index))
    names(cfDNA_1stTrimester_hypo) <- fh[fh$Material=="cfDNA_1stTrimester" & fh$Type=="Hypomethylated",]$Index
    cfDNA_1stTrimester_hyper <- rep(TRUE,length(fh[fh$Material=="cfDNA_1stTrimester" & fh$Type=="Hypermethylated",]$Index))
    names(cfDNA_1stTrimester_hyper) <- fh[fh$Material=="cfDNA_1stTrimester" & fh$Type=="Hypermethylated",]$Index

    cfDNA_atDiagnosis_hypo <- rep(TRUE,length(fh[fh$Material=="cfDNA_atDiagnosis" & fh$Type=="Hypomethylated",]$Index))
    names(cfDNA_atDiagnosis_hypo) <- fh[fh$Material=="cfDNA_atDiagnosis" & fh$Type=="Hypomethylated",]$Index
    cfDNA_atDiagnosis_hyper <- rep(TRUE,length(fh[fh$Material=="cfDNA_atDiagnosis" & fh$Type=="Hypermethylated",]$Index))
    names(cfDNA_atDiagnosis_hyper) <- fh[fh$Material=="cfDNA_atDiagnosis" & fh$Type=="Hypermethylated",]$Index

    vmat1 <- rbind.fill.matrix(t(FreshTissue_hypo), t(Blood_hypo))
    vmat1[is.na(vmat1)] <- FALSE
    eulmat1 <- t(vmat1)
    colnames(eulmat1) <- c("Placenta","Blood")
    eufit1 <- eulerr::euler(eulmat1,shape="ellipse")

    vmat2 <- rbind.fill.matrix(t(FreshTissue_hyper), t(Blood_hyper))
    vmat2[is.na(vmat2)] <- FALSE
    eulmat2 <- t(vmat2)
    colnames(eulmat2) <- c("Placenta","Blood")
    eufit2 <- eulerr::euler(eulmat2,shape="ellipse")

    outfig1 <- paste0(outfigprefix,".2way_fresh_blood.hypo.overlap.euler.pdf")
    pdf(outfig1,height=5.6,width=5.6)
    p1 <- plot(eufit1, fills = c("#a6cee3","#bdbdbd"), quantities = TRUE)
    print(p1)
    dev.off()

    outfig2 <- paste0(outfigprefix,".2way_fresh_blood.hyper.overlap.euler.pdf")
    pdf(outfig2,height=5.6,width=5.6)
    p2 <- plot(eufit2, fills = c("#a6cee3","#bdbdbd"), quantities = TRUE)
    print(p2)
    dev.off()

    vmat3 <- rbind.fill.matrix(t(cfDNA_atDiagnosis_hypo),t(cfDNA_1stTrimester_hypo))
    vmat3[is.na(vmat3)] <- FALSE
    eulmat3 <- t(vmat3)
    colnames(eulmat3) <- c("cfDNA at diagnosis","cfDNA 1st trimester")
    eufit3 <- eulerr::euler(eulmat3,shape="ellipse")
    outfig3 <- paste0(outfigprefix,".2way_cfDNA.hypo.overlap.euler.pdf")
    pdf(outfig3,height=5.6,width=5.6)
    p3 <- plot(eufit3,fills=c("#b2df8a", "#009E73"),quantities=TRUE)
    print(p3)
    dev.off()

    vmat4 <- rbind.fill.matrix(t(cfDNA_atDiagnosis_hyper),t(cfDNA_1stTrimester_hyper))
    vmat4[is.na(vmat4)] <- FALSE
    eulmat4 <- t(vmat4)
    colnames(eulmat4) <- c("cfDNA at diagnosis","cfDNA 1st trimester")
    eufit4 <- eulerr::euler(eulmat4,shape="ellipse")
    outfig4 <- paste0(outfigprefix,".2way_cfDNA.hyper.overlap.euler.pdf")
    pdf(outfig4,height=5.6,width=5.6)
    p4 <- plot(eufit4,fills=c("#b2df8a", "#009E73"),quantities=TRUE)
    print(p4)
    dev.off()

    vmat5 <- rbind.fill.matrix(t(FreshTissue_hypo),t(cfDNA_atDiagnosis_hypo))
    vmat5[is.na(vmat5)] <- FALSE
    eulmat5 <- t(vmat5)
    colnames(eulmat5) <- c("Placenta","cfDNA at diagnosis")
    eufit5 <- eulerr::euler(eulmat5,shape="ellipse")
    outfig5 <- paste0(outfigprefix,".2way_fresh_cfd.hypo.overlap.euler.pdf")
    pdf(outfig5,height=5.6,width=5.6)
    p5 <- plot(eufit5,fills=c("#a6cee3","#b2df8a"),quantities=TRUE)
    print(p5)
    dev.off()

    vmat6 <- rbind.fill.matrix(t(FreshTissue_hyper),t(cfDNA_atDiagnosis_hyper))
    vmat6[is.na(vmat6)] <- FALSE
    eulmat6 <- t(vmat6)
    colnames(eulmat6) <- c("Placenta","cfDNA at diagnosis")
    eufit6 <- eulerr::euler(eulmat6,shape="ellipse")
    outfig6 <- paste0(outfigprefix,".2way_fresh_cfd.hyper.overlap.euler.pdf")
    pdf(outfig6,height=5.6,width=5.6)
    p6 <- plot(eufit6,fills=c("#a6cee3","#b2df8a"),quantities=TRUE)
    print(p6)
    dev.off()

    vmat7 <- rbind.fill.matrix(t(FreshTissue_hypo), t(Blood_hypo), t(cfDNA_atDiagnosis_hypo))
    vmat7[is.na(vmat7)] <- FALSE
    eulmat7 <- t(vmat7)
    colnames(eulmat7) <- c("Placenta", "Blood", "cfDNA at diagnosis")
    eufit7 <- eulerr::euler(eulmat7,shape="ellipse")
    outfig7 <- paste0(outfigprefix,".3way_fresh_blood_cfd.hypo.overlap.euler.pdf")
    pdf(outfig7,height=5.6,width=5.6)
    p7 <- plot(eufit7, fills = c("#a6cee3","#bdbdbd","#b2df8a"), quantities = TRUE)
    print(p7)
    dev.off()

    vmat8 <- rbind.fill.matrix(t(FreshTissue_hyper), t(Blood_hyper), t(cfDNA_atDiagnosis_hyper))
    vmat8[is.na(vmat8)] <- FALSE
    eulmat8 <- t(vmat8)
    colnames(eulmat8) <- c("Placenta", "Blood", "cfDNA at diagnosis")
    eufit8 <- eulerr::euler(eulmat8,shape="ellipse")
    outfig8 <- paste0(outfigprefix,".3way_fresh_blood_cfd.hyper.overlap.euler.pdf")
    pdf(outfig8,height=5.6,width=5.6)
    p8 <- plot(eufit8, fills = c("#a6cee3","#bdbdbd","#b2df8a"), quantities = TRUE)
    print(p8)
    dev.off()
}

ploteuler(sumfh,outfigprefix)
