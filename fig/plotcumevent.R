library(data.table)
library(survival)
library(survminer)

args <- commandArgs(TRUE)
infh <- args[1]
outfigprefix <- args[2]

plotcumevent <- function(infh,outfigprefix) {
    fh <- fread(infh,header=TRUE,sep="\t",data.table=FALSE)
    #fh$PEpred <- ifelse(fh$PEscorelabel=="Ctrl",0,1)
    fh$PEpred <- factor(fh$PEscorelabel,levels=c("Ctrl","Case"))
    fh$Metapred <- ifelse(fh$CombinedScore < 0.5,"Ctrl","Case")
    fh$Metapred <- factor(fh$Metapred,levels=c("Ctrl","Case"))
    fh$GADwks <- as.numeric(fh$GADdays/7)
    #fh$CaseCtrl <- factor(fh$CaseCtrl,level=c("Ctrl","Case"))
    fh$CaseCtrl <- ifelse(fh$CaseCtrl=="Ctrl",0,1)
    trainfh <- subset(fh,TrainValidation=="Training")
    validfh <- subset(fh,TrainValidation=="Validation")

    cfDNAmeTrainfit <- survfit(Surv(GADwks, CaseCtrl)~PEpred, data=trainfh)
    metaTrainfit <- survfit(Surv(GADwks, CaseCtrl)~Metapred, data=trainfh)
    trainfit <- list(cfDNAme=cfDNAmeTrainfit,cfDNAme_FMF=metaTrainfit)
    
    cfDNAmeValidfit <- survfit(Surv(GADwks, CaseCtrl)~PEpred, data=validfh)
    metaValidfit <- survfit(Surv(GADwks, CaseCtrl)~Metapred, data=validfh)
    validfit <- list(cfDNAme=cfDNAmeValidfit,cfDNAme_FMF=metaValidfit)

    outfig0train <- paste0(outfigprefix,".cfme.GAD.train.pdf")
    pdf(outfig0train,width=7.8,height=5.5)
    p1 <- ggsurvplot(cfDNAmeTrainfit,data=trainfh,fun="event",size=1,palette=c("gray70","royalblue4"),legend.labs=c("cfDNAme - low risk of PE", "cfDNAme - high risk of PE"),conf.int=TRUE,pval=TRUE,pval.method=TRUE,pval.method.size=6,pval.size=6,pval.method.coord=c(25,0.7),pval.coord=c(25,0.65),font.x = c(14, "bold", "black"),font.y=c(14,"bold","black"),font.tickslab = c(12, "bold", "black"),xlim=c(24,42),break.time.by=4,xlab="Gestational age at delivery (weeks)",ylab="Cummulative event - risk of developing PE")
    p1$plot <- p1$plot+theme(legend.title=element_text(size=14),legend.text = element_text(size=14))
    print(p1,newpage=FALSE)
    dev.off()

    outfig0valid <- paste0(outfigprefix,".cfme.GAD.valid.pdf")
    pdf(outfig0valid,width=7.8,height=5.5)
    p2 <- ggsurvplot(cfDNAmeValidfit,data=validfh,fun="event",size=1,palette=c("gray70","royalblue4"),legend.labs=c("cfDNAme - low risk of PE", "cfDNAme - high risk of PE"),conf.int=TRUE,pval=TRUE,pval.method=TRUE,pval.method.size=6,pval.size=6,pval.method.coord=c(25,0.7),pval.coord=c(25,0.65),font.x = c(14, "bold", "black"),font.y=c(14,"bold","black"),font.tickslab = c(12, "bold", "black"),xlim=c(24,42),break.time.by=4,xlab="Gestational age at delivery (weeks)",ylab="Cummulative event - risk of developing PE")
    p2$plot <- p2$plot+theme(legend.title=element_text(size=14),legend.text = element_text(size=14))
    print(p2,newpage=FALSE)
    dev.off()

    outfig1train <- paste0(outfigprefix,".cfme_fmf.GAD.train.pdf")
    pdf(outfig1train,width=7.8,height=5.5)
    p3 <- ggsurvplot(trainfit,data=trainfh,combine=TRUE,fun="event",size=1.5,palette=c("gray70","royalblue4","gray80","royalblue3"),linetype=c(2,2,1,1),legend.labs=c("cfDNAme - low risk of PE", "cfDNAme - high risk of PE","cfDNAme + FMF - low risk of PE", "cfDNAme + FMF - high risk of PE"),conf.int=FALSE,pval=TRUE,pval.method=TRUE,pval.method.size=7,pval.size=7,pval.method.coord=c(25,0.7),pval.coord=c(25,0.65),font.x = c(14, "bold", "black"),font.y=c(14,"bold","black"),font.tickslab = c(12, "bold", "black"),xlim=c(24,42),break.time.by=4,xlab="Gestational age at delivery (weeks)",ylab="Cummulative event - risk of developing PE")
    p3$plot <- p3$plot+theme(legend.title=element_text(size=14),legend.text = element_text(size=14))+guides(color=guide_legend(nrow=2))
    print(p3,newpage=FALSE)
    dev.off()

    outfig1valid <- paste0(outfigprefix,".cfme_fmf.GAD.valid.pdf")
    pdf(outfig1valid,width=7.8,height=5.5)
    p4 <- ggsurvplot(validfit,data=validfh,combine=TRUE,fun="event",size=1.5,palette=c("gray70","royalblue4","gray80","royalblue3"),linetype=c(2,2,1,1),legend.labs=c("cfDNAme - low risk of PE", "cfDNAme - high risk of PE","cfDNAme + FMF - low risk of PE", "cfDNAme + FMF - high risk of PE"),conf.int=FALSE,pval=TRUE,pval.method=TRUE,pval.method.size=7,pval.size=7,pval.method.coord=c(25,0.7),pval.coord=c(25,0.65),font.x = c(14, "bold", "black"),font.y=c(14,"bold","black"),font.tickslab = c(12, "bold", "black"),xlim=c(24,42),break.time.by=4,xlab="Gestational age at delivery (weeks)",ylab="Cummulative event - risk of developing PE")
    p4$plot <- p4$plot+theme(legend.title=element_text(size=14),legend.text = element_text(size=14))+guides(color=guide_legend(nrow=2))
    print(p4,newpage=FALSE)
    dev.off()
}

plotcumevent(infh,outfigprefix)
