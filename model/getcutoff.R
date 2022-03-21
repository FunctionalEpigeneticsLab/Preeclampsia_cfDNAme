library(data.table)
library(pROC)

args <- commandArgs(TRUE)
infh <- args[1]
outfh <- args[2]
outfig <- args[3]

getcvprob <- function(infh, outfh, outfig) {
    fh <- fread(infh, header=FALSE, sep="\t", data.table=FALSE)
    colnames(fh) <- c("SubjectID","Phenotype","defaultcut","Prob")
    glmperfm <- roc(response=fh$Phenotype, predictor=fh$Prob, ci=TRUE,levels=c("Ctrl","Case"), direction="<")
    repauc <- paste0("AUC: ", glmperfm$auc)   
    repci <- paste("CI: ", glmperfm$ci)
    youdenbest <- coords(glmperfm, "best", ret="threshold", transpose = FALSE, best.method="youden")
    topleftbest <- coords(glmperfm, "best", ret="threshold", transpose = FALSE, best.method="closest.topleft")
    print(paste0("Youden best cutoff: ", youdenbest))
    print(paste0("Topleft best cutoff: ", topleftbest))
    ssupper <- 1.0
    sslower <- 0.5
    spupper <- 1.0
    splower <- 0.5
    #glmcoords <- coords(roc=glmperfm, x = "all", transpose = FALSE)
    #coordacc <- glmcoords[(glmcoords$specificity >= splower & glmcoords$specificity <= spupper & glmcoords$sensitivity >= sslower & glmcoords$sensitivity <= ssupper),]
    glmcoords <- coords(roc=glmperfm, x="all", ret=c("threshold","specificity","sensitivity","accuracy","precision","recall"), transpose=FALSE)
    glmcoords$F1 <- 2*glmcoords$precision*glmcoords$recall/(glmcoords$precision+glmcoords$recall)
    coordacc <- glmcoords[(glmcoords$specificity >= splower & glmcoords$specificity <= spupper & glmcoords$sensitivity >= sslower & glmcoords$sensitivity <= ssupper),]
    coordacc <- coordacc[order(coordacc$accuracy),]
    write.table(coordacc, outfh, sep="\t", row.names=FALSE, quote=FALSE)
    print(coordacc)

    pdf(outfig, height=9, width=9)
    plot(glmperfm,print.auc=TRUE,print.auc.y=0.6,print.auc.x=0.3)
    dev.off()
}

getcvprob(infh, outfh, outfig)
