library(data.table)
library(pROC)
library(reshape2)
library(ggplot2)

# Combine eRPM/FMF scores with PE scores
## Combine scores using the training dataset and evaluate performance with leave-one-out analysis

TrainingCombineFMF <- function(infh,FMFscore,trainrocfig, trainoutfh, traincombrocfig) {
    print(FMFscore)
    fhall <- fread(infh, header=TRUE, sep="\t", data.table=FALSE)
    fh <- subset(fhall, Train.Validation=="Training")
    pescoreperfm <- roc(response=fh$Case.Ctrl,predictor=fh$PE.score,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    fmfperfm <- roc(response=fh$Case.Ctrl,predictor=fh[,c(FMFscore)],ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    print("ok")
    #fmfperfm <- roc(response=fh$Case.Ctrl,predictor=fh$FMFscore,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    pdf(trainrocfig, height=8, width=8)
    plot(pescoreperfm,col="#045a8d",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.5,ci=TRUE,lwd=4,cex.lab=1.5,cex.axis=1.5)
    plot(fmfperfm,col="#a6bddb",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    legend("bottomright",legend=c("cfDNAme score","FMF risk score"),col=c("#045a8d","#a6bddb"),lwd=3,cex=1.5)
    dev.off()
    predt <- data.frame()
    fh$Pheno <- ifelse(fh$Case.Ctrl=="Ctrl",0,1)
    fh$log10FMFSCORE <- log10(fh[,c(FMFscore)])

    trainallglm <- glm(Pheno ~ PE.score+log10FMFSCORE+SeqFF, family=binomial(link='logit'), data=fh)

    n_train <- nrow(fh)
    for (i in 1:n_train) {
        traindt <- fh[-i,]
	testvars <- as.data.frame(fh[i,c("PE.score","log10FMFSCORE")])
	trainglm <- glm(Pheno ~ PE.score+log10FMFSCORE, family=binomial(link='logit'), data=traindt)

	testpred <- predict(trainglm,newdata=testvars,type="response")
	testpredclass <- ifelse(testpred > 0.5,"Case","Ctrl")
	predres <- data.frame(testpred,testpredclass)
	colnames(predres) <- c("CombinedScore","CombinedPred")

	predt <- rbind(predt, predres)
    }
    outdt <- cbind(fh[,c("SubjectID","Case.Ctrl","PE.score","PE.score.label","FMFRISK34","Center")],predt)
    write.table(outdt, trainoutfh, row.names=FALSE, quote=FALSE, sep="\t")

    combineperfm <- roc(response=outdt$Case.Ctrl,predictor=outdt$CombinedScore,ci=TRUE,levels=c("Ctrl","Case"),direction="<")

    pdf(traincombrocfig, height=8, width=8)
    plot(combineperfm,col="#1b7837",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.45,ci=TRUE,lwd=4,cex.lab=1.5,cex.axis=1.5)
    plot(pescoreperfm,col="#762a83",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    plot(fmfperfm,col="#9970ab",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.35,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    legend("bottomright",legend=c("cfDNAme score+FMF risk score","cfDNAme score","FMF risk score"),col=c("#1b7837","#762a83","#c2a5cf"),lwd=3,cex=1.5)
    dev.off()
}

## Use the combined model to perform prediction on the validation dataset

ValidCombineFMF <- function(infh,FMFscore,validoutfh, validrocfig, valid1rocfig, valid2rocfig, classoutfig) {
    fh <- fread(infh, header=TRUE, sep="\t", data.table=FALSE)
    fh$Pheno <- ifelse(fh$Case.Ctrl=="Ctrl",0,1)
    fh$log10FMFSCORE <- log10(fh[,c(FMFscore)])
    traindt <- subset(fh, Train.Validation=="Training")
    testdt <- subset(fh, Train.Validation=="Validation")

    trainglm <- glm(Pheno ~ PE.score+log10FMFSCORE, family=binomial(link='logit'), data=traindt)
    print(summary(trainglm))
    print("odds ratio")
    print(exp(coef(trainglm)))
    #sink()

    testvars <- testdt[,c("PE.score","log10FMFSCORE")]

    testpred <- predict(trainglm,newdata=testvars,type="response")
    testpredclass <- ifelse(testpred > 0.5,"Case","Ctrl")
    predt <- cbind(testdt[,c("SubjectID","Case.Ctrl","PE.score","PE.score.label",FMFscore,"FMFRISK34","log10FMFSCORE","Center")],testpred,testpredclass)
    write.table(predt,validoutfh,row.names=FALSE,quote=FALSE,sep="\t")

    validpespf <- roc(response=predt$Case.Ctrl,predictor=predt$PE.score,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    validlfmfpf <- roc(response=predt$Case.Ctrl,predictor=predt$log10FMFSCORE,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    validpeadlfmfpf <- roc(response=predt$Case.Ctrl,predictor=predt$testpred,ci=TRUE,levels=c("Ctrl","Case"),direction="<")

    print("validpeadlfmfpf vs. validpespf")
    print(roc.test(validpeadlfmfpf, validpespf, method="bootstrap",boot.n=10000,boot.stratified=FALSE,alternative="greater"))
    print("validlfmfpf vs. validpespf")
    print(roc.test(validlfmfpf,validpespf, method="bootstrap",boot.n=10000,boot.stratified=FALSE,alternative="greater"))
    print("validpeadlfmfpf vs. validlfmfpf")
    print(roc.test(validpeadlfmfpf, validlfmfpf, method="bootstrap",boot.n=10000,boot.stratified=FALSE,alternative="greater"))

    print("validpeadlfmf")
    validpeadlfmfcoords90sp <- coords(roc=validpeadlfmfpf, x=0.9,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"), transpose = FALSE)
    print(validpeadlfmfcoords90sp)
    validpeadlfmfcoords80sp <- coords(roc=validpeadlfmfpf, x=0.8,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"),transpose = FALSE)
    print(validpeadlfmfcoords80sp)
    validpeadlfmfcoordsAll <- coords(roc=validpeadlfmfpf, "all", ret=c("threshold","specificity","sensitivity","accuracy","ppv","npv"),transpose=FALSE)
    print(validpeadlfmfcoordsAll)

    print("validlfmf")
    validlfmfcoords90sp <- coords(roc=validlfmfpf, x=0.9,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"), transpose = FALSE)
    print(validlfmfcoords90sp)
    validlfmfcoords80sp <- coords(roc=validlfmfpf, x=0.8,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"),transpose = FALSE)
    print(validlfmfcoords80sp)
    validlfmfcoordsAll <- coords(roc=validlfmfpf, "all", ret=c("threshold","specificity","sensitivity","accuracy","ppv","npv"),transpose=FALSE)
    print(validlfmfcoordsAll)

    print("validpes")
    validpescoords90sp <- coords(roc=validpespf, x=0.9,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"), transpose = FALSE)
    print(validpescoords90sp)
    validpescoords80sp <- coords(roc=validpespf, x=0.8,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"),transpose = FALSE)
    print(validpescoords80sp)
    validpescoordsAll <- coords(roc=validpespf, "all", ret=c("threshold","specificity","sensitivity","accuracy","ppv","npv"),transpose=FALSE)
    print(validpescoordsAll)

    pdf(validrocfig, height=8, width=8)
    plot(validpeadlfmfpf,col="#1b7837",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.45,ci=TRUE,lwd=4,cex.lab=1.5,cex.axis=1.5)
    plot(validpespf,col="#762a83",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    plot(validlfmfpf,col="#9970ab",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.35,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    legend("bottomright",legend=c("cfDNAme score+FMF risk score","cfDNAme score","FMF risk score"),col=c("#1b7837","#762a83","#c2a5cf"),lwd=3,cex=1.5)
    dev.off()

    ZOLpredt <- subset(predt,Center=="ZOL")
    ZOLvalidpespf <- roc(response=ZOLpredt$Case.Ctrl,predictor=ZOLpredt$PE.score,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    ZOLvalidlfmfpf <- roc(response=ZOLpredt$Case.Ctrl,predictor=ZOLpredt$log10FMFSCORE,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    ZOLvalidpeadlfmfpf <- roc(response=ZOLpredt$Case.Ctrl,predictor=ZOLpredt$testpred,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    pdf(valid1rocfig, height=8, width=8)
    plot(ZOLvalidpeadlfmfpf,col="#1b7837",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.45,ci=TRUE,lwd=4,cex.lab=1.5,cex.axis=1.5)
    plot(ZOLvalidpespf,col="#762a83",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    plot(ZOLvalidlfmfpf,col="#9970ab",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.35,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    legend("bottomright",legend=c("cfDNAme score+FMF risk score","cfDNAme score","FMF risk score"),col=c("#1b7837","#762a83","#c2a5cf"),lwd=3,cex=1.5)
    dev.off()

    SJBpredt <- subset(predt,Center=="SJB")
    SJBvalidpespf <- roc(response=SJBpredt$Case.Ctrl,predictor=SJBpredt$PE.score,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    SJBvalidlfmfpf <- roc(response=SJBpredt$Case.Ctrl,predictor=SJBpredt$log10FMFSCORE,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    SJBvalidpeadlfmfpf <- roc(response=SJBpredt$Case.Ctrl,predictor=SJBpredt$testpred,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    pdf(valid2rocfig, height=8, width=8)
    plot(SJBvalidpeadlfmfpf,col="#1b7837",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.45,ci=TRUE,lwd=4,cex.lab=1.5,cex.axis=1.5)
    plot(SJBvalidpespf,col="#762a83",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    plot(SJBvalidlfmfpf,col="#9970ab",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.35,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    legend("bottomright",legend=c("cfDNAme score+FMF risk score","cfDNAme score","FMF risk score"),col=c("#1b7837","#762a83","#c2a5cf"),lwd=3,cex=1.5)
    dev.off()

    classdf <- predt[,c("SubjectID","Case.Ctrl","PE.score.label","FMFRISK34","testpredclass")]
    colnames(classdf) <- c("SubjectID","Disease","cfDNAme","FMF risk score","cfDNAmeFMFriskscore")
    classdf_long <- melt(classdf,id.vars=c("SubjectID"),variable.name="Type",value.name="Class")

    pdf(classoutfig,height=9,width=6.3)
    p1 <- ggplot(classdf_long,aes(x=Type,y=SubjectID,fill=Class))+geom_tile()+theme_bw()+theme(axis.text.x = element_text(angle = 45))+scale_fill_manual(values=c("#b2182b", "#2166ac","#d6604d","#4393c3"))
    print(p1)
    dev.off()
}

### Combine PE, ePRM/FMF, and seqFF for prediction

TrainingCombineSeqFF_FMF <- function(infh,trainrocfig, trainoutfh, traincombrocfig) {
    fhall <- fread(infh, header=TRUE, sep="\t", data.table=FALSE)
    fh <- subset(fhall, Train.Validation=="Training")
    pescoreperfm <- roc(response=fh$Case.Ctrl,predictor=fh$PE.score,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    fmf37perfm <- roc(response=fh$Case.Ctrl,predictor=fh$FMFSCORE37,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    predt <- data.frame()
    fh$Pheno <- ifelse(fh$Case.Ctrl=="Ctrl",0,1)
    fh$log10FMFSCORE37 <- log10(fh$FMFSCORE37)

    trainallglm <- glm(Pheno ~ PE.score+log10FMFSCORE37+SeqFF, family=binomial(link='logit'), data=fh)
    print(summary(trainallglm))

    n_train <- nrow(fh)
    for (i in 1:n_train) {
        traindt <- fh[-i,]
	testvars <- as.data.frame(fh[i,c("PE.score","log10FMFSCORE37","SeqFF")])
	trainglm <- glm(Pheno ~ PE.score+log10FMFSCORE37+SeqFF, family=binomial(link='logit'), data=traindt)

	testpred <- predict(trainglm,newdata=testvars,type="response")
	testpredclass <- ifelse(testpred > 0.5,"Case","Ctrl")
	predres <- data.frame(testpred,testpredclass)
	colnames(predres) <- c("CombinedScore","CombinedPred")

	predt <- rbind(predt, predres)
    }
    print(head(fh))
    outdt <- cbind(fh[,c("SubjectID","Case.Ctrl","PE.score","PE.score.label","FMFSCORE37","FMFRISK37","Center","SeqFF")],predt)
    write.table(outdt, trainoutfh, row.names=FALSE, quote=FALSE, sep="\t")

    combineperfm <- roc(response=outdt$Case.Ctrl,predictor=outdt$CombinedScore,ci=TRUE,levels=c("Ctrl","Case"),direction="<")

    pdf(traincombrocfig, height=8, width=8)
    plot(combineperfm,col="#1b7837",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.45,ci=TRUE,lwd=4,cex.lab=1.5,cex.axis=1.5)
    plot(pescoreperfm,col="#762a83",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    plot(fmf37perfm,col="#9970ab",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.35,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    legend("bottomright",legend=c("cfDNAme score+FMF risk score+SeqFF","cfDNAme score","FMF risk score"),col=c("#1b7837","#762a83","#c2a5cf"),lwd=3,cex=1.5)
    dev.off()
}

#TrainingCombineSeqFF_FMF(infh,trainrocfig, trainoutfh, traincombrocfig)

ValidCombineSeqFF_FMF <- function(infh, validoutfh, validrocfig, valid1rocfig, valid2rocfig) {
    fh <- fread(infh, header=TRUE, sep="\t", data.table=FALSE)
    fh$Pheno <- ifelse(fh$Case.Ctrl=="Ctrl",0,1)
    fh$log10FMFSCORE37 <- log10(fh$FMFSCORE37)
    traindt <- subset(fh, Train.Validation=="Training")
    testdt <- subset(fh, Train.Validation=="Validation")

    trainglm <- glm(Pheno ~ PE.score+log10FMFSCORE37+SeqFF, family=binomial(link='logit'), data=traindt)
    print(summary(trainglm))
    print("odds ratio")
    print(exp(coef(trainglm)))
    #sink()

    testvars <- testdt[,c("PE.score","log10FMFSCORE37","SeqFF")]

    testpred <- predict(trainglm,newdata=testvars,type="response")
    testpredclass <- ifelse(testpred > 0.5,"Case","Ctrl")
    predt <- cbind(testdt[,c("SubjectID","Case.Ctrl","PE.score","PE.score.label","FMFSCORE37","log10FMFSCORE37","FMFRISK37","Center","SeqFF")],testpred,testpredclass)
    write.table(predt,validoutfh,row.names=FALSE,quote=FALSE,sep="\t")

    validpespf <- roc(response=predt$Case.Ctrl,predictor=predt$PE.score,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    validlfmf37pf <- roc(response=predt$Case.Ctrl,predictor=predt$log10FMFSCORE37,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    validpeadSeqFFlfmf37pf <- roc(response=predt$Case.Ctrl,predictor=predt$testpred,ci=TRUE,levels=c("Ctrl","Case"),direction="<")

    pdf(validrocfig, height=8, width=8)
    plot(validpeadSeqFFlfmf37pf,col="#1b7837",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.45,ci=TRUE,lwd=4,cex.lab=1.5,cex.axis=1.5)
    plot(validpespf,col="#762a83",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    plot(validlfmf37pf,col="#9970ab",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.35,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    legend("bottomright",legend=c("cfDNAme score+FMF risk score+SeqFF","cfDNAme score","FMF risk score"),col=c("#1b7837","#762a83","#c2a5cf"),lwd=3,cex=1.5)
    dev.off()

    ZOLpredt <- subset(predt,Center=="ZOL")
    ZOLvalidpespf <- roc(response=ZOLpredt$Case.Ctrl,predictor=ZOLpredt$PE.score,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    ZOLvalidlfmf37pf <- roc(response=ZOLpredt$Case.Ctrl,predictor=ZOLpredt$log10FMFSCORE37,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    ZOLvalidpeadSeqFFlfmf37pf <- roc(response=ZOLpredt$Case.Ctrl,predictor=ZOLpredt$testpred,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    pdf(valid1rocfig, height=8, width=8)
    plot(ZOLvalidpeadSeqFFlfmf37pf,col="#1b7837",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.45,ci=TRUE,lwd=4,cex.lab=1.5,cex.axis=1.5)
    plot(ZOLvalidpespf,col="#762a83",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    plot(ZOLvalidlfmf37pf,col="#9970ab",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.35,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    legend("bottomright",legend=c("cfDNAme score+FMF risk score+SeqFF","cfDNAme score","FMF risk score"),col=c("#1b7837","#762a83","#c2a5cf"),lwd=3,cex=1.5)
    dev.off()

    SJBpredt <- subset(predt,Center=="SJB")
    SJBvalidpespf <- roc(response=SJBpredt$Case.Ctrl,predictor=SJBpredt$PE.score,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    SJBvalidlfmf37pf <- roc(response=SJBpredt$Case.Ctrl,predictor=SJBpredt$log10FMFSCORE37,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    SJBvalidpeadSeqFFlfmf37pf <- roc(response=SJBpredt$Case.Ctrl,predictor=SJBpredt$testpred,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    pdf(valid2rocfig, height=8, width=8)
    plot(SJBvalidpeadSeqFFlfmf37pf,col="#1b7837",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.45,ci=TRUE,lwd=4,cex.lab=1.5,cex.axis=1.5)
    plot(SJBvalidpespf,col="#762a83",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    plot(SJBvalidlfmf37pf,col="#9970ab",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.35,ci=TRUE,lwd=4,add=TRUE,cex.lab=1.5,cex.axis=1.5)
    legend("bottomright",legend=c("cfDNAme score+FMF risk score+SeqFF","cfDNAme score","FMF risk score"),col=c("#1b7837","#762a83","#c2a5cf"),lwd=3,cex=1.5)
    dev.off()
}

#ValidCombineSeqFF_FMF(infh, validoutfh, validrocfig, valid1rocfig, valid2rocfig)


obsolete_CombineFMF <- function(infh,allrocoutfig, validrocoutfig,modeloutfh,predoutprefix,classoutfig) {
    fh <- fread(infh, header=TRUE, sep="\t", data.table=FALSE)
    pescoreperfm <- roc(response=fh$Case.Ctrl,predictor=fh$PE.score,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    fmf34perfm <- roc(response=fh$Case.Ctrl,predictor=fh$FMFSCORE34,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    fmf37perfm <- roc(response=fh$Case.Ctrl,predictor=fh$FMFSCORE37,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    pdf(allrocoutfig, height=8, width=8)
    plot(pescoreperfm,col="#045a8d",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.5,ci=TRUE,lwd=4)
    plot(fmf34perfm,col="#3690c0",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.45,ci=TRUE,lwd=4,add=TRUE)
    plot(fmf37perfm,col="#a6bddb",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE)
    legend("bottomright",legend=c("PEscore","FMFscore34","FMFscore37"),col=c("#045a8d","#3690c0","#a6bddb"),lwd=3)
    dev.off()
    
    fh$Pheno <- ifelse(fh$Case.Ctrl=="Ctrl",0,1)
    fh$log10FMFSCORE34 <- log10(fh$FMFSCORE34)
    fh$log10FMFSCORE37 <- log10(fh$FMFSCORE37)
    traindt <- subset(fh, Train.Validation=="Training")
    testdt <- subset(fh, Train.Validation=="Validation")
    trainglm1 <- glm(Pheno ~ PE.score+FMFSCORE34, family=binomial(link='logit'), data=traindt)
    sink(modeloutfh)
    print(summary(trainglm1))
    print("odds ratio")
    print(exp(coef(trainglm1)))
    
    trainglm2 <- glm(Pheno ~ PE.score+log10FMFSCORE34, family=binomial(link='logit'), data=traindt)
    print(summary(trainglm2))
    print("odds ratio")
    print(exp(coef(trainglm2)))

    trainglm3 <- glm(Pheno ~ PE.score+FMFSCORE37, family=binomial(link='logit'), data=traindt)
    print(summary(trainglm3))
    print("odds ratio")
    print(exp(coef(trainglm3)))

    trainglm4 <- glm(Pheno ~ PE.score+log10FMFSCORE37, family=binomial(link='logit'), data=traindt)
    print(summary(trainglm4))
    print("odds ratio")
    print(exp(coef(trainglm4)))
    sink()

    testvars1 <- testdt[,c("PE.score","FMFSCORE34")]
    testvars2 <- testdt[,c("PE.score","log10FMFSCORE34")]
    testvars3 <- testdt[,c("PE.score","FMFSCORE37")]
    testvars4 <- testdt[,c("PE.score","log10FMFSCORE37")]

    testpred1 <- predict(trainglm1,newdata=testvars1,type="response")
    testpredclass1 <- ifelse(testpred1 > 0.5,"Case","Ctrl")
    predt1 <- cbind(testdt[,c("SubjectID","Case.Ctrl","PE.score","PE.score.label","FMFSCORE34","FMFRISK34","FMFSCORE37","FMFRISK37")],testpred1,testpredclass1)
    outpred1 <- paste0(predoutprefix,".PEscore.FMFscore34.pred.out.tsv")
    write.table(predt1,outpred1,row.names=FALSE,quote=FALSE,sep="\t")

    testpred2 <- predict(trainglm2,newdata=testvars2,type="response")
    testpredclass2 <- ifelse(testpred2 > 0.5,"Case","Ctrl")
    predt2 <- cbind(testdt[,c("SubjectID","Case.Ctrl","PE.score","PE.score.label","FMFSCORE34","FMFRISK34","FMFSCORE37","FMFRISK37")],testpred2,testpredclass2)
    outpred2 <-	paste0(predoutprefix,".PEscore.logFMFscore34.pred.out.tsv")
    write.table(predt2,outpred2,row.names=FALSE,quote=FALSE,sep="\t")

    testpred3 <- predict(trainglm3,newdata=testvars3,type="response")
    testpredclass3 <- ifelse(testpred3 > 0.5,"Case","Ctrl")
    predt3 <- cbind(testdt[,c("SubjectID","Case.Ctrl","PE.score","PE.score.label","FMFSCORE34","FMFRISK34","FMFSCORE37","FMFRISK37")],testpred3,testpredclass3)
    outpred3 <-	paste0(predoutprefix,".PEscore.FMFscore37.pred.out.tsv")
    write.table(predt3,outpred3,row.names=FALSE,quote=FALSE,sep="\t")

    testpred4 <- predict(trainglm4,newdata=testvars4,type="response")
    testpredclass4 <- ifelse(testpred4 > 0.5,"Case","Ctrl")
    predt4 <- cbind(testdt[,c("SubjectID","Case.Ctrl","PE.score","PE.score.label","FMFSCORE34","FMFRISK34","FMFSCORE37","FMFRISK37")],testpred4,testpredclass4)
    outpred4 <- paste0(predoutprefix,".PEscore.logFMFscore37.pred.out.tsv")
    write.table(predt4,outpred4,row.names=FALSE,quote=FALSE,sep="\t")


    validpespf <- roc(response=testdt$Case.Ctrl,predictor=testdt$PE.score,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    validlfmf34pf <- roc(response=testdt$Case.Ctrl,predictor=testdt$log10FMFSCORE34,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    validlfmf37pf <- roc(response=testdt$Case.Ctrl,predictor=testdt$log10FMFSCORE37,ci=TRUE,levels=c("Ctrl","Case"),direction=">")
    validpeadlfmf34pf <- roc(response=predt2$Case.Ctrl,predictor=predt2$testpred2,ci=TRUE,levels=c("Ctrl","Case"),direction="<")
    validpeadlfmf37pf <- roc(response=predt4$Case.Ctrl,predictor=predt4$testpred4,ci=TRUE,levels=c("Ctrl","Case"),direction="<")

    print("validpeadlfmf37pf vs. validpespf")
    print(roc.test(validpeadlfmf37pf, validpespf))
    print(roc.test(validpeadlfmf37pf, validpespf, method="bootstrap",boot.n=5000))
    print("validpeadlfmf34pf vs. validpespf")
    print(roc.test(validpeadlfmf34pf, validpespf))
    print(roc.test(validpeadlfmf34pf, validpespf, method="bootstrap",boot.n=5000))
    print("validpeadlfmf37pf vs. validlfmf37pf")
    print(roc.test(validpeadlfmf37pf, validlfmf37pf))
    print(roc.test(validpeadlfmf37pf, validlfmf37pf, method="bootstrap",boot.n=5000))
    print("validpeadlfmf34pf vs. validlfmf34pf")
    print(roc.test(validpeadlfmf34pf, validlfmf34pf))
    print(roc.test(validlfmf34pf, validpeadlfmf34pf))

    print("validpeadlfmf37")
    validpeadlfmf37coords90sp <- coords(roc=validpeadlfmf37pf, x=0.9,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"), transpose = FALSE)
    print(validpeadlfmf37coords90sp)
    validpeadlfmf37coords80sp <- coords(roc=validpeadlfmf37pf, x=0.8,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"),transpose = FALSE)
    print(validpeadlfmf37coords80sp)

    print("validpeadlfmf34")
    validpeadlfmf34coords90sp <- coords(roc=validpeadlfmf34pf, x=0.9,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"), transpose = FALSE)
    print(validpeadlfmf34coords90sp)
    validpeadlfmf34coords80sp <- coords(roc=validpeadlfmf34pf, x=0.8,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"),transpose = FALSE)
    print(validpeadlfmf34coords80sp)

    print("validlfmf37")
    validlfmf37coords90sp <- coords(roc=validlfmf37pf, x=0.9,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"), transpose = FALSE)
    print(validlfmf37coords90sp)
    validlfmf37coords80sp <- coords(roc=validlfmf37pf, x=0.8,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"),transpose = FALSE)
    print(validlfmf37coords80sp)

    print("validlfmf34")
    validlfmf34coords90sp <- coords(roc=validlfmf34pf, x=0.9,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"), transpose = FALSE)
    print(validlfmf34coords90sp)
    validlfmf34coords80sp <- coords(roc=validlfmf34pf, x=0.8,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"),transpose = FALSE)
    print(validlfmf34coords80sp)

    print("validpes")
    validpescoords90sp <- coords(roc=validpespf, x=0.9,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"), transpose = FALSE)
    print(validpescoords90sp)
    validpescoords80sp <- coords(roc=validpespf, x=0.8,input="specificity", ret=c("sensitivity","accuracy","ppv","npv","threshold"),transpose = FALSE)
    print(validpescoords80sp)

    pdf(validrocoutfig, height=8, width=8)
    plot(validpeadlfmf37pf,col="#1b7837",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.5,ci=TRUE,lwd=4)
    plot(validpeadlfmf34pf,col="#5aae61",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.45,ci=TRUE,lwd=4,add=TRUE)
    plot(validpespf,col="#762a83",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.4,ci=TRUE,lwd=4,add=TRUE)
    plot(validlfmf37pf,col="#9970ab",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.35,ci=TRUE,lwd=4,add=TRUE)
    plot(validlfmf34pf,col="#c2a5cf",print.auc=TRUE,print.auc.x=0.3,print.auc.y=0.3,ci=TRUE,lwd=4,add=TRUE)
    legend("bottomright",legend=c("PEscore+FMFscore37","PEscore+FMFscore34","PEscore","FMFscore37","FMFscore34"),col=c("#1b7837","#5aae61","#762a83","#9970ab","#c2a5cf"),lwd=2)
    dev.off()

    classdt <- merge(predt2,predt4,by=c("SubjectID","Case.Ctrl","PE.score","PE.score.label","FMFSCORE34","FMFRISK34","FMFSCORE37","FMFRISK37"))
    print(head(classdt))
    classdf <- classdt[,c("SubjectID","Case.Ctrl","PE.score.label","FMFRISK34","FMFRISK37","testpredclass2","testpredclass4")]
    colnames(classdf) <- c("SubjectID","Disease","PEscoreMeth","FMFRISK34","FMFRISK37","MethPlus34","MethPlus37")
    classdf_long <- melt(classdf,id.vars=c("SubjectID"),variable.name="Type",value.name="Class")
    
    pdf(classoutfig,height=9,width=6.3)
    p1 <- ggplot(classdf_long,aes(x=Type,y=SubjectID,fill=Class))+geom_tile()+theme_bw()+theme(axis.text.x = element_text(angle = 45))+scale_fill_manual(values=c("#b2182b", "#2166ac","#d6604d","#4393c3"))
    print(p1)
    dev.off()
}

#CombineFMF(infh,allrocoutfig, validrocoutfig,modeloutfh,predoutprefix,classoutfig)
