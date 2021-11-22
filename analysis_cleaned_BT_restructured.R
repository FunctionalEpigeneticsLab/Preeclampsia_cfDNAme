source('~/Downloads/data_process.R', chdir = TRUE)
setwd("~/data_processing/preeclampsia/Freeze Sept 2021/Archive/")
library(grDevices)
##########################
# load in all count table
##########################

count <- read.table('./11N1_forward_count.txt',header = T,sep = '\t')
colnames(count) <- sub('.cov.gz','',colnames(count))
count$Probe_id <- paste(count$Probe, count$Chromosome, count$Start, count$End, sep = '_')
count <- count[,-c(1:13)]
N1sample11_forward_count <- count
rm(count)
count <- read.table('./11N1_reverse_count.txt',header = T,sep = '\t')
colnames(count) <- sub('.cov.gz','',colnames(count))
count$Probe_id <- paste(count$Probe, count$Chromosome, count$Start, count$End, sep = '_')
count <- count[,-c(1:13)]
N1sample11_reverse_count <- count
rm(count)
count <- read.table('./16N1_forward_count.txt',header = T,sep = '\t')
colnames(count) <- sub('.cov.gz','',colnames(count))
count$Probe_id <- paste(count$Probe, count$Chromosome, count$Start, count$End, sep = '_')
count <- count[,-c(1:13)]
N1sample16_forward_count <- count
rm(count)
count <- read.table('./16N1_reverse_count.txt',header = T,sep = '\t')
colnames(count) <- sub('.cov.gz','',colnames(count))
count$Probe_id <- paste(count$Probe, count$Chromosome, count$Start, count$End, sep = '_')
count <- count[,-c(1:13)]
N1sample16_reverse_count <- count
rm(count)
count <- read.table('./48N1_forward_count.txt',header = T,sep = '\t')
colnames(count) <- sub('.cov.gz','',colnames(count))
count$Probe_id <- paste(count$Probe, count$Chromosome, count$Start, count$End, sep = '_')
count <- count[,-c(1:13)]
N1sample48_forward_count <- count
rm(count)

count <- read.table('./48N1_reverse_count.txt',header = T,sep = '\t')
colnames(count) <- sub('.cov.gz','',colnames(count))
count$Probe_id <- paste(count$Probe, count$Chromosome, count$Start, count$End, sep = '_')
count <- count[,-c(1:13)]
N1sample48_reverse_count <- count
rm(count)

count <- read.table('./16N1_second_forward_count.txt',header = T,sep = '\t')
colnames(count) <- sub('.cov.gz','',colnames(count))
count$Probe_id <- paste(count$Probe, count$Chromosome, count$Start, count$End, sep = '_')
count <- count[,-c(1:13)]
N1sample_second_forward_count <- count
rm(count)

count <- read.table('./16N1_second_reverse_count.txt',header = T,sep = '\t')
colnames(count) <- sub('.cov.gz','',colnames(count))
count$Probe_id <- paste(count$Probe, count$Chromosome, count$Start, count$End, sep = '_')
count <- count[,-c(1:13)]
N1sample_second_reverse_count <- count
rm(count)

count <- read.table('./F_N120210827.txt',header = T,sep = '\t')
colnames(count) <- sub('.cov.gz','',colnames(count))
count$Probe_id <- paste(count$Probe, count$Chromosome, count$Start, count$End, sep = '_')
count <- count[,-c(1:13)]
N1sample_large_forward_count <- count
rm(count)

count <- read.table('./R_N120210827.txt',header = T,sep = '\t')
colnames(count) <- sub('.cov.gz','',colnames(count))
count$Probe_id <- paste(count$Probe, count$Chromosome, count$Start, count$End, sep = '_')
count <- count[,-c(1:13)]
N1sample_large_reverse_count <- count
rm(count)


colnames(N1sample_large_forward_count)<-c(as.matrix(MyNames[match(colnames(N1sample_large_forward_count),MyNames[,1]),2]))
colnames(N1sample_large_forward_count)[94]<-"Probe_id"
colnames(N1sample_large_reverse_count)<-colnames(N1sample_large_forward_count)

forward_count <- cbind(
    N1sample11_forward_count[,-12],
    N1sample16_forward_count[,-17],
    N1sample48_forward_count[,-49],
    N1sample_second_forward_count[,-17], 
    N1sample_large_forward_count[,-94]
)
forward_count$Probe_id <- N1sample11_forward_count$Probe_id

reverse_count <- cbind(
    N1sample11_reverse_count[,-12],
    N1sample16_reverse_count[,-17],
    N1sample48_reverse_count[,-49],
    N1sample_second_reverse_count[,-17],
    N1sample_large_reverse_count[,-94]
)
reverse_count$Probe_id <- N1sample11_reverse_count$Probe_id

################
# remove count files of samples with poor BS-conversion
################
minimalConversion <- 0.985
QCData <- as.data.frame(read.table("BSseqQC.txt", header = TRUE, sep = "\t"))
PoorBSconversionSamples<-QCData$Sample.BT2[QCData$Bisulfite.Conv < minimalConversion]
PoorBSconversionSamples
reverse_count<-reverse_count[! names(reverse_count) %in% PoorBSconversionSamples]
forward_count <-forward_count[! names(forward_count) %in% PoorBSconversionSamples]


################
# merge count files of samples analysed more than once
################
# create matrix with data file info
DataFiles <- c(names(reverse_count[,-ncol(reverse_count)]))

#the last experiment has different file names, so convert to conventional ones
LastExpConversion<-as.matrix(read.csv(file = "Code15092021.csv"))
LastExpConversion<-LastExpConversion[LastExpConversion[,1]  %in%  DataFiles,]
rownames(LastExpConversion)<-LastExpConversion[,1]

DataFiles2<-DataFiles
DataFiles2[DataFiles %in% LastExpConversion[,1]]<-LastExpConversion[DataFiles[DataFiles %in% LastExpConversion[,1]],2]

DataFileMatrix<-as.data.frame(t(matrix(data = unlist(strsplit(DataFiles2, "_")), nrow = 2, ncol = length(DataFiles))))
DataFileMatrix[,2]<-sub("N1", "", DataFileMatrix[,2])
DataFileMatrix[,3]<-DataFiles

RepeatAnalysedSamples<-unique(DataFileMatrix[duplicated(DataFileMatrix[,2]),2])

for (RepSample in RepeatAnalysedSamples){
forward_count$temp<-rowSums(forward_count[,DataFileMatrix[,2] == RepSample])
names(forward_count)[names(forward_count) == "temp"]<-paste0("Repeated_",RepSample)
reverse_count $temp<-rowSums(reverse_count[,DataFileMatrix[,2] == RepSample])
names(reverse_count)[names(reverse_count) == "temp"]<-paste0("Repeated_",RepSample)
}
forward_count <-forward_count[, !colnames(forward_count) %in% c(DataFiles[DataFileMatrix[,2] %in% RepeatAnalysedSamples], "Probe_id")]
reverse_count <-reverse_count[, !colnames(reverse_count) %in% c(DataFiles[DataFileMatrix[,2] %in% RepeatAnalysedSamples], "Probe_id")]

reverse_count$Probe_id <- N1sample11_reverse_count$Probe_id
forward_count$Probe_id <- N1sample11_forward_count$Probe_id


DataFiles <- c(names(reverse_count[,-ncol(reverse_count)]))

#recreate the DataFileMatrix, but with repeat samples included
LastExpConversion<-as.matrix(read.csv(file = "Code15092021.csv"))
LastExpConversion<-LastExpConversion[LastExpConversion[,1]  %in%  DataFiles,]
rownames(LastExpConversion)<-LastExpConversion[,1]

DataFiles2<-DataFiles
DataFiles2[DataFiles %in% LastExpConversion[,1]]<-LastExpConversion[DataFiles[DataFiles %in% LastExpConversion[,1]],2]

DataFileMatrix<-as.data.frame(t(matrix(data = unlist(strsplit(DataFiles2, "_")), nrow = 2, ncol = length(DataFiles))))
DataFileMatrix[,2]<-sub("N1", "", DataFileMatrix[,2])
DataFileMatrix[,3]<-DataFiles



# replace(DataFiles, LastExpConversion[,1], LastExpConversion[,2])

# DataFileMatrix[,2]<-sub("N1", "", DataFileMatrix[,2])
# # DataFileMatrix[,2]<-sub("PE", "", DataFileMatrix[,2])
# DataFileMatrix[,3]<-DataFiles

################
## Remove samples without sufficient coverage, create new matrices
################

total_count <- forward_count[, !colnames(forward_count) %in% 'Probe_id'] + reverse_count[,!colnames(reverse_count) %in% 'Probe_id']
# MaxFraction = 25
# MinCount = 40
MaxFraction = 30
MinCount = 50

# # plot(colSums(total_count < MinCount) / nrow(total_count), colSums(total_count)/ nrow(total_count), xlab = paste0("% of probes which have total_count < ", MinCount), ylab = "average count per probe", log = "y", pch = 16)
# points(colSums(total_count < 30) / nrow(total_count), colSums(total_count)/ nrow(total_count), col = "green", pch = 16)
# points(colSums(total_count < 40) / nrow(total_count), colSums(total_count)/ nrow(total_count), col = "red", pch = 16)
# points(colSums(total_count < 70) / nrow(total_count), colSums(total_count)/ nrow(total_count), col = "blue", pch = 16)
# lines(c(MaxFraction, MaxFraction)/100, c(1,500))
# legend(0.7,500, c(30,40,50,70), fill = c("green",  "red","black", "blue"))


keep_samples <- DataFileMatrix[!(colSums(total_count < MinCount) / nrow(total_count)*100) > MaxFraction,3]
total_count<-total_count[names(reverse_count) %in% keep_samples]
reverse_count2 <-reverse_count[names(reverse_count) %in% keep_samples]
forward_count2 <-forward_count[names(forward_count) %in% keep_samples]
DataFileMatrix  <-DataFileMatrix[DataFileMatrix[,3] %in% keep_samples,]

colnames(forward_count2)<-DataFileMatrix[,2]
colnames(reverse_count2)<-DataFileMatrix[,2]
colnames(total_count)<-DataFileMatrix[,2]

forward_count2$Probe_id <- N1sample11_forward_count$Probe_id
reverse_count2$Probe_id <- N1sample11_forward_count$Probe_id


meth <- forward_count2[,!colnames(forward_count2) %in% 'Probe_id']/(forward_count2[,!colnames(forward_count2) %in% 'Probe_id'] + reverse_count2[,!colnames(reverse_count2) %in% 'Probe_id'])*100

######
# clean up data matrix: remove probes without sufficient measurements, replace NAs with Median
######

meth[total_count < 20] <- NA
# ProbesOI <- N1sample11_forward_count$Probe_id
ProbesOI <- N1sample11_forward_count$Probe_id[rowSums(is.na(meth)) < 0.1*ncol(total_count)]
meth <- meth[rowSums(is.na(meth)) < 0.1*ncol(total_count),]
forward_count2<-forward_count2[]

hist(rowSums(meth == 0, na.rm = T),breaks = 50)

MeanMeth<-colMeans(meth, na.rm = T)
MedianMeth<-apply(meth,2,median, na.rm = T)
MedianMethMatrix<-t(matrix(MedianMeth , ncol = nrow(meth), nrow = ncol(meth)))
tmp <- meth
tmp[is.na(tmp)] <- MedianMethMatrix[is.na(tmp)] 
rownames(tmp)<-ProbesOI
meth.noNA<-tmp #storing for further down the line



################
### load sample annotations
################

SampleData <- as.data.frame(read.table("characteristics20211109.txt", header = TRUE, sep = "\t"))
length(SampleData$StudyID)
## number of samples
length(DataFileMatrix[,2])
## number of samples with annotation
sum(DataFileMatrix[,2] %in% SampleData$StudyID)
## samples without annotation
DataFileMatrix[,2][!DataFileMatrix[,2] %in% SampleData$StudyID]


################
## Remove samples without annotation
################
NoAnnotationSamples<-DataFileMatrix[,2][!DataFileMatrix[,2] %in% SampleData$StudyID]
reverse_count2<-reverse_count2[! names(reverse_count2) %in% NoAnnotationSamples]
forward_count2 <-forward_count2[! names(forward_count2) %in% NoAnnotationSamples]
meth<-meth[! names(meth) %in% NoAnnotationSamples]
meth.noNA <-meth.noNA[! names(meth.noNA) %in% NoAnnotationSamples]

################
## Remove samples with excessively late NIPT1
################

#remove sample data without sequencing data 
SampleData<-SampleData[SampleData$StudyID %in% DataFileMatrix[,2],]
rownames(SampleData)<-SampleData$StudyID
SampleData<-SampleData[names(meth),]

table(SampleData$Capture.Date.N1)
hist(SampleData$FF)
SampleData$Capture.Date.N1<-as.factor(SampleData$Capture.Date.N1)


boxplot(SampleData$GA.NIPT.Days/7 ~ SampleData$Case.Ctrl)
maxWeekNIPT1<-15
LateNIPT1_Samples<-as.character(SampleData$StudyID[SampleData$GA.NIPT.Days/7> maxWeekNIPT1])
LateNIPT1_Samples
reverse_count2<-reverse_count2[! names(reverse_count2) %in% LateNIPT1_Samples]
forward_count2 <-forward_count2[! names(forward_count2) %in% LateNIPT1_Samples]
meth<-meth[! names(meth) %in% LateNIPT1_Samples]
meth.noNA <-meth.noNA[! names(meth.noNA) %in% LateNIPT1_Samples]
SampleData<-SampleData[! rownames(SampleData) %in% LateNIPT1_Samples,]
total_count<-total_count[,rownames(SampleData)]


################################################
#                                   ANALYSIS                                        #
################################################

################
## PCA plots
################
tmp<-meth.noNA
Meth.pca <- prcomp(t(tmp), scale = T)
summary(Meth.pca)

pdf("PCAplots.pdf")
plot(Meth.pca $x[,1:2], col = c("red", "black")[SampleData$Case.Ctrl], pch = 16)
text(Meth.pca $x[,1:2], rownames(Meth.pca $x), cex = 0.5, pos = 3, col = c("red", "black")[SampleData$Case.Ctrl])

plot(Meth.pca $x[,3:4], col = c("red", "black")[SampleData$Case.Ctrl], pch = 16)
text(Meth.pca $x[,3:4], rownames(Meth.pca $x), cex = 0.5, pos = 3, col = c("red", "black")[SampleData$Case.Ctrl])
plot(Meth.pca $x[,5:6], col = c("red", "black")[SampleData$Case.Ctrl], pch = 16)
text(Meth.pca $x[,5:6], rownames(Meth.pca $x), cex = 0.5, pos = 3, col = c("red", "black")[SampleData$Case.Ctrl])
plot(Meth.pca $x[,1:2], col = c("red", "black", "green", "grey70", "royalblue")[SampleData$Capture.Date.N1], pch = 16)
text(Meth.pca $x[,1:2], rownames(Meth.pca $x), cex = 0.5, pos = 3, col = c("red", "black", "green", "blue")[SampleData$Capture.Date.N1])
rbPal <- colorRampPalette(c('red',"yellow",'blue'))
plot(Meth.pca $x[,1:2], col = rbPal(100)[as.numeric(cut(SampleData$GA.NIPT.Days, breaks = 100))], pch = 16)
text(Meth.pca $x[,1:2], rownames(Meth.pca $x), cex = 0.5, pos = 3, col = "black")
plot(Meth.pca $x[,1:2], col = rbPal(100)[as.numeric(cut(colMeans(tmp), breaks = 100))], pch = 16)
cor(Meth.pca$x, y = colMeans(tmp))


SampleData$bisulfite.conv.NIPT1<-droplevels(SampleData$bisulfite.conv.NIPT1)
SampleData$FF<-as.numeric(as.character(SampleData$FF))
SampleData$bisulfite.conv.NIPT1<-as.numeric(levels(SampleData$bisulfite.conv.NIPT1))[SampleData$bisulfite.conv.NIPT1]

plot(Meth.pca $x[,1:2], col = rbPal(100)[cut(SampleData$bisulfite.conv.NIPT1, breaks = 100)], pch = 16)
text(Meth.pca $x[,1:2], rownames(Meth.pca $x), cex = 0.5, pos = 3, col = "black")
text(Meth.pca $x[,1:2], labels = round(SampleData$bisulfite.conv.NIPT1,3), cex = 0.5, pos = 1, col = "black")



plot(Meth.pca $x[,1:2], col = rbPal(100)[cut(colMeans(total_count<1), breaks = 100)], pch = 16)
par(mfrow=c(6,1), mar = c(3,4,1,0), mgp = c(1.8,0.6,0))
barplot(summary(Meth.pca)$importance[2,], beside = T, ylim = c(0,0.12), main = "variance explained per PC", las = 2)
barplot(cor(Meth.pca $x, y = colMeans(total_count<1))^2, beside = T, ylim = c(0,1), main = "fraction of regions without coverage", las = 2)
barplot(cor(Meth.pca $x, y = colMeans(meth, na.rm = T))^2, beside = T, ylim = c(0,1), main = "average methylation", las = 2)
barplot(cor(Meth.pca $x, y = SampleData$GA.NIPT.Days, use = "pairwise.complete")^2, beside = T, ylim = c(0,1), main = "gestational age (days) of NIPT" , las = 2)
barplot(cor(Meth.pca $x, y = c(1:2)[ifelse(SampleData$Case.Ctrl=="Ctrl",1,2)], use = "pairwise.complete")^2, beside = T, ylim = c(0,1), main = "case or control", las = 2)
barplot(cor(Meth.pca $x, y = SampleData$bisulfite.conv.NIPT1, use = "pairwise.complete")^2, beside = T, ylim = c(0,1), main = "bisulfite non-conversion estimate" , las = 2)

boxplot(colMeans(tmp)~ SampleData$Case.Ctrl, pch = 16)
t.test(colMeans(tmp)~ SampleData$Case.Ctrl, pch = 16)
dev.off()


probe_gene <- read.csv('probe_gene_annotation.csv',sep = ',',header = T)
#head(probe_gene)
probe_gene$Probe_id <- paste(probe_gene$Probe, probe_gene$Chromosome, probe_gene$Start, probe_gene$End, sep = '_')

################
# moderated t statistics
################

tmp <- as.data.frame(t(meth.noNA))
names(tmp)<-ProbesOI
tmp$condition1 <- factor(SampleData$Case.Ctrl, levels = c('Ctrl','Case'))
tmp$exp <- factor(SampleData$Capture.Date.N1)
tmp$GA_NIPT <- SampleData$GA.NIPT.Days
tmp$FF <- SampleData$FF
tmp$MeanMethylation <-colMeans(meth, na.rm = T)
tmp$coverageRate <- colMeans(total_count<1, na.rm = T)
design <- model.matrix(~ tmp$GA_NIPT)

# tmp <- meth.t.FF.filtered.B[,!colnames(meth.t.FF.filtered.B) %in% c('exp','condition2','FF','sample')]
# tmp$condition1 <- factor(tmp$condition1, levels = c('Ctrl','Case'))
# design <- model.matrix(~ tmp$condition1)
# design <- model.matrix(~ tmp$FF)

new_tmp <- t(tmp[,!colnames(tmp) %in% c('condition1', "exp", "MeanMethylation", "coverageRate", "GA_NIPT", "FF", "GAD.days")])
# new_tmp <- new_tmp[,tmp$condition1 == "Ctrl"]
# design <- model.matrix(~ tmp$GA_NIPT[tmp$condition1 == "Ctrl"]+tmp$MeanMethylation[tmp$condition1 == "Ctrl"]+tmp$coverageRate[tmp$condition1 == "Ctrl"])
design <- model.matrix(~ tmp$condition1+ tmp$MeanMethylation + tmp$coverageRate + tmp$FF)
# design <- model.matrix(~ tmp$GA_NIPT[tmp$condition1 == "Ctrl"])
fit <- lmFit(new_tmp, design)  # Column 1 contains row-names
fit2 <- eBayes(fit)
# fit3 <- eBayes(fit)

# cor(fit3$t[,2],fit2$t[,2], pch = 16, cex = 0.5)
sum(fit2$p.value[,2] < 0.05) # 
top200 <- rownames(head(fit2$p.value[order(fit2$p.value[,2], decreasing = F),],200))
# top200_GA_NIPT_associated <- rownames(head(fit2$p.value[order(fit2$p.value[,2], decreasing = F),],200))
# top.sig_GA_NIPT_associated <- rownames(fit2$p.value[which(fit2$p.value[,2] < 0.05),])
top100 <- rownames(head(fit2$p.value[order(fit2$p.value[,2], decreasing = F),],100))
top50 <- rownames(head(fit2$p.value[order(fit2$p.value[,2], decreasing = F),],50))
top.sig <- rownames(fit2$p.value[which(fit2$p.value[,2] < 0.05),])

MyCols  <- densCols(-log10(fit2$p.value[,2])~fit2$coefficients[,2], nbin=200, colramp = colorRampPalette(c("grey85", "black")))
plot(-log10(fit2$p.value[,2])~fit2$coefficients[,2], cex = 0.4, pch = 16, xlim = c(-5,5), col = MyCols)


par(mfrow=c(4,10), mar = c(4,4,0,0), mgp = c(1.8,0.6,0))
for (i in top200[1:40]){boxplot(tmp[,i] ~tmp$condition1, las = 2, ylab = "% methylation", xlab = NA)}
for (i in top50[1:40]){boxplot(tmp[,i] ~c(ceiling(tmp$GA_NIPT/7)), las = 2, ylab = "% methylation", xlab = NA)}



par(mfrow=c(4,10), mar = c(4,4,0,0), mgp = c(1.8,0.6,0))
for (i in top200_GA_NIPT_ass[1:40]){
	plot(tmp$GA_NIPT,tmp[,i], pch = 16, col = ifelse(tmp$condition1=="Ctrl", 2,1), cex = 1)
	abline(lm(tmp[tmp$condition1=="Ctrl",i]~tmp$GA_NIPT[tmp$condition1=="Ctrl"]), col = 2)
	abline(lm(tmp[tmp$condition1!="Ctrl",i]~tmp$GA_NIPT[tmp$condition1!="Ctrl"]), col = 1)
}



result.limma <- topTable(fit2, number = ncol(new_tmp))
nrow(subset(result.limma,subset =  P.Value < 0.05 & logFC < 0))
nrow(subset(result.limma,subset =  P.Value < 0.05 & logFC > 0))
nrow(subset(result.limma,subset =  adj.P.Val < 0.05 & logFC < 0))
nrow(subset(result.limma,subset =  adj.P.Val < 0.05 & logFC > 0))
ggplot(result.limma, aes(x=logFC, y=-log10(P.Value)))+geom_point()+theme_bw()+
geom_vline(xintercept=c(-1, 1), col="red",linetype="dashed") +
geom_hline(yintercept=-log10(0.05), col="red",linetype="dashed") 
result.limma <- cbind(result.limma, probe_gene[row.names(result.limma),c('Feature','ID')])




# order of samples has been arranged as Case then Control
tmp <- data.matrix(tmp[,top200])
#tmp <- scale(tmp) #?
tmpp <- data.matrix(tmp[,top200])
hclust.ave <- function(x) hclust(x, method="ward.D2")

col1 <- brewer.pal(11, "Spectral")
heatmap.2(scale(tmpp),scale = 'none', hclustfun = hclust.ave,trace = 'none', col = col1, cexCol = 0.4, RowSideColors=col1[as.numeric(factor(tmp$condition1))], margins=c(15,8)) #

train_ratio<-4/5
top_n<-200
alpha<-0.2

#####################
# Calculate AUCs
######################

tmp <- as.matrix(t(meth.noNA))
tmp<-tmp/rowMeans(tmp, na.rm = T)
# tmp<-tmp/FFmeth
# tmp<-tmp/apply(tmp, 1, median, na.rm= T)
# tmp<-tmp/as.numeric(SampleData$FF)
tmp <- as.data.frame(tmp)
names(tmp)<-ProbesOI
tmp$condition1 <- factor(SampleData$Case.Ctrl, levels = c('Case','Ctrl'))
tmp$exp <- factor(SampleData$Capture.Date.N1)
tmp$GA_NIPT <- SampleData$GA.NIPT.Days
tmp$GAD.days <- SampleData$GAD.days
tmp$FF <- as.numeric(SampleData$FF)
tmp$FFmeth <- FFmeth
tmp$MeanMethylation <-colMeans(meth, na.rm = T)
tmp$coverageRate <- colMeans(total_count<1, na.rm = T)

## remove experiment 200904
tmp <- tmp[tmp$exp != "200904", ]
tmp <- tmp[!is.na(tmp$GAD.days),]
EarlyGAD<-rownames(tmp)[tmp$GAD.days<217 & tmp$condition1=="Case"]
MyControls<-rownames(tmp)[tmp$GAD.days>250 & tmp$condition1=="Ctrl"]

paste0("we have ", length(EarlyGAD)," cases and ",length(MyControls)," controls")
tmp <- tmp[c(EarlyGAD, MyControls), ]
meth.t.FF.filtered <- tmp
AUC.collection = data.frame(matrix(ncol = 10, nrow = 1))
nonzero.collection = list()
nonzero.collection[[1]] <- list(); nonzero.collection[[2]] <- list(); nonzero.collection[[3]] <- list(); nonzero.collection[[4]] <- list(); nonzero.collection[[5]] <- list(); nonzero.collection[[6]] <- list(); nonzero.collection[[7]] <- list(); nonzero.collection[[8]] <- list(); nonzero.collection[[9]] <- list();nonzero.collection[[10]] <- list()
false.positive <- nonzero.collection
true.positive <- nonzero.collection

AUC_Reps = 250
PE_scores1<-matrix(data = NA, ncol = nrow(meth.t.FF.filtered), nrow= AUC_Reps, dimnames = list(1: AUC_Reps, rownames(meth.t.FF.filtered)))
PE_scores<-list(PE_scores1,PE_scores1,PE_scores1,PE_scores1,PE_scores1,PE_scores1,PE_scores1,PE_scores1,PE_scores1, PE_scores1)
for(i in seq(1:AUC_Reps)){
    tmp <- apply_model(meth.t.FF.filtered, 9/10, 1000, 0.2)
    AUC.collection[i,1] <- tmp[1]
    nonzero.collection[[1]] <- append(nonzero.collection[[1]],tmp[2])
    false.positive[[1]] <- append(false.positive[[1]], tmp[[3]]@x.values)
    true.positive[[1]] <- append(true.positive[[1]], tmp[[3]]@y.values)
	PE_scores[[1]][i,names(tmp[[4]])]<-tmp[[4]]

    tmp <- apply_model(meth.t.FF.filtered, 9/10, 800, 0.2)
    AUC.collection[i,2] <- tmp[1]
    nonzero.collection[[2]] <- append(nonzero.collection[[2]],tmp[2])
    false.positive[[2]] <- append(false.positive[[2]], tmp[[3]]@x.values)
    true.positive[[2]] <- append(true.positive[[2]], tmp[[3]]@y.values)
	PE_scores[[2]][i,names(tmp[[4]])]<-tmp[[4]]

    tmp <- apply_model(meth.t.FF.filtered, 9/10, 600, 0.2)
    AUC.collection[i,3] <- tmp[1]
    nonzero.collection[[3]] <- append(nonzero.collection[[3]],tmp[2])
    false.positive[[3]] <- append(false.positive[[3]], tmp[[3]]@x.values)
    true.positive[[3]] <- append(true.positive[[3]], tmp[[3]]@y.values)    
	PE_scores[[3]][i,names(tmp[[4]])]<-tmp[[4]]

    tmp <- apply_model(meth.t.FF.filtered, 9/10, 400, 0.15)
    AUC.collection[i,4] <- tmp[1]
    nonzero.collection[[4]] <- append(nonzero.collection[[4]],tmp[2])
    false.positive[[4]] <- append(false.positive[[4]], tmp[[3]]@x.values)
    true.positive[[4]] <- append(true.positive[[4]], tmp[[3]]@y.values)    
	PE_scores[[4]][i,names(tmp[[4]])]<-tmp[[4]]

    tmp <- apply_model(meth.t.FF.filtered, 9/10, 400, 0.3)
    AUC.collection[i,5] <- tmp[1]
    nonzero.collection[[5]] <- append(nonzero.collection[[5]],tmp[2])
    false.positive[[5]] <- append(false.positive[[5]], tmp[[3]]@x.values)
    true.positive[[5]] <- append(true.positive[[5]], tmp[[3]]@y.values)    
	PE_scores[[5]][i,names(tmp[[4]])]<-tmp[[4]]

    tmp <- apply_model(meth.t.FF.filtered, 9/10, 400, 0.2)
    AUC.collection[i,6] <- tmp[1]
    nonzero.collection[[6]] <- append(nonzero.collection[[6]],tmp[2])
    false.positive[[6]] <- append(false.positive[[6]], tmp[[3]]@x.values)
    true.positive[[6]] <- append(true.positive[[6]], tmp[[3]]@y.values)    
	PE_scores[[6]][i,names(tmp[[4]])]<-tmp[[4]]

   tmp <- apply_model(meth.t.FF.filtered, 9/10, 400, 0.1)
    AUC.collection[i,7] <- tmp[1]
    nonzero.collection[[7]] <- append(nonzero.collection[[7]],tmp[2])
    false.positive[[7]] <- append(false.positive[[7]], tmp[[3]]@x.values)
    true.positive[[7]] <- append(true.positive[[7]], tmp[[3]]@y.values)    
	PE_scores[[7]][i,names(tmp[[4]])]<-tmp[[4]]

    tmp <- apply_model(meth.t.FF.filtered, 9/10, 100, 0.1)
    AUC.collection[i,8] <- tmp[1]
    nonzero.collection[[8]] <- append(nonzero.collection[[8]],tmp[2])
    false.positive[[8]] <- append(false.positive[[8]], tmp[[3]]@x.values)
    true.positive[[8]] <- append(true.positive[[8]], tmp[[3]]@y.values)    
	PE_scores[[8]][i,names(tmp[[4]])]<-tmp[[4]]

    tmp <- apply_model(meth.t.FF.filtered, 9/10, 200, 0.2)
    AUC.collection[i,9] <- tmp[1]
    nonzero.collection[[9]] <- append(nonzero.collection[[9]],tmp[2])
    false.positive[[9]] <- append(false.positive[[9]], tmp[[3]]@x.values)
    true.positive[[9]] <- append(true.positive[[9]], tmp[[3]]@y.values)    
	PE_scores[[9]][i,rownames(tmp[[4]])]<-tmp[[4]]

    tmp <- apply_model(meth.t.FF.filtered, 4/5, 200, 0.2)
    AUC.collection[i,10] <- tmp[1]
    nonzero.collection[[10]] <- append(nonzero.collection[[10]],tmp[2])
    false.positive[[10]] <- append(false.positive[[10]], tmp[[3]]@x.values)
    true.positive[[10]] <- append(true.positive[[10]], tmp[[3]]@y.values)   
	PE_scores[[10]][i,rownames(tmp[[4]])]<-tmp[[4]]

print(i)
}

boxplot(AUC.collection, col = c("grey70", "royalblue"))
dim(AUC.collection)
colMeans(AUC.collection, na.rm=T)
apply(AUC.collection, 2,median,na.rm=T)


LassoToPlot<-7
MyColVar<-ifelse(meth.t.FF.filtered$condition1=="Case", 1,2)[order(colMeans(PE_scores[[LassoToPlot]], na.rm = T))]
# pdf("PE_scores_model1.pdf")
boxplot(PE_scores[[LassoToPlot]][,order(colMeans(PE_scores[[LassoToPlot]], na.rm = T))], col = c("red", "royalblue")[MyColVar], las = 1, cex.names = 0.2, cex.axis = 0.3, horizontal = T, med.col = c("pink", "lightblue")[MyColVar])
# dev.off()
# pdf("PE_scores_model5.pdf")
plot(colMeans(PE_scores[[LassoToPlot]][,order(colMeans(PE_scores[[LassoToPlot]], na.rm = T))], na.rm = T), col = c("red", "royalblue")[MyColVar], las = 1, cex.names = 0.3, cex.axis = 0.5, beside = T, pch = 16)
text(1:ncol(PE_scores[[LassoToPlot]]), colMeans(PE_scores[[LassoToPlot]][,order(colMeans(PE_scores[[LassoToPlot]], na.rm = T))], na.rm = T), labels = meth.t.FF.filtered$GAD.days[order(colMeans(PE_scores[[LassoToPlot]], na.rm = T))], pos = 3, col = c("red", "royalblue")[MyColVar], las = 1, cex.axis = 0.5,  pch = 16, cex = 0.3)
text(1:ncol(PE_scores[[LassoToPlot]]), colMeans(PE_scores[[LassoToPlot]][,order(colMeans(PE_scores[[LassoToPlot]], na.rm = T))], na.rm = T), labels = rownames(meth.t.FF.filtered)[order(colMeans(PE_scores[[LassoToPlot]], na.rm = T))], pos = 4, col = c("red", "royalblue")[MyColVar], cex.axis = 0.5, pch = 16, cex = 0.3)
# dev.off()

MyColVar<-ifelse(meth.t.FF.filtered$condition1=="Case", 1,2)[order(apply(PE_scores[[LassoToPlot]], 2,median,na.rm = T))]
boxplot(log10(PE_scores[[LassoToPlot]][,order(apply(PE_scores[[LassoToPlot]], 2,median,na.rm = T))]), col = c("red", "royalblue")[MyColVar], las = 1, cex.names = 0.3, cex.axis = 0.3, horizontal = T, med.col = c("pink", "lightblue")[MyColVar])
boxplot(PE_scores[[LassoToPlot]][,order(apply(PE_scores[[LassoToPlot]], 2,median,na.rm = T))], col = c("red", "royalblue")[MyColVar], las = 1, cex.names = 0.3, cex.axis = 0.3, horizontal = T, med.col = c("pink", "lightblue")[MyColVar])
boxplot(apply(PE_scores[[LassoToPlot]], 2,median,na.rm = T)~ meth.t.FF.filtered$condition1)
t.test(apply(PE_scores[[LassoToPlot]], 2,median,na.rm = T)~ meth.t.FF.filtered$condition1)


plot(fa)

rm(meth.t.FF.filtered, AUC.collection, nonzero.collection, true.positive, false.positive)
load("N1_M_excludeSecond.Rda")

 pdf(file = "N1 M-value exclude Second average ROC with random ROC.pdf",  
     width = 8, 
     height = 8) 

# top 200, alpha = 0
colnumber <- nrow(meth.t.FF.filtered) - length(sample(1:nrow(meth.t.FF.filtered), nrow(meth.t.FF.filtered)*4/5)) + 1


x<-matrix(unlist(false.positive[[4]]), ncol = 50, nrow = 11)
y<-matrix(unlist(true.positive[[4]]), ncol = 50, nrow = 11)
i.tmp = 4
false.positive.tmp = false.positive[[i.tmp]]
true.positive.tmp = true.positive[[i.tmp]]
false.positive.reassign.tmp = false.positive.reassign[[i.tmp]]
true.positive.reassign.tmp = true.positive.reassign[[i.tmp]]
maintitle = paste('Top 200 probes, Alpha 0, average AUC = ',mean(AUC.collection[,i.tmp]))

x=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,las=1)
lines(x=c(0,1),y=c(0,1), type="l", pch=22, lty=5,col = '#275eb8')
x.sem=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02)
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02)

par(new=TRUE)

x=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,
     las=1,col ='#6990cf')
x.sem=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02,col = '#6990cf')
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02,col = '#6990cf')

rm(i.tmp,false.positive.tmp,true.positive.tmp,false.positive.reassign.tmp,true.positive.reassign.tmp,maintitle)

# top 100, alpha = 0
i.tmp = 2
false.positive.tmp = false.positive[[i.tmp]]
true.positive.tmp = true.positive[[i.tmp]]
false.positive.reassign.tmp = false.positive.reassign[[i.tmp]]
true.positive.reassign.tmp = true.positive.reassign[[i.tmp]]
maintitle = paste('Top 100 probes, Alpha 0, average AUC = ',mean(AUC.collection[,i.tmp]))

x=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,las=1)
lines(x=c(0,1),y=c(0,1), type="l", pch=22, lty=5,col = '#275eb8')
x.sem=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02)
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02)

par(new=TRUE)

x=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,
     las=1,col ='#6990cf')
x.sem=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02,col = '#6990cf')
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02,col = '#6990cf')

rm(i.tmp,false.positive.tmp,true.positive.tmp,false.positive.reassign.tmp,true.positive.reassign.tmp,maintitle)

# top 200, alpha = 0.2
i.tmp = 3
false.positive.tmp = false.positive[[i.tmp]]
true.positive.tmp = true.positive[[i.tmp]]
false.positive.reassign.tmp = false.positive.reassign[[i.tmp]]
true.positive.reassign.tmp = true.positive.reassign[[i.tmp]]
maintitle = paste('Top 200 probes, Alpha 0.2, average AUC = ',mean(AUC.collection[,i.tmp]))

x=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,las=1)
lines(x=c(0,1),y=c(0,1), type="l", pch=22, lty=5,col = '#275eb8')
x.sem=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02)
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02)

par(new=TRUE)

x=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,
     las=1,col ='#6990cf')
x.sem=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02,col = '#6990cf')
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02,col = '#6990cf')

rm(i.tmp,false.positive.tmp,true.positive.tmp,false.positive.reassign.tmp,true.positive.reassign.tmp,maintitle)

# top 100, alpha = 0.2
i.tmp = 4
false.positive.tmp = false.positive[[i.tmp]]
true.positive.tmp = true.positive[[i.tmp]]
false.positive.reassign.tmp = false.positive.reassign[[i.tmp]]
true.positive.reassign.tmp = true.positive.reassign[[i.tmp]]
maintitle = paste('Top 100 probes, Alpha 0.2, average AUC = ',mean(AUC.collection[,i.tmp]))

x=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,las=1)
lines(x=c(0,1),y=c(0,1), type="l", pch=22, lty=5,col = '#275eb8')
x.sem=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02)
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02)

par(new=TRUE)

x=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,
     las=1,col ='#6990cf')
x.sem=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02,col = '#6990cf')
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02,col = '#6990cf')

rm(i.tmp,false.positive.tmp,true.positive.tmp,false.positive.reassign.tmp,true.positive.reassign.tmp,maintitle)

# top 200, alpha = 0.5
i.tmp = 5
false.positive.tmp = false.positive[[i.tmp]]
true.positive.tmp = true.positive[[i.tmp]]
false.positive.reassign.tmp = false.positive.reassign[[i.tmp]]
true.positive.reassign.tmp = true.positive.reassign[[i.tmp]]
maintitle = paste('Top 200 probes, Alpha 0.5, average AUC = ',mean(AUC.collection[,i.tmp]))

x=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,las=1)
lines(x=c(0,1),y=c(0,1), type="l", pch=22, lty=5,col = '#275eb8')
x.sem=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02)
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02)

par(new=TRUE)

x=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,
     las=1,col ='#6990cf')
x.sem=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02,col = '#6990cf')
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02,col = '#6990cf')

rm(i.tmp,false.positive.tmp,true.positive.tmp,false.positive.reassign.tmp,true.positive.reassign.tmp,maintitle)

# top 100, alpha = 0.5
i.tmp = 6
false.positive.tmp = false.positive[[i.tmp]]
true.positive.tmp = true.positive[[i.tmp]]
false.positive.reassign.tmp = false.positive.reassign[[i.tmp]]
true.positive.reassign.tmp = true.positive.reassign[[i.tmp]]
maintitle = paste('Top 100 probes, Alpha 0.5, average AUC = ',mean(AUC.collection[,i.tmp]))

x=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,las=1)
lines(x=c(0,1),y=c(0,1), type="l", pch=22, lty=5,col = '#275eb8')
x.sem=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02)
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02)

par(new=TRUE)

x=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,
     las=1,col ='#6990cf')
x.sem=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02,col = '#6990cf')
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02,col = '#6990cf')

rm(i.tmp,false.positive.tmp,true.positive.tmp,false.positive.reassign.tmp,true.positive.reassign.tmp,maintitle)

# top 200, alpha = 0.7
i.tmp = 7
false.positive.tmp = false.positive[[i.tmp]]
true.positive.tmp = true.positive[[i.tmp]]
false.positive.reassign.tmp = false.positive.reassign[[i.tmp]]
true.positive.reassign.tmp = true.positive.reassign[[i.tmp]]
maintitle = paste('Top 200 probes, Alpha 0.7, average AUC = ',mean(AUC.collection[,i.tmp]))

x=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,las=1)
lines(x=c(0,1),y=c(0,1), type="l", pch=22, lty=5,col = '#275eb8')
x.sem=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02)
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02)

par(new=TRUE)

x=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,
     las=1,col ='#6990cf')
x.sem=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02,col = '#6990cf')
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02,col = '#6990cf')

rm(i.tmp,false.positive.tmp,true.positive.tmp,false.positive.reassign.tmp,true.positive.reassign.tmp,maintitle)

# top 100, alpha = 0.7
i.tmp = 8
false.positive.tmp = false.positive[[i.tmp]]
true.positive.tmp = true.positive[[i.tmp]]
false.positive.reassign.tmp = false.positive.reassign[[i.tmp]]
true.positive.reassign.tmp = true.positive.reassign[[i.tmp]]
maintitle = paste('Top 100 probes, Alpha 0.7, average AUC = ',mean(AUC.collection[,i.tmp]))

x=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,las=1)
lines(x=c(0,1),y=c(0,1), type="l", pch=22, lty=5,col = '#275eb8')
x.sem=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02)
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02)

par(new=TRUE)

x=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,
     las=1,col ='#6990cf')
x.sem=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02,col = '#6990cf')
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02,col = '#6990cf')

rm(i.tmp,false.positive.tmp,true.positive.tmp,false.positive.reassign.tmp,true.positive.reassign.tmp,maintitle)

# top 200, alpha = 1
i.tmp = 9
false.positive.tmp = false.positive[[i.tmp]]
true.positive.tmp = true.positive[[i.tmp]]
false.positive.reassign.tmp = false.positive.reassign[[i.tmp]]
true.positive.reassign.tmp = true.positive.reassign[[i.tmp]]
maintitle = paste('Top 200 probes, Alpha 1, average AUC = ',mean(AUC.collection[,i.tmp]))

x=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,las=1)
lines(x=c(0,1),y=c(0,1), type="l", pch=22, lty=5,col = '#275eb8')
x.sem=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02)
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02)

par(new=TRUE)

x=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,
     las=1,col ='#6990cf')
x.sem=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02,col = '#6990cf')
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02,col = '#6990cf')

rm(i.tmp,false.positive.tmp,true.positive.tmp,false.positive.reassign.tmp,true.positive.reassign.tmp,maintitle)

# top 100, alpha = 1
i.tmp = 10
false.positive.tmp = false.positive[[i.tmp]]
true.positive.tmp = true.positive[[i.tmp]]
false.positive.reassign.tmp = false.positive.reassign[[i.tmp]]
true.positive.reassign.tmp = true.positive.reassign[[i.tmp]]
maintitle = paste('Top 200 probes, Alpha 1, average AUC = ',mean(AUC.collection[,i.tmp]))

x=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,las=1)
lines(x=c(0,1),y=c(0,1), type="l", pch=22, lty=5,col = '#275eb8')
x.sem=apply(matrix(unlist(false.positive.tmp[which(lapply(false.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.tmp[which(lapply(true.positive.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02)
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02)

par(new=TRUE)

x=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
y=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,mean)
plot(x,y,type='l',
     xlab = 'False positive rate',ylab = 'True positive rate', main = maintitle,
     las=1,col ='#6990cf')
x.sem=apply(matrix(unlist(false.positive.reassign.tmp[which(lapply(false.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
y.sem=apply(matrix(unlist(true.positive.reassign.tmp[which(lapply(true.positive.reassign.tmp,length) ==colnumber)]), ncol = colnumber, byrow = TRUE),2,sd)/sqrt(100)
arrows(x0=x, y0=y-y.sem, x1=x, y1=y+y.sem, code=3, angle=90, length=0.02,col = '#6990cf')
arrows(x0=x-x.sem, y0=y, x1=x+x.sem, y1=y, code=3, angle=90, length=0.02,col = '#6990cf')

rm(i.tmp,false.positive.tmp,true.positive.tmp,false.positive.reassign.tmp,true.positive.reassign.tmp,maintitle)

dev.off()


