####downstream analysis using seqmonk output

library(data.table)

LoadProbeIndex <- function(indexfh) {
    idx <- fread(indexfh, header=TRUE, sep="\t", data.table=FALSE)
    return(idx)
}

GetAvgCount <- function(frlist, indexfh) {
    samlist <- fread(samples.file, header=FALSE, sep="\t", colClasses=c("character", "character"), data.table=FALSE)
    idx <- LoadProbeIndex(indexfh)
    
    for (i in 1:nrow(samlist)) {
         print(paste("reading in", samlist[i,1]))
         fwdfh <- fread(samlist[i,1], header=TRUE, sep="\t", data.table=FALSE)
	 colnames(fwdfh) <- gsub(".cov.gz","",colnames(fwdfh))
	 idxfwdfh <- merge(idx, fwdfh, by=c("Probe","Chromosome","Start","End"), all.x=TRUE)
	 idxfwdfh <- idxfwdfh[order(idxfwdfh$Index),]
	 idxfwdfh <- idxfwdfh[,-c(1:4,6:14)]
	 idxfwdfh_long <- melt(idxfwdfh, id.vars=c("Index"))
	 colnames(idxfwdfh_long) <- c("Index","Sample","Fwdval")
	 
	 print(paste("reading in", samlist[i,2]))
	 revfh <- fread(samlist[i,2], header=TRUE, sep="\t", data.table=FALSE)
	 colnames(revfh) <- gsub(".cov.gz","",colnames(revfh))
	 idxrevfh <- merge(idx, revfh, by=c("Probe","Chromosome","Start","End"), all.x=TRUE)
	 idxrevfh <- idxrevfh[order(idxrevfh$Index),]
	 idxrevfh <- idxrevfh[,-c(1:4,6:14)]
	 idxrevfh_long <- melt(idxrevfh, id.vars=c("Index"))
	 colnames(idxrevfh_long) <- c("Index","Sample","Revval")
	 
	 dt_long <- merge(idxfwdfh_long, idxrevfh_long, by=c("Index","Sample"))
	 dt_long$Total <- dt_long$Fwdval+dt$Revval
	 dt_long$Avgmc <- dt_long$Fwdval/dt_long$Total
	 
	 #dtrt <- dcast(melt(dt_long, id.vars=c("Index","Sample")), Index~Sample+variable)
	 dt <- dcast(dt_long, Index ~ Sample, value.var="Avgmc")
    }
}