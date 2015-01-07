loadFiles <- function(path, group=rep(1:2, each=3)){
	each.rep <- as.vector(table(group))[1]

	# load simulated truth (1000 DS genes)
	annofile <- paste0(path, "true_exons_simulation.txt")
	ds <- read.columns(annofile, required.col=c("exon", "status_exon"))
	ds <- ds$exon[ds$status_exon==1]
	ds <- as.character(ds[!is.na(ds)])

	# load counts (and filter low expressed features)
	countfile <- paste0(path, "/featureCounts/featureCounts.txt")
	counts <- read.table(countfile)
	exonid <- counts[,1]
	counts <- counts[,2:7]
	keep <- rowSums(cpm(counts)>1)>= each.rep # standard filter
	# keep <- rep(TRUE, nrow(counts)) # no filter
	
	# load voom results 
	voomfile <- paste0(path, "voom_exon_results.txt")
	voom <- read.table(voomfile, header=TRUE)[,-3]
	colnames(voom)[1:2] <- c("GeneID", "ExonID")
	
	out <- list(ds=ds, counts=counts, exonid=exonid, keep=keep, voom=voom)
}

makePlot <- function(checkid, results, cutoff=0.05){
	ds <- results$ds
	counts <- results$counts
	exonid <- results$exonid
	keep <- results$keep
	voom <- results$voom
	
	# subset for specific gene
	checkid <-  strsplit2(checkid, split=":")[,1]
	ii <- grep(checkid, exonid)
	counts.i <- counts[ii,]
	exonid.i <- exonid[ii]
	keep.i <- keep[ii]
	o <- order(exonid.i)
	counts.i <- counts.i[o,]
	exonid.i <- exonid.i[o]
	keep.i <- keep.i[o]
	
	# plot counts before filtering
	par(mfrow=c(1,2))
	filtered.notds <- which(!keep.i & !exonid.i %in% ds)
	filtered.ds <- which(!keep.i & exonid.i %in% ds)
	for (i in 1) plot(1:nrow(counts.i), log2((counts.i+0.05)[,i]), type="l", col=group[i], ylab="log2-count", xlab="", main="Filtering", ylim=c(-6,15), xaxt="n") #, ...)
	axis(side=1, at=1:nrow(counts.i), labels=exonid.i, las=2)
	for (i in 2:ncol(counts.i)) lines(1:nrow(counts.i), log2((counts.i+0.05)[,i]), col=group[i])
	abline(h=log2(0.05), col=max(group)+1)
	text(0, log2(0.05), "count of zero", pos=4, col=max(group)+1)
	for (i in 1:ncol(counts.i)) points(filtered.notds, log2((counts.i+0.05)[filtered.notds,i]), pch=15, col=max(group)+1)
	for (i in 1:ncol(counts.i)) points(filtered.ds, log2((counts.i+0.05)[filtered.ds,i]), pch=15, col=max(group)+2)
	legend("topright", col=c(max(group)+1:2), pch=15:16, c("filtered out", "filtered out and generated as DS"), bty="n")
	
	# plot counts after filtering, and with final results
	counts.i <- counts.i[keep.i,]
	exonid.i <- exonid.i[keep.i]
	m <- match(exonid.i, voom$ExonID)
	detected.i <- voom[m,"FDR"]<cutoff
	truepos <- which(detected.i & exonid.i %in% ds)
	falsepos <- which(detected.i & !exonid.i %in% ds)
	falseneg <- which(!detected.i & exonid.i %in% ds)
	for (i in 1) plot(1:nrow(counts.i), log2((counts.i+0.05)[,i]), type="l", col=group[i], ylab="log2-count", xlab="", main="voom results", ylim=c(-6,15), xaxt="n") #, ...)
	axis(side=1, at=1:nrow(counts.i), labels=exonid.i, las=2)
	for (i in 2:ncol(counts.i)) lines(1:nrow(counts.i), log2((counts.i+0.05)[,i]), col=group[i])
	for (i in 1:ncol(counts.i)) points(truepos, log2((counts.i+0.05)[truepos,i]), pch=15, col=max(group)+3)
	for (i in 1:ncol(counts.i)) points(falsepos, log2((counts.i+0.05)[falsepos,i]), pch=16, col=max(group)+4)
	for (i in 1:ncol(counts.i)) points(falseneg, log2((counts.i+0.05)[falseneg,i]), pch=17, col=max(group)+5)
	legend("topright", col=c(max(group)+3:5), pch=15:17, c("true positive", "false positive", "false negative"), bty="n")
}