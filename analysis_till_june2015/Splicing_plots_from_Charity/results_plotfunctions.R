# 18 June 2014 (Last updated 23 June 2014)

# proportion of raw p-value below a cut-off
propP <- function(pvalues, p.cutoff=0.05, plot.points=TRUE, pch=21, col="black", at=1, ...){
	if(is.data.frame(pvalues)) pvalues <- pvalues$P.Value
	pvalues <- pvalues[!is.na(pvalues)]
	pp <- NULL
	for (i in 1:length(p.cutoff)) {
	pp <- c(pp, sum(pvalues<p.cutoff[i]))
	}
	pp <- pp/length(pvalues)
	if(plot.points) {
	highlight <- pp > p.cutoff
	bg <- col
	bg[highlight] <- "white"
	points(pp, rep(at, length(pp)), col=col, pch=pch, bg=bg, ...)	
	}
	out <- pp
}



# proportion of false discoveries below a FDR cut-off versus true positive rate
propFD.TP <- function(table, truth, level, fdr.cutoff=0.05, plot.points=TRUE, pch=21, col="black", ...){
	id <- table[,grep(level, colnames(table))]
	fdr <- table$FDR
	keep <- !is.na(fdr) 
	id <- id[keep]
	fdr <- fdr[keep]
	if(is.list(truth)) truth <- truth[[level]]
	truth <- truth[!is.na(truth)]
	truth <- unique(truth) 
	tp <- NULL
	fd <- NULL
	for (i in 1:length(fdr.cutoff)) {
		ndiscoveries <- sum(fdr<fdr.cutoff[i])
		nfalsediscoveries <- sum(fdr<fdr.cutoff[i] & !id %in% truth) 
		detected <- as.character(id[fdr<fdr.cutoff[i]])
		detected <- detected[!is.na(detected)]
		ntruepos <- sum(truth %in% detected)
		tp <- c(tp, ntruepos/length(truth))
		fd <- c(fd, nfalsediscoveries/ndiscoveries)	
	}	
	if(plot.points) {
		lines(fd, tp, type="l", col="black")	
		highlight <- fd > fdr.cutoff
		bg <- col
		bg[highlight] <- "white"
		points(fd, tp, col=col, pch=pch, bg=bg, ...)	
	}
	out <- list(tp=tp, fd=fd)
	out
}


# # plotFD function(table, truth, level, plot.lines=TRUE, ...){
	# o <- order(table$FDR)
	# table <- table[o,]
	# id <- table[,grep(level, colnames(table))]
	
# }


