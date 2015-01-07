# 18 June 2014 (Last updated 20 June 2014)
source("results_plotfunctions.R")
cutoff <- c(0.01, 0.05, 0.1)
col.cutoff <- c("blue", "red", "limegreen")

# pdist
plotname <- paste0("./results/pdist_", sim, ".pdf")
pdf(plotname, height=5, width=5)
par(mfrow=c(2,2))
hist(voom.exon$P.Value, xlab="P-value", main="voom-diffSplice (exon)", col="pink")
hist(voom.gene$P.Value, xlab="P-value", main="voom-diffSplice (gene)", col="deeppink")
hist(dexseq.exon$P.Value, xlab="P-value", main="DEXSeq (exon)", col="turquoise")
hist(cuffdiff.iso$P.Value, xlab="P-value", main="CuffDiff2 (transcript)", col="royalblue")
dev.off()



# FOR NULL SIM
if (length(grep("null", sim))!=0) {
	
	# proportion of raw p-value below a cut-off
	plotname <- paste0("./results/size_", sim, ".pdf")
	pdf(plotname, height=4, width=6)
	par(mar=c(4,10,2,1))
	plot(1:4, xlim=c(0,0.15), ylim=c(0.8,4.2),
		type="n", ylab="", yaxt="n", xlab="Proportion of features below p-value cutoff", xaxt="n",
		main="")
	axis(side=2, at=1:4, labels=c("DEXSeq (exon)", "voom-diffSplice (exon)", "voom-diffSplice (gene)", "CuffDiff2 (transcript)"),
	 	las=1, tick=F)
	axis(side=1, at=cutoff, labels=cutoff, las=2)
	abline(v=cutoff, lty=2, col=col.cutoff)
	propP(dexseq.exon, cutoff, at=1, col=col.cutoff, cex=2)
	propP(voom.exon, cutoff, at=2, col=col.cutoff, cex=2)
	propP(voom.gene, cutoff, at=3, col=col.cutoff, cex=2)
	propP(cuffdiff.iso, cutoff, at=4, col=col.cutoff, cex=2)
	dev.off()

}


# FOR NONNULL SIM
if (length(grep("null", sim))==0) {
	
	# proportion of false discoveries below a FDR cut-off versus true positive rate
	plotname <- paste0("./results/fd_tp_", sim, ".pdf")
	pdf(plotname, height=5, width=5)
	par(mar=c(4,4,2,1))
	plot(1:5, xlim=c(0,1), ylim=c(0,1),
		type="n", ylab="True positive rate", xlab="False discovery rate", xaxt="n")
	axis(side=1, at=cutoff, labels=cutoff, las=2)
	abline(v=cutoff, lty=2, lwd=2, col=col.cutoff)
	propFD.TP(cuffdiff.iso, ds, level="Transcript", cutoff, col=col.cutoff, pch=21, cex=1.5)
	propFD.TP(dexseq.exon, ds, level="Exon", cutoff, col=col.cutoff, pch=22, cex=1.5)
	propFD.TP(dexseq.gene, ds, level="Gene", cutoff, col=col.cutoff, pch=23, cex=1.5)
	propFD.TP(voom.exon, ds, level="Exon", cutoff, col=col.cutoff, pch=24, cex=1.5)
	propFD.TP(voom.gene, ds, level="Gene", cutoff, col=col.cutoff, pch=25, cex=1.5)
	legend("topright", legend=c("CuffDiff2 (gene)", "DEXSeq (exon)", "DEXSeq (gene)", "voom-diffSplice (exon)", "voom-diffSplice (gene)"),
		bty="n", cex=0.7, pch=21:25)
	dev.off()

}