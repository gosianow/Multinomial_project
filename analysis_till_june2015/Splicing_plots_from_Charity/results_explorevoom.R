# What do voomex false discoveries look like
# 30 June 2014 (Last updated 2 July 2014)
options(digits=2)
options(width=200)

library(limma)
library(edgeR)
source("results_explorevoomfunctions.R")

# set simulation type
sim <- "hsapiens_V2"
path <- paste0("./simulation_", sim, "/")

results <- loadFiles(path)
results$voom$FD <- !results$voom$ExonID %in% results$ds
checkids <- as.character(results$voom$ExonID[results$voom$FD])

# most significant exons come from genes with greater no. of exons
nexons <- rowsum(rep(1, nrow(results$voom)), results$voom$GeneID)
m <- match(results$voom$GeneID, rownames(nexons))
results$voom$NExons <- nexons[m]
pdf("./results/voom_significance_by_nexons.pdf", width=10)
bin.size <- 10
nbins <- 50
nexons.bin <- results$voom$NExons[1:(bin.size*nbins)]
nexons.bin <- matrix(nexons.bin, nrow=bin.size)
boxplot(nexons.bin, 
	ylab="Total no. of exons within associated gene",
	xlab=paste0("Top ranking exons by adj p-value (bin size = ", bin.size, ")"), 
	las=2)
dev.off()

for (i in 1:20) {
	checkid <- checkids[i]
	checkid <- strsplit2(checkid, split=":")[1]
	plotname <- paste0("./results/", sim, "_check_", checkid, ".pdf")
	pdf(plotname, width=20, height=8)
	par(mar=c(12,4,4,0))
	makePlot(checkid, results)
	dev.off()
}


#%#%#%#%#%#%#%#%##%#%#%#%#%#% EXTRA STUFF #%#%#%#%#%#%#%#%##%#%#%#%#%#% 

# A fairly large proportion of exons generated as DS have zero/low counts (and are filtered out downstream)
# -> loss of power to detect 'true positives' downstream
ds <- results$ds
counts <- results$counts
exonid <- results$exonid
rownames(counts) <- exonid
# zero count DS exons
# load featureCounts counts and "truth"
ds.counts <- counts[exonid %in% ds,]
rownames(ds.counts) <- exonid[exonid %in% ds]
ds.zero <- rowSums(ds.counts==0)
table(ds.zero==6) # 434 exons out of 10256 are all zero counts (4% of DS exons are not expressed at all)
ds.low <- rowSums(ds.counts<=3)
table(ds.low==6) # 1029 exons out of 10256 are all low counts (approx 10% are filtered out)

# Are these lowly expressed DS exons detected as significant by dexseq?
# Result: only about 6% of them are (loss of power)
ds.low.exons <- names(ds.low)[ds.low==6]
dexseq.exon <- read.columns("/home/Shared_taupo/tmp/Simulation/simulation_hsapiens_V2/dexseq_exon_results.txt", 
	required.col=c("groupID", "featureID", "pvalue", "padj", "log2fold_Condition2_Condition1"))
dexseq.exon$featureID <- gsub(":E", ":", dexseq.exon$groupID)
dexseq.exon$groupID <- strsplit2(dexseq.exon$groupID, split=":")[,1]
colnames(dexseq.exon) <- c("GeneID", "ExonID", "P.Value", "FDR", "logFC")
dexseq.ds.low.exons <- dexseq.exon[dexseq.exon$ExonID %in% ds.low.exons,]
dim(dexseq.ds.low.exons) # 1029 exons with low expression in all samples
table(is.na(dexseq.ds.low.exons$FDR)) # 889 exons are  not testable by dexseq 
table(dexseq.ds.low.exons$FDR<0.05) # 64 exons are significant at FDR of 0.05 

# Many ds genes have 40% or more exons set as "ds" so that it is hard to fit a "pattern" to the gene. 
# Is this realistic? 
nexons <- strsplit2(exonid, split=":")[,1] 
nexons <- rowsum(rep(1, length(nexons)), nexons)
nexons.ds <- strsplit2(ds, split=":")[,1] 
nexons.ds <- rowsum(rep(1, length(nexons.ds)), nexons.ds)
m <- match(rownames(nexons.ds), rownames(nexons))
nexons.ds <- as.data.frame(cbind(nexons[m], nexons.ds))
colnames(nexons.ds) <- c("all", "ds")
nexons.ds$prop.ds <- nexons.ds$ds/nexons.ds$all 
pdf("./results/proportion_of_DS_exons_within_gene.pdf", width=10, height=5)
par(mfrow=c(1,2))
hist(nexons.ds$prop.ds, xlab="Proportion of exons simulated as DS within a gene", main="")
plot(nexons.ds$all, nexons.ds$prop.ds, xlab="Total number of exons within a gene", ylab="Proportion of exons simulated as DS within a gene")
dev.off()

# 8 genes aren't actually differentially spliced, as they don't have ds exons.
# (These genes should be removed from list of ds genes -- what about transcripts, should any be removed?)
missinggenes <- setdiff(ds$Gene, rownames(nexons.ds))
for (i in 1:length(missinggenes)) cat(missinggenes[i], "\n", ds$Exon[grep(missinggenes[i], ds$Exon)], "\n")

# Check what transcripts look like where "pattern inversion" takes place.
annofile <- paste0(path, "true_exons_simulation.txt")
ds.anno <- read.columns(annofile, required.col=c("Gene", "Transcript", "exon", "status", "status_exon"))
# for one particular gene:
tmp <- ds.anno[grep("ENSG00000197299", ds.anno$Gene),]
unique(tmp$Transcript)
nexons.ds[grep("ENSG00000197299", rownames(nexons.ds)),]
tmp
tmp <- ds.anno[grep("ENSG00000061273", ds.anno$Gene),]
unique(tmp$Transcript)
nexons.ds[grep("ENSG00000061273", rownames(nexons.ds)),]
tmp

# Are false positives detected by voom, also detected by dexseq?
voom <- results$voom
dexseq <- dexseq.exon
falsepos <- as.character(voom$ExonID[voom$FDR<0.05 & voom$FD==TRUE]) # 287 false positives detected by voom
table(dexseq$FDR[dexseq$ExonID %in% falsepos]<0.05) # 106 also detected, 10 not detected, 171 are not estimable 

# FDs sorted by percentage DS 
m <- match(voom$GeneID, rownames(nexons.ds))
voom$PropDS <- nexons.ds$prop.ds[m] 
voom$PropDS[is.na(voom$PropDS)] <- 0
pdf("./results/voom_prop_ds_in_FD_exons_ranked_by_significance.pdf")
bin.size <- 10
nbins <- 20
tmp <- voom$PropDS[voom$FD]
tmp <- tmp[1:(bin.size*nbins)]
tmp <- matrix(tmp, nrow=bin.size)
boxplot(tmp, 
	ylab="Associated proportion of diff. spliced exons within gene",
	xlab=paste0("Top ranking falsely discovered exons by adj p-value (bin size = ", bin.size, ")"),
	las=2) 
dev.off()


source("results_plotfunctions.R")

# Size plot for nonDS exons in non-null sim (results similar to null sim)
table <- voom[voom$PropDS==0,]
pdf("./results/size_for_nonDS_exons.pdf")
plot(0:1, type="n", ylim=c(0,2), xlim=c(0,1))
cutoff <- c(0.01, 0.05, 0.1)
abline(v=cutoff, col=2:4, lty=2)
propP(pvalues=table$P.Value, p.cutoff=cutoff, col=2:4)
dev.off()


pdf("./results/fd_by_prop_ds.pdf")
plot(0:1, type="n", ylim=c(0,0.2), xlim=c(0,0.6))
cutoff <- c(0.01, 0.05, 0.1)
abline(v=cutoff, col=2:4, type=2)
i <- (0:10)/10
for (j in 1:5) {
table <- voom[voom$PropDS>i[j] & voom$PropDS<=i[j+1],]
propFD.TP(table=table, truth=ds, level="ExonID", fdr.cutoff=cutoff, col=2:4, pch=20+j)
}
for (j in 6:10) {
table <- voom[voom$PropDS>i[j] & voom$PropDS<=i[j+1],]
propFD.TP(table=table, truth=ds, level="ExonID", fdr.cutoff=cutoff, col=5:7, pch=15+j)
}
legend("topright", pch=rep(21:25, 2), legend=c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0"), col=rep(c(2,5), each=5))
dev.off()


pdf("./results/prop_ds_dist_of_fds.pdf")
hist(voom$PropDS[voom$FD & voom$FDR<0.05])
dev.off()



# Plot true and false postives according to %DS
a <- matrix(nrow=11, ncol=5)
colnames(a) <- c("TP", "FP", "TN", "FN", "TotalExons")
rownames(a) <- c("0.0", "0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
cutoff <- 0.05
j <- (0:10)/10

i <- 1
tt <- voom[voom$PropDS==j[i],]
a[i,"TotalExons"] <- nrow(tt)
a[i,"TP"] <- sum(tt$FDR<cutoff & tt$ExonID %in% ds)/nrow(tt)
a[i,"FP"] <- sum(tt$FDR<cutoff & !tt$ExonID %in% ds)/nrow(tt)
a[i,"TN"] <- sum(tt$FDR>cutoff & !tt$ExonID %in% ds)/nrow(tt)
a[i,"FN"] <- sum(tt$FD>cutoff & tt$ExonID %in% ds)/nrow(tt)

for (i in 2:11) {
	tt <- voom[voom$PropDS>j[i-1] & voom$PropDS<=j[i],]
	a[i,"TotalExons"] <- nrow(tt)
	a[i,"TP"] <- sum(tt$FDR<cutoff & tt$ExonID %in% ds)/nrow(tt)
	a[i,"FP"] <- sum(tt$FDR<cutoff & !tt$ExonID %in% ds)/nrow(tt)
	a[i,"TN"] <- sum(tt$FDR>cutoff & !tt$ExonID %in% ds)/nrow(tt)
	a[i,"FN"] <- sum(tt$FDR>cutoff & tt$ExonID %in% ds)/nrow(tt)
}

# Good and bad false discoveries
ds.genes <- unique(strsplit2(ds, split=":")[,1])
a2 <- cbind(a, FPgood=NA, FPbad=NA)
i <- 1
tt <- voom[voom$PropDS==j[i],]
a2[i,"FPgood"] <- sum(tt$FDR<cutoff & !tt$ExonID %in% ds & tt$GeneID %in% ds.genes)/nrow(tt)
a2[i,"FPbad"] <- sum(tt$FDR<cutoff & !tt$ExonID %in% ds & !tt$GeneID %in% ds.genes)/nrow(tt)

for (i in 2:11) {
	tt <- voom[voom$PropDS>j[i-1] & voom$PropDS<=j[i],]
	a2[i,"FPgood"] <- sum(tt$FDR<cutoff & !tt$ExonID %in% ds & tt$GeneID %in% ds.genes)/nrow(tt)
	a2[i,"FPbad"] <- sum(tt$FDR<cutoff & !tt$ExonID %in% ds & !tt$GeneID %in% ds.genes)/nrow(tt)
}

pdf("./results/fd_tp_by_prop_ds.pdf", width=21)
par(mfrow=c(1,3))
barplot(t(a[,-5]), beside=TRUE, col=1:4, ylim=c(0,1.1), las=2)
legend("topright", col=1:4, c("TP", "FP", "TN", "FN"), pch=15)
abline(h=cutoff, col="red", lty=2)
bp <- barplot(t(a[,1:2]), beside=TRUE, col=1:2, las=2)
legend("topright", col=1:2, c("TP", "FP"), pch=15, bty="n")
text(bp[2,], a[,"FP"], round(a[,"FP"]/rowSums(a[,c("FP", "TP")]),2), col="red", pos=3)
bp <- barplot(t(a2[,c(1,6,7)]), beside=TRUE, col=c(1,5,6), las=2)
legend("topright", col=c(1,5,6), c("TP", "FPgood", "FPbad"), pch=15, bty="n")
dev.off()



# Plot number of exons for true and false postives according to %DS
pdf("./results/fp_tp_by_prop_ds_and_nexons.pdf", width=15)
plot(1:11, ylim=c(0,80), xlim=c(1,11), type="n",
	ylab="Number of exons in associated gene", xlab="", xaxt="n")
axis(side=1, at=1:11, las=2, 
	labels=c("0.0", "0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0"))	
for (i in 1:11) {
	if(i==1) tt <- voom[voom$PropDS==j[i],]
	if(i>1)	tt <- voom[voom$PropDS>j[i-1] & voom$PropDS<=j[i],]
	boxplot(tt$NExons[tt$FDR<cutoff & tt$ExonID %in% ds],
		add=TRUE, at=(i-0.1), boxwex=0.3, col="grey")
	boxplot(tt$NExons[tt$FDR<cutoff & !tt$ExonID %in% ds],
		add=TRUE, at=(i+0.1), boxwex=0.3, col="red")
}
legend("topright", col=c("grey", "red"), c("True positives", "False positives"), pch=15, bty="n")
dev.off()



# Find how many genes and exons suffer from pattern inversion.
# For each gene, find the proportion of FP and FN (and TP, TN)
# b <- voom[,c("GeneID", "NExons", "PropDS")][!duplicated(voom$GeneID),]
# b <- as.data.frame(cbind(b, tp=NA, fp=NA, tn=NA, fn=NA))
# for (i in 1:nrow(b)) {
	# b$tp[i] <- sum(voom$GeneID==b$GeneID[i] & voom$FDR<cutoff & voom$ExonID %in% ds)
	# b$fp[i] <- sum(voom$GeneID==b$GeneID[i] & voom$FDR<cutoff & !voom$ExonID %in% ds)
	# b$tn[i] <- sum(voom$GeneID==b$GeneID[i] & voom$FDR>=cutoff & !voom$ExonID %in% ds)
	# b$fn[i] <- sum(voom$GeneID==b$GeneID[i] & voom$FDR>=cutoff & voom$ExonID %in% ds)
# }
# b[,c("tp", "fp", "tn", "fn")] <- b[,c("tp", "fp", "tn", "fn")]/b$NExons 
# save(b, file="./results/genewise_stats.RData")
load("./results/genewise_stats.RData")

o <- order(b$PropDS, decreasing=TRUE)
b <- b[o,]
bb <- b$fp+b$fn

bbi <- bb>=0.5

pdf("./results/pattern_inversion.pdf", width=10, height=10)
par(mfrow=c(2,2))
ii <- b$PropDS>0
plot(b$NExons[ii], b$PropDS[ii], cex=2*(b$fp[ii]/max(b$fp[ii])), 
	col="red", main="Size as false pos rate")
plot(b$NExons[ii], b$PropDS[ii], cex=2*(b$fn[ii]/max(b$fn[ii])),
	col="blue", main="Size as false neg rate")
plot(b$NExons[ii], b$PropDS[ii], cex=2*(b$tp[ii]/max(b$tp[ii])), 
	col="black", main="Size as true pos rate")
plot(b$NExons[ii], b$PropDS[ii], cex=2*(b$tn[ii]/max(b$tn[ii])),
	col="green", main="Size as true neg rate")
dev.off()


pdf("./results/pattern_inversion_2.pdf", width=10, height=10)
par(mfrow=c(2,2))
plot(b$NExons[ii], b$PropDS[ii], cex=(1-bb[ii])^2, main="Size as 1 - (false pos and neg rate)")
plot(b$NExons[ii], b$PropDS[ii], cex=((b$NExons*b$fp)[ii])/2, main="Size as number of false positives")
dev.off()

#%$%#%$%#%$%#%$%#%$%#%$%#%$%#%$%#%$%#%$%#%$%#%$%#%$%
geneid <- strsplit2(rownames(counts), split=":")[,1]
amean.exon <- rowMeans(counts)
amean <- rowsum(amean.exon, geneid, sort=FALSE)/rowsum(rep(1, length(geneid)), geneid, sort=FALSE) 

cutoff <- 0.05
voom$Cat <- "TP"
# voom$Cat[!voom$ExonID %in% ds & voom$FDR<cutoff] <- "FP"
voom$Cat[!voom$ExonID %in% ds & voom$FDR<cutoff & voom$GeneID %in% ds.genes] <- "FPgood"
voom$Cat[!voom$ExonID %in% ds & voom$FDR<cutoff & !voom$GeneID %in% ds.genes] <- "FPbad"
voom$Cat[!voom$ExonID %in% ds & voom$FDR>=cutoff] <- "TN"
voom$Cat[voom$ExonID %in% ds & voom$FDR>=cutoff] <- "FN"
genes.fd <- unique(voom$GeneID[voom$FD])

ngenes <- 100
tmp <- voom[voom$GeneID %in% genes.fd[1:ngenes],]
TP <- rowsum(as.numeric(tmp$Cat=="TP"), tmp$GeneID, reorder=FALSE)
FP <- rowsum(as.numeric(tmp$Cat=="FP"), tmp$GeneID, reorder=FALSE)
TN <- rowsum(as.numeric(tmp$Cat=="TN"), tmp$GeneID, reorder=FALSE)
FN <- rowsum(as.numeric(tmp$Cat=="FN"), tmp$GeneID, reorder=FALSE)
cat <- t(cbind(TP, FP, TN, FN))
legend <- c("TP", "FP", "TN", "FN")

width <- tmp$PropDS[!duplicated(tmp$GeneID)]

width1 <- (width+0.5)^4

m <- match(rownames(FP), rownames(amean))
width2 <- amean[m,]
width2 <- width2^0.25

pdf("./results/cat_top100geneswithFDs.pdf", height=14, width=14)
par(mfrow=c(1,3))
par(mar=c(4,7,2,1))
barplot(cat, col=1:4, horiz=TRUE, legend=legend, las=2, cex.names=0.6, xlab="Number of exons", 
	main="Top 100 genes with false positives")
barplot(cat, col=1:4, horiz=TRUE, legend=legend, las=2, cex.names=0.6, xlab="Number of exons", 
	main="(bar width prop. to % of postives within gene)", width=width1)
barplot(cat, col=1:4, horiz=TRUE, legend=legend, las=2, cex.names=0.6, xlab="Number of exons", 
	main="(bar width prop. to average gene expression)",	width=width2)
dev.off()
