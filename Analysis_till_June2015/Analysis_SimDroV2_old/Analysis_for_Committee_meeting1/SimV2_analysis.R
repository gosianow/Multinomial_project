#######################################################
# for Committee meeting
# Created 14 July 2014 / 
# Last updated 16 July 2014 (02 Sep 2014)

#######################################################
# BioC 2.14


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2/")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata

#######################################################
# compare featureCounts & htseq 
#######################################################

library(DEXSeq)

ecs <- DEXSeqDataSetFromHTSeq(countfiles = paste0("htseqCounts/htseq_samp_", metadata$SampleName1, ".txt"), sampleData= data.frame(row.names = metadata$SampleName, condition = metadata$condition), design = ~sample + exon + condition:exon)

counts <- featureCounts(ecs)
colnames(counts) <- paste0("htseq", 1:6)
htseq <- data.frame(E= gsub("E", "",rownames(counts)), counts)


fc <- read.table("featureCounts/featureCounts.txt")
colnames(fc) <- c("E", paste0("fc", 1:6))


cnts <- merge(htseq, fc, by="E", all=TRUE)

write.table(cnts, "PLOTS/Compare_htseq_and_fc.xls", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )


pdf("PLOTS/Compare_htseq_and_fc.pdf")

smoothScatter(log(cnts$htseq1+1), log(cnts$fc1+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")
smoothScatter(log(cnts$htseq2+1), log(cnts$fc2+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")
smoothScatter(log(cnts$htseq3+1), log(cnts$fc3+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")
smoothScatter(log(cnts$htseq4+1), log(cnts$fc4+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")
smoothScatter(log(cnts$htseq5+1), log(cnts$fc5+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")
smoothScatter(log(cnts$htseq6+1), log(cnts$fc6+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")

dev.off()

#######################################################
# compare featureCounts_v3 & htseq 
#######################################################

library(DEXSeq)

ecs <- DEXSeqDataSetFromHTSeq(countfiles = paste0("htseqCounts/htseq_samp_", metadata$SampleName1, ".txt"), sampleData= data.frame(row.names = metadata$SampleName, condition = metadata$condition), design = ~sample + exon + condition:exon)

counts <- featureCounts(ecs)
colnames(counts) <- paste0("htseq", 1:6)
htseq <- data.frame(E= gsub("E", "",rownames(counts)), counts)


fc <- read.table("featureCounts_v3/fc.txt", header = T)
colnames(fc) <- c("E", paste0("fc", 1:6))


cnts <- merge(htseq, fc, by="E")

write.table(cnts, "PLOTS/Compare_htseq_and_fc_v3.xls", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )


pdf("PLOTS/Compare_htseq_and_fc_v3.pdf")

smoothScatter(log(cnts$htseq1+1), log(cnts$fc1+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")
smoothScatter(log(cnts$htseq2+1), log(cnts$fc2+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")
smoothScatter(log(cnts$htseq3+1), log(cnts$fc3+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")
smoothScatter(log(cnts$htseq4+1), log(cnts$fc4+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")
smoothScatter(log(cnts$htseq5+1), log(cnts$fc5+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")
smoothScatter(log(cnts$htseq6+1), log(cnts$fc6+1), nrpoints = Inf, xlab="htseq", ylab="fc")
abline(a=0, b=1, col="red")

dev.off()



#######################################################
# run DM on  HTSeq/DEXSeq counts:
#######################################################

library(DEXSeq)

ecs <- DEXSeqDataSetFromHTSeq(countfiles = paste0("htseqCounts/htseq_samp_", metadata$SampleName1, ".txt"), sampleData= data.frame(row.names = metadata$SampleName, condition = metadata$condition), design = ~sample + exon + condition:exon)

counts <- featureCounts(ecs)
gene.id <- unlist(lapply(seq(rownames(counts)), function(i){e <- rownames(counts)[i];strsplit(e, ":")[[1]][1] }))
ete.id <-  rownames(counts)


library("edgeR")
# create DGEList object
dge <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))
# filter counts by log-cpm!!
keep <- rowSums(cpm(dge)>1) >= 3
dge <- dge[keep,]

## add 1 to counts
# dge$counts <- dge$counts + 1
dge$counts[ dge$counts == 0 ] <- 1

## save dge
dir.create("DM/htseq/", showWarnings=F, recursive=T)
name1 <- "htseq"

write.table(data.frame(dge$genes,dge$counts), paste0("DM/htseq/dge_counts_",name1,".xls"), quote=F, sep="\t", row.names=F, col.names=T)




# run DM G2

design <- model.matrix(~condition, data=metadata)
design

source("/home/gosia/R/R_Multinomial_project/Analysis_SimDroV2/Analysis_for_Committee_meeting1/dirmult.R")
source("/home/gosia/R/R_Multinomial_project/Analysis_SimDroV2/Analysis_for_Committee_meeting1/glmFitDM.R")


g2fit <- g2FitDM(dge)

g2lrt <- g2LRTDM(g2fit)


head(g2lrt$table)

save(dge, g2fit, g2lrt, file=paste0("DM/htseq/DM_",name1,"_results.RData"))

write.table(g2lrt$table, paste0("DM/htseq/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)



#######################################################
# run DM on featureCounts data
#######################################################

library(limma)

fc <- read.table("featureCounts/featureCounts.txt")

counts <- fc[,2:7]
colnames(counts) <- metadata$SampleName1
gene.id <- strsplit2(fc[,1], ":")[,1]
ete.id <- fc[,1]


library("edgeR")
# create DGEList object
dge.org <- dge <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))

name1 <- "fc"
dge.org$counts[ dge.org$counts == 0 ] <- 1
write.table(data.frame(dge.org$genes,dge.org$counts), paste0("DM/fc/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)

# filter counts by log-cpm 
keep <- rowSums(cpm(dge)>1) >= 3
dge <- dge[keep,]

## add 1 to counts
# dge$counts <- dge$counts + 1
dge$counts[ dge$counts == 0 ] <- 1

# normalisation
dge <- calcNormFactors(dge)


## save dge
dir.create("DM/fc", showWarnings=F, recursive=T)
name1 <- "fc"

write.table(data.frame(dge$genes,dge$counts), paste0("DM/fc/dge_counts_",name1,".xls"), quote=F, sep="\t", row.names=F, col.names=T)



# run DM G2

design <- model.matrix(~condition, data=metadata)
design

source("/home/gosia/R/R_Multinomial_project/Analysis_SimV2/Analysis_for_Committee_meeting1/dirmult.R")
source("/home/gosia/R/R_Multinomial_project/Analysis_SimV2/Analysis_for_Committee_meeting1/glmFitDM.R")


g2fit <- g2FitDM(dge)

g2lrt <- g2LRTDM(g2fit)


head(g2lrt$table)

write.table(g2lrt$table, paste0("DM/fc/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)

save(dge, g2fit, g2lrt, file=paste0("DM/fc/DM_",name1,"_results.RData"))


#######################################################
# run DM on MISO data
#######################################################


library(limma)

miso <- read.table("MISOCounts/MISO_counts_all_samps.txt", header = T)

counts <- miso[,2:7]
gene.id <- strsplit2(miso[,1], ":")[,1]
ete.id <- miso[,1]


library("edgeR")
# create DGEList object
dge <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))
# filter counts by log-cpm 
keep <- rowSums(cpm(dge)>1) >= 3
dge <- dge[keep,]


# for MISO: filter events with only 0s: 0,0,0,0,0
keep2 <- unlist( lapply(strsplit(strsplit2(dge$genes$ete.id, ":")[,2], ","), function(ev){ sum(as.numeric(ev)) > 0 }) )
dge <- dge[keep2,]


## add 1 to counts
dge$counts[ dge$counts == 0 ] <- 1

# normalisation
dge <- calcNormFactors(dge)

## save dge
dir.create("DM/miso", showWarnings=F, recursive=T)
name1 <- "miso"

write.table(data.frame(dge$genes,dge$counts), paste0("DM/miso/dge_counts_",name1,".xls"), quote=F, sep="\t", row.names=F, col.names=T)




# run DM G2

design <- model.matrix(~condition, data=metadata)
design

source("/home/gosia/R/R_Multinomial_project/Analysis_SimV2/Analysis_for_Committee_meeting1/dirmult.R")
source("/home/gosia/R/R_Multinomial_project/Analysis_SimV2/Analysis_for_Committee_meeting1/glmFitDM.R")


g2fit <- g2FitDM(dge)

g2lrt <- g2LRTDM(g2fit)

head(g2lrt$table)

write.table(g2lrt$table, paste0("DM/miso/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)




#######################################################
# calculate gene-level p-values for DEXSeq with Simes' method
#######################################################

dexseq.ex <- read.table("Results_from_Katarina/dexseq_exon_results.txt", header = T, stringsAsFactors = FALSE)

dexseq.ex <- dexseq.ex[!is.na(dexseq.ex$pvalue) & (duplicated(dexseq.ex$geneID, fromLast = T) | duplicated(dexseq.ex$geneID, fromLast = F) ) ,c("geneID", "pvalue")]

dexseq.ex.spl <- split(dexseq.ex, dexseq.ex$geneID)

dexseq.g <- data.frame(geneID=names(dexseq.ex.spl), pvalue=0, padjust=0)

dexseq.g$pvalue <- sapply(dexseq.ex.spl, function(g){
  
  n <- nrow(g)
  pv <- min(sort(g$pvalue,decreasing = FALSE)*n/(1:n))
  
})


dexseq.g$padjust <- p.adjust(dexseq.g$pvalue, method = "BH")


write.table(dexseq.g, "Results_from_Katarina/dexseq_gene_Simes_results.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



# comapre with DEXSeq results

dexseq.g.org <- read.table("Results_from_Katarina/dexseq_gene_results.txt", header = T, stringsAsFactors = FALSE)


d <- merge(dexseq.g.org, dexseq.g, by = 1, all=TRUE)

sum(is.na(d$padjust.x))
sum(is.na(d$padjust.y))
sum(is.na(d$padjust.y) && is.na(d$padjust.x))

d$padjust.x[is.na(d$padjust.x)] <- 1.1
d$padjust.y[is.na(d$padjust.y)] <- 1.1


pdf("Results_from_Katarina/dexseq_gene_results.pdf")
smoothScatter(d$padjust.x, d$padjust.y, nrpoints = Inf, pch=19, xlab="DEXSeq", ylab="Simes")
smoothScatter(log(d$padjust.x), log(d$padjust.y), nrpoints = Inf, xlab="DEXSeq", ylab="Simes")
dev.off()


#######################################################
# load information about simulation 
#######################################################

simu_info.e <- read.table("Simu_info/true_exons_simulation.txt", header=T)
simu_info.g <- read.table("Simu_info/true_genes_simulation.txt", header=T)


head(simu_info.g)
head(simu_info.e)

simu_info1 <- simu_info.g[!duplicated(simu_info.g$Gene),c("Gene", "status", "num", "rpk.mean")]

# count proportion of defferentially spliced exons
simu_info.e.spl <- split(simu_info.e, simu_info.e$Gene)

prop <- sapply(simu_info.e.spl, function(g){
  
  paste0(sum(g$status_exon==1), "/", sum(g$status_exon==0 | g$status_exon==1))

})


library(limma)
prop2 <- strsplit2(prop, "/")


simu_info2 <- data.frame(Gene = names(prop), num.ex=prop2[,2], num.diff.ex=prop2[,1])


simu_info <- merge(simu_info1, simu_info2, by = "Gene", all=TRUE)


write.table(simu_info, "Simu_info/simu_info.xls", quote = F, sep = "\t", row.names = F, col.names = T)


#######################################################
# merge Diff Expr results with simu_info
#######################################################

results <- list()

# results produced by Katarina

rt <- read.table("Results_from_Katarina/dexseq_gene_results.txt", header = T, stringsAsFactors = F)
head(rt)
rt <- rt[,c("Gene", "padjust")]
colnames(rt) <- c("Gene",  "adjPValue_htseq_dexseq")

results[["htseq_dexseq"]] <- rt



rt <- read.table("Results_from_Katarina/dexseq_gene_Simes_results.txt", header = T, stringsAsFactors = F)
head(rt)
rt <- rt[,c("geneID","pvalue" ,"padjust")]
colnames(rt) <- c("Gene", "PValue_htseq_dexseq_simons", "adjPValue_htseq_dexseq_simons")

results[["htseq_dexseq_simons"]] <- rt


rt <- read.table("Results_from_Katarina/voom_gene_results.txt", header = T, stringsAsFactors = F)
head(rt)
rt <- rt[,c("gene.id","P.Value" ,"FDR")]
colnames(rt) <- c("Gene", "PValue_fc_voomex", "adjPValue_fc_voomex")

results[["fc_voomex"]] <- rt


## problems when adding this results (same genes have multiple p-values)
# rt <- read.table("Results_from_Katarina/cuffdiff_results.txt", header = T, stringsAsFactors = F)
# head(rt)
# rt <- rt[,c("ensembl_gene_id","p_value" ,"q_value")]
# colnames(rt) <- c("Gene", "PValue_cuffdiff", "adjPValue_cuffdiff")

# results[["cuffdiff"]] <- rt


# DM results


rt <- read.table("DM/fc/DM_fc_results.xls", header = T, stringsAsFactors = F)
head(rt)
rt <- rt[,c("GeneID","PValue" ,"FDR")]
colnames(rt) <- c("Gene", "PValue_fc_DM", "adjPValue_fc_DM")

results[["fc_DM"]] <- rt



rt <- read.table("DM/htseq/DM_htseq_results.xls", header = T, stringsAsFactors = F)
head(rt)
rt <- rt[,c("GeneID","PValue" ,"FDR")]
colnames(rt) <- c("Gene", "PValue_htseq_DM", "adjPValue_htseq_DM")

results[["htseq_DM"]] <- rt



rt <- read.table("DM/miso/DM_miso_results.xls", header = T, stringsAsFactors = F)
head(rt)
rt <- rt[,c("GeneID","PValue" ,"FDR")]
colnames(rt) <- c("Gene", "PValue_miso_DM", "adjPValue_miso_DM")

results[["miso_DM"]] <- rt


## merge results into one table


table <- simu_info
for(i in 1:length(results)){
  table <- merge(table, results[[i]], by="Gene", all=TRUE)  
}

table <- unique(table)

write.table(table, "PLOTS/Table_all_results.xls",  quote = F, sep = "\t", row.names = F, col.names = T)

## check table 

dim(table[!complete.cases(table), ])



# # double check when adding cuffdiff results 
# table <- table[!is.na(table$Gene), ]
# table.dupl <- table[duplicated(table$Gene, fromLast = T) | duplicated(table$Gene, fromLast = F), ]
# write.table(table.dupl, "PLOTS/Table_all_results_dupl.xls",  quote = F, sep = "\t", row.names = F, col.names = T)
# table <- unique(table)
# table.dupl <- table[duplicated(table$Gene, fromLast = T) | duplicated(table$Gene, fromLast = F), ]
# write.table(table.dupl, "PLOTS/Table_all_results_dupl.xls",  quote = F, sep = "\t", row.names = F, col.names = T)



sum(!is.na(table$PValue_fc_voomex))
sum(!is.na(table$PValue_fc_DM))


          



#######################################################
# ->>>>>> load results
#######################################################


table <- read.table("PLOTS/Table_all_results.xls", header = T, stringsAsFactors = F)

colors <- c("hotpink", "magenta", "dodgerblue3",  "orange", "brown", "brown1")
names(colors) <- c("htseq_dexseq", "htseq_dexseq_simons", "fc_voomex",  "fc_DM", "htseq_DM", "miso_DM")                                  

name <- ""



#######################################################
# histograms of p-values
#######################################################


pdf(paste0("PLOTS/", name, "hist_pvalues.pdf"))


hist(table$PValue_fc_DM, col=colors["fc_DM"], breaks=50, main="fc_DM", cex.main=3, cex.lab=1.5, cex.axis = 1.5)
hist(table$PValue_htseq_DM, col=colors["htseq_DM"], breaks=50, main="htseq_DM", cex.main=3, cex.lab=1.5, cex.axis = 1.5)
hist(table$PValue_miso_DM, col=colors["miso_DM"], breaks=50, main="miso_DM", cex.main=3, cex.lab=1.5, cex.axis = 1.5)
hist(table$PValue_htseq_dexseq_simons, col=colors["htseq_dexseq_simons"], breaks=50, main="htseq_dexseq_simons", cex.main=3, cex.lab=1.5, cex.axis = 1.5)
hist(table$PValue_fc_voomex, col=colors["fc_voomex"], breaks=50, main="fc_voomex", cex.main=3, cex.lab=1.5, cex.axis = 1.5)


dev.off()


#######################################################
# ROCx plot from benchTools -- DOES NOT WORK!!!
#######################################################


# install.packages("/home/gosia/R/packages/benchTools_0.0-2.tar.gz", dependencies=TRUE, lib="/home/gosia/R/libraries/3.1.0")

library(benchTools)

# DOES NOT WORK!!!
# # plots with X
# 
# table.tmp <- table[complete.cases(table[,c("PValue_fc_DM","PValue_htseq_DM", "PValue_miso_DM", "adjPValue_fc_DM","adjPValue_htseq_DM", "adjPValue_miso_DM")]), ]
# 
# score <- as.matrix(1-table.tmp[,c("PValue_fc_DM","PValue_htseq_DM", "PValue_miso_DM")])
# label <- table.tmp$status
# score.X <- as.matrix(1-table.tmp[,c("adjPValue_fc_DM","adjPValue_htseq_DM", "adjPValue_miso_DM")])
# 
# r <- rocX(score, label,  score.X = score.X, threshold.X = 0.95, label.ordering = c(FALSE, TRUE))
# 
# pdf(paste0("PLOTS/",name,"roc_all_DM.pdf"))
# plot(r, col=c("orange", "orange", "orange"), lty= c( 3), cex.main=3, cex.lab=1.5, cex.axis = 2, lwd=3)
# legend("bottomright", c("fc_DM", "htseq_DM", "miso_DM"), col=c("orange", "orange", "orange"), lty=c(1, 2, 3), lwd=2.5)
# dev.off()



score <- as.matrix(1-table[,c("PValue_fc_DM","PValue_htseq_DM", "PValue_miso_DM")])
label <- table$status

r <- rocX(score, label)

pdf(paste0("PLOTS/",name,"roc_all_DM.pdf"))
plot(r, col=c("orange", "brown", "yellow"), cex.main=3, cex.lab=1.5, cex.axis = 2, lwd=3)
legend("bottomright", c("fc_DM", "htseq_DM", "miso_DM"), col=c("orange", "brown", "yellow"), lty=1, lwd=2.5)
dev.off()


score <- as.matrix(1-table[,c("PValue_fc_DM","PValue_htseq_DM", "PValue_miso_DM")])
score[is.na(score)] <- 0
label <- table$status

r <- rocX(score, label)

pdf(paste0("PLOTS/",name,"roc_all_DMp1s.pdf"))
plot(r, col=c("orange", "brown", "yellow"), cex.main=3, cex.lab=1.5, cex.axis = 2, lwd=3)
legend("bottomright", c("fc_DM", "htseq_DM", "miso_DM"), col=c("orange", "brown", "yellow"), lty=1, lwd=2.5)
dev.off()



#######################################################
# ROCx plot  --- my function
#######################################################


ROCx_version1 <- function(pvs, apvs, status, col=2, pch=4, cex=2, add=FALSE, xlim=c(0,1)){
  # status==1, means true DS, status!=1, means not DS
  
  xx <- sort(unique(c( 0 , pvs , 1 )))
  sens <- sapply(xx, function(x) mean(pvs[status==1] <= x , na.rm = T))
  spec <- sapply(xx, function(x) mean(pvs[status==0] > x , na.rm = T))
  
  if(add==FALSE){
    plot( 1-spec, sens, xlab="FPR", ylab="TPR", type="l", col= col, cex=cex, xlim=xlim)
  }else{   
    lines( 1-spec, sens, xlab="FPR", ylab="TPR", type="l", col= col, cex=cex, xlim=xlim)
  }
  
  TPR <- sum( status==1 & apvs < 0.05 , na.rm = T) / sum(status==1, na.rm = T)
  points((approx(sens, 1-spec, xout=TPR)$y), TPR, col=col, pch=pch, cex=cex)
  
  
}



ROCx_notNorm_version1 <- function(pvs, apvs, status, add=FALSE, ...){
  # status==1, means true DS, status!=1, means not DS
  
  xx <- sort(unique(c( 0 , pvs , 1 )))
  P <- sum(status==1)
  N <- sum(status==0)
  
  NAs <- sum(is.na(pvs) & status==0)
  
  sens <- sapply(xx, function(x) sum(pvs[status==1] < x , na.rm = T)) / P
  spec <- sapply(xx, function(x) sum(pvs[status==0] >= x , na.rm = T)) / N + NAs/N
  
  if(add==FALSE){
    plot( 1 - spec, sens, xlab="False positive rate", ylab="True positive rate", type="l", ...)
  }else{   
    lines( 1 - spec, sens, type="l", ...)
  }
  
  TPR <- sum(apvs[status==1] < 0.05 , na.rm = T) / P
  points((approx(sens, 1-spec, xout=TPR)$y), TPR, ...)  
  
}



####


ROCx <- function(pvs, apvs, status, add=FALSE, ...){
  # status==1, means true DS, status!=1, means not DS
  
  xx <- sort(unique(c( 0 , pvs , 1 )))
  sens <- sapply(xx, function(x) mean(pvs[status==1] <= x , na.rm = T))
  spec <- sapply(xx, function(x) mean(pvs[status==0] > x , na.rm = T))
  
  if(add==FALSE){
    plot( 1-spec, sens, xlab="False positive rate", ylab="True positive rate", type="l", ...)
  }else{   
    lines( 1-spec, sens, type="l", ...)
  }
  
  TPR <- sum( status==1 & apvs < 0.05 , na.rm = T) / sum(status==1, na.rm = T)
  points((approx(sens, 1-spec, xout=TPR)$y), TPR, ...)  
  
}




ROCx_notNorm <- function(pvs, apvs, status, add=FALSE, ...){
  # status==1, means true DS, status!=1, means not DS
  
  xx <- sort(unique(c( 0 , pvs , 1 )))
  P <- sum(status==1)
  N <- sum(status==0)
  
  NAs <- sum(is.na(pvs) & status==0)
  
  TPRv <- sapply(xx, function(x) sum(pvs[status==1] < x , na.rm = T)) / P
  FPRv <- sapply(xx, function(x) sum(pvs[status==0] < x , na.rm = T)) / N 
  
  if(add==FALSE){
    plot(FPRv, TPRv, xlab="False positive rate", ylab="True positive rate", type="l", ...)
  }else{   
    lines(FPRv, TPRv, type="l", ...)
  }
  
  TPR <- sum(apvs[status==1] < 0.05 , na.rm = T) / P
  points((approx(TPRv, FPRv, xout=TPR)$y), TPR, ...)  
  
}



#######################################################
# generate different ROCx plots 
#######################################################



pdf(paste0("PLOTS/",name,"rocX_all_DM.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx(pvs, apvs, status, col="orange", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table[,c("PValue_htseq_DM")]
apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
ROCx(pvs, apvs, status,  add=TRUE, col="brown", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table[,c("PValue_miso_DM")]
apvs <- table[,c("adjPValue_miso_DM")]
status <- table$status
ROCx(pvs, apvs, status,  add=TRUE, col="yellow", pch=4, xlim=c(0,1), lwd=2, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "miso_DM"), col=c("orange", "brown", "yellow"), lty=1, lwd=2.5)
dev.off()




table.tmp <- table[complete.cases(table[,c("PValue_fc_DM","PValue_htseq_DM", "PValue_miso_DM", "adjPValue_fc_DM","adjPValue_htseq_DM", "adjPValue_miso_DM")]), ]


pdf(paste0("PLOTS/",name,"rocX_all_DM_complete_cases.pdf"))

pvs <- table.tmp[,c("PValue_fc_DM")]
apvs <- table.tmp[,c("adjPValue_fc_DM")]
status <- table.tmp$status
ROCx(pvs, apvs, status, col="orange", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table.tmp[,c("PValue_htseq_DM")]
apvs <- table.tmp[,c("adjPValue_htseq_DM")]
status <- table.tmp$status
ROCx(pvs, apvs, status,  add=TRUE, col="brown", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table.tmp[,c("PValue_miso_DM")]
apvs <- table.tmp[,c("adjPValue_miso_DM")]
status <- table.tmp$status
ROCx(pvs, apvs, status,  add=TRUE, col="yellow", pch=4, xlim=c(0,1), lwd=2, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "miso_DM"), col=c("orange", "brown", "yellow"), lty=1, lwd=2.5)
dev.off()




table.tmp <- table

pdf(paste0("PLOTS/",name,"rocX_all_DM_noNAs.pdf"))

pvs <- table.tmp[,c("PValue_fc_DM")]
apvs <- table.tmp[,c("adjPValue_fc_DM")]
pvs[is.na(pvs)] <- 1
apvs[is.na(apvs)] <- 1
status <- table.tmp$status

ROCx(pvs, apvs, status, col="orange", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table.tmp[,c("PValue_htseq_DM")]
apvs <- table.tmp[,c("adjPValue_htseq_DM")]
pvs[is.na(pvs)] <- 1
apvs[is.na(apvs)] <- 1
status <- table.tmp$status
ROCx(pvs, apvs, status,  add=TRUE, col="brown", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table.tmp[,c("PValue_miso_DM")]
apvs <- table.tmp[,c("adjPValue_miso_DM")]
pvs[is.na(pvs)] <- 1
apvs[is.na(apvs)] <- 1
status <- table.tmp$status
ROCx(pvs, apvs, status,  add=TRUE, col="yellow", pch=4, xlim=c(0,1), lwd=2, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "miso_DM"), col=c("orange", "brown", "yellow"), lty=1, lwd=2.5)
dev.off()





pdf(paste0("PLOTS/",name,"rocX_all_DM_notNorm.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col="orange", pch=4, xlim=c(0,1), ylim=c(0,1), lwd=2, cex=3)

pvs <- table[,c("PValue_htseq_DM")]
apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col="brown", pch=4,  lwd=2, cex=3)

pvs <- table[,c("PValue_miso_DM")]
apvs <- table[,c("adjPValue_miso_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col="yellow", pch=4, lwd=2, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "miso_DM"), col=c("orange", "brown", "yellow"), lty=1, lwd=2.5)
dev.off()





pdf(paste0("PLOTS/",name,"rocX_all.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx(pvs, apvs, status, col="orange", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table[,c("PValue_htseq_DM")]
apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
ROCx(pvs, apvs, status,  add=TRUE, col="brown", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table[,c("PValue_fc_voomex")]
apvs <- table[,c("adjPValue_fc_voomex")]
status <- table$status
ROCx(pvs, apvs, status,  add=TRUE, col="blue", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table[,c("PValue_htseq_dexseq_simons")]
apvs <- table[,c("adjPValue_htseq_dexseq_simons")]
status <- table$status
ROCx(pvs, apvs, status,  add=TRUE, col="magenta", pch=4, xlim=c(0,1), lwd=2, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons"), col=c("orange", "brown", "blue", "magenta"), lty=1, lwd=2.5)
dev.off()





pdf(paste0("PLOTS/",name,"rocX_all_notNorm.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col="orange", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table[,c("PValue_htseq_DM")]
apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col="brown", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table[,c("PValue_fc_voomex")]
apvs <- table[,c("adjPValue_fc_voomex")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col="blue", pch=4, xlim=c(0,1), lwd=2, cex=3)

pvs <- table[,c("PValue_htseq_dexseq_simons")]
apvs <- table[,c("adjPValue_htseq_dexseq_simons")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col="magenta", pch=4, xlim=c(0,1), lwd=2, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons"), col=c("orange", "brown", "blue", "magenta"), lty=1, lwd=2.5)
dev.off()



#######################################################
# generate FINAL ROCx plots 
#######################################################



pdf(paste0("PLOTS/",name,"rocX_all_DM.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx(pvs, apvs, status, col=colors["fc_DM"], pch=4, xlim=c(0,1), lwd=3, cex=3, cex.lab=1.5, cex.axis=1.5)

pvs <- table[,c("PValue_htseq_DM")]
apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
ROCx(pvs, apvs, status,  add=TRUE, col=colors["htseq_DM"], pch=4, xlim=c(0,1), lwd=3, cex=3)

pvs <- table[,c("PValue_miso_DM")]
apvs <- table[,c("adjPValue_miso_DM")]
status <- table$status
ROCx(pvs, apvs, status,  add=TRUE, col=colors["miso_DM"], pch=4, xlim=c(0,1), lwd=3, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "miso_DM"), col=colors[c("fc_DM", "htseq_DM", "miso_DM")], lty=1, lwd=3, cex=1.5)
dev.off()





pdf(paste0("PLOTS/",name,"rocX_all_DM_notNorm.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col=colors["fc_DM"], pch=4, xlim=c(0,1), ylim=c(0,1), lwd=3, cex=3, cex.lab=1.5, cex.axis=1.5)

pvs <- table[,c("PValue_htseq_DM")]
apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["htseq_DM"], pch=4,  lwd=3, cex=3)

pvs <- table[,c("PValue_miso_DM")]
apvs <- table[,c("adjPValue_miso_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["miso_DM"], pch=4, lwd=3, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "miso_DM"), col=colors[c("fc_DM", "htseq_DM", "miso_DM")], lty=1, lwd=3, cex=1.5)
dev.off()





pdf(paste0("PLOTS/",name,"rocX_all.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx(pvs, apvs, status, col=colors["fc_DM"], pch=4, xlim=c(0,1), lwd=3, cex=3, cex.lab=1.5, cex.axis=1.5)

pvs <- table[,c("PValue_htseq_DM")]
apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
ROCx(pvs, apvs, status,  add=TRUE, col=colors["htseq_DM"], pch=4, xlim=c(0,1), lwd=3, cex=3)

pvs <- table[,c("PValue_fc_voomex")]
apvs <- table[,c("adjPValue_fc_voomex")]
status <- table$status
ROCx(pvs, apvs, status,  add=TRUE, col=colors["fc_voomex"], pch=4, xlim=c(0,1), lwd=3, cex=3)

pvs <- table[,c("PValue_htseq_dexseq_simons")]
apvs <- table[,c("adjPValue_htseq_dexseq_simons")]
status <- table$status
ROCx(pvs, apvs, status,  add=TRUE, col=colors["htseq_dexseq_simons"], pch=4, xlim=c(0,1), lwd=3, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons"), col=colors[c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons")], lty=1, lwd=3, cex=1.5)
dev.off()





pdf(paste0("PLOTS/",name,"rocX_all_notNorm.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col=colors["fc_DM"], pch=4, xlim=c(0,1), lwd=4, cex=3, cex.lab=1.5, cex.axis=1.5)

pvs <- table[,c("PValue_htseq_DM")]
apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["htseq_DM"], pch=4, xlim=c(0,1), lwd=4, cex=3)

pvs <- table[,c("PValue_fc_voomex")]
apvs <- table[,c("adjPValue_fc_voomex")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["fc_voomex"], pch=4, xlim=c(0,1), lwd=4, cex=3)

pvs <- table[,c("PValue_htseq_dexseq_simons")]
apvs <- table[,c("adjPValue_htseq_dexseq_simons")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["htseq_dexseq_simons"], pch=4, xlim=c(0,1), lwd=4, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons"), col=colors[c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons")], lty=1, lwd=4, cex=1.5)
dev.off()





#######################################################
# TPR vs acieved FDR
#######################################################



TPRvsFDR <- function(status, apvs, FDR.cut.off=c(0.01, 0.05, 0.1), pch=c( 22, 23, 24), col="red", cex.axis=1.5, cex=2.5, add=FALSE, ...){
  
  apvs[is.na(apvs)] <- 1
  
  n <- length(status)
  q <- length(FDR.cut.off)
  TPR <- rep(0, q)
  FDR <- rep(0, q)
  
  for(i in 1:q){
    # i=1
    status.est <- as.numeric(apvs <= FDR.cut.off[i])
    
    TP <- sum(status==1 & status.est==1)
    FP <- sum(status==0 & status.est==1)
    FN <- sum(status==1 & status.est==0)
    
    TPR[i] <- TP/(TP+FN)
    FDR[i] <- FP/(FP+TP)
    
  }  
  tf <- cbind( FDR , TPR)
  
  bg <- rep(col,q)
  bg[FDR > FDR.cut.off] <- "white"
  
  if(add==FALSE){
   
    plot(1:2, type="n", xlab="Achieved FDR", ylab="TPR", xaxt="n", col=col, cex.axis=cex.axis, ...)
    axis(side=1, at=(2:10)/10, labels=(2:10)/10, las=1, col.ticks="grey", col.axis="grey", cex.axis=cex.axis)
    axis(side=1, at=FDR.cut.off, labels=FDR.cut.off, las=1, cex.axis=cex.axis)

    for(i in 1:q)
    lines(rep(FDR.cut.off[i], 50), seq(-0.1,1.1,length.out=50), type="b", pch=pch[i], cex=0.5, bg=1)
    
    lines(tf, type="l", col=col, ...)
    points(tf, pch=pch, bg=bg, col=col, cex=cex, ...)
    
  }else{    
    lines(tf, type="l", col=col, ...)
    points(tf, pch=pch, bg=bg, col=col, cex=cex,...)
  }
  
  
  
}




pdf(paste0("PLOTS/", name, "TPRachievedFDR.pdf"))

apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["fc_DM"], xlim=c(0,0.3), ylim=c(0,1), lwd=4, cex=2.5, cex.lab=1.5, cex.axis=1.5)

apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["htseq_DM"], add=TRUE, lwd=4, cex=2.5)

apvs <- table[,c("adjPValue_fc_voomex")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["fc_voomex"], add=TRUE, lwd=4, cex=2.5)

apvs <- table[,c("adjPValue_htseq_dexseq_simons")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["htseq_dexseq_simons"], add=TRUE, lwd=4, cex=2.5)

legend("topright", c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons"), col=colors[ c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons")], lty=1, lwd=4, cex=1.5)
dev.off()



pdf(paste0("PLOTS/", name, "TPRachievedFDR2.pdf"))

apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["fc_DM"], xlim=c(0,0.3), ylim=c(0,1), lwd=3, cex=2.5, cex.lab=1.5, cex.axis=1.5)

apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["htseq_DM"], add=TRUE, lwd=3, cex=2.5)

apvs <- table[,c("adjPValue_fc_voomex")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["fc_voomex"], add=TRUE, lwd=3, cex=2.5)

apvs <- table[,c("adjPValue_htseq_dexseq_simons")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["htseq_dexseq_simons"], add=TRUE, lwd=3, cex=2.5)

apvs <- table[,c("adjPValue_htseq_dexseq")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["htseq_dexseq"], add=TRUE, lwd=3, cex=2.5)


legend("topright", c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons", "htseq_dexseq"), col=colors[ c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons", "htseq_dexseq")], lty=1, lwd=3, cex=1.5)
dev.off()



pdf(paste0("PLOTS/", name, "TPRachievedFDR_DM.pdf"))

apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["fc_DM"], xlim=c(0,0.25), ylim=c(0,1), lwd=3, cex=2.5, cex.lab=1.5, cex.axis=1.5)

apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["htseq_DM"], add=TRUE, lwd=3, cex=2.5)

apvs <- table[,c("adjPValue_miso_DM")]
status <- table$status
TPRvsFDR(status, apvs, col=colors["miso_DM"], add=TRUE, lwd=3, cex=2.5)


legend("topright", c("fc_DM", "htseq_DM", "miso_DM"), col=colors[ c("fc_DM", "htseq_DM", "miso_DM")], lty=1, lwd=3, cex=1.5)
dev.off()



#######################################################
# venn diagrams
#######################################################
library(VennDiagram)

FDR.cutoff <- 0.05


plot.venn <- function(venne.genes, colors, venn.methods, name="", cex=1, cat.cex=1.4, lwd=2, lty=1, alpha=0.5, margin=0.05){
  
  colors.col <- colors[venn.methods]
  
  if(length(venn.methods)==4)
    colors.col <- colors.col[c(1, 3, 4 ,2)]
  
  venn.d <- venn.diagram(venne.genes[venn.methods], filename=NULL, fill = colors[venn.methods], col=colors.col, cex=cex, cat.cex=cat.cex, lwd=lwd, lty=lty, alpha=alpha, margin=margin)
  
  pdf(paste0("PLOTS/", name, "Venn_Diagr.pdf"))
  grid.draw(venn.d)
  dev.off()
  
  
}


## TP as a separate circle

venne.genes <- list()

venne.genes$fc_DM <- na.omit(table[ table$adjPValue_fc_DM < FDR.cutoff, "Gene"])
venne.genes$htseq_DM <- na.omit(table[ table$adjPValue_htseq_DM < FDR.cutoff, "Gene"])
venne.genes$miso_DM <- na.omit(table[ table$adjPValue_miso_DM < FDR.cutoff, "Gene"])
venne.genes$fc_voomex <- na.omit(table[ table$adjPValue_fc_voomex < FDR.cutoff, "Gene"])
venne.genes$htseq_dexseq_simons <- na.omit(table[ table$adjPValue_htseq_dexseq_simons < FDR.cutoff, "Gene"])
venne.genes$htseq_dexseq <- na.omit(table[ table$adjPValue_htseq_dexseq < FDR.cutoff, "Gene"])
venne.genes$True <- table[table$status == 1, "Gene"]



plot.venn(venne.genes, colors=c(colors[c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons")], True="grey"), venn.methods = c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simons", "True"), name=paste0(name, "All_"), margin=0.1, cat.cex=1.4, cex=1.7)


plot.venn(venne.genes, colors=c(colors[c("fc_DM", "htseq_DM", "htseq_dexseq", "htseq_dexseq_simons")], True="grey"), venn.methods = c("fc_DM", "htseq_DM", "htseq_dexseq", "htseq_dexseq_simons", "True"), name=paste0(name, "All2_"), margin=0.1, cat.cex=1.4, cex=1.7)



plot.venn(venne.genes, colors=c(colors[c("fc_DM", "htseq_DM", "miso_DM")], True="grey"), venn.methods = c("fc_DM", "htseq_DM", "miso_DM", "True"), name=paste0(name, "All_DM_"), margin=0.1, cat.cex=1.4, cex=1.7)




























































































