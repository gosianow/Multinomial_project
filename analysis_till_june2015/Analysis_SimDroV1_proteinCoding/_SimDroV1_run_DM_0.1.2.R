##############################################################################

# BioC 3.0
# Created 28 Jan 2015:

# Run DM_0.1.2 on 
# - htseq couts
# - BitSeq counts

# Modified 5 Feb 2015
# - Use source files that will be wrapped in DM_0.1.3

##############################################################################


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding")

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata



out.dir <- "DM_0.1.2/"
dir.create(out.dir, showWarnings=F, recursive=T)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)


### Source all R files in DM package / DM_0.1.3
Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])


##############################################################################################################
# htseq counts
##############################################################################################################

out.dir <- "DM_0.1.2/htseq/"
dir.create(out.dir, showWarnings=F, recursive=T)


#### load htseq counts
# htseq <- list()
# 
# for(i in 1:6){
#   # i = 1
#   htseq[[i]] <- read.table(paste0("htseq/htseq", i, ".txt"), header = FALSE, as.is = TRUE)
#   colnames(htseq[[i]]) <- c("exon", paste0("counts", i))  
# }
# 
# htseq.m <- htseq[[1]]
# 
# for(i in 2:6){  
#   htseq.m <- merge(htseq.m, htseq[[i]][,c("exon", paste0("counts", i))], by = "exon", all = TRUE, sort = FALSE)
# }
# 
# write.table(htseq.m, paste0("htseq/htseq_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


htseq <- read.table("htseq/htseq_counts.txt", header = TRUE, as.is = TRUE)
htseq <- htseq[grepl(pattern = "FB", htseq$exon),]

expr <- htseq[,2:7]
colnames(expr) <- metadata$SampleName1
gene_id <- strsplit2(htseq[,1], ":")[,1]
ete_id <- htseq[,1]

dgeOrg <- DGEList(counts=expr, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

name1 <- "htseq"
dge <- dgeOrg
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)


######################################
# htseq counts: filtering
######################################

# out.dir <- "DM_0.1.2/htseq/Filtering005/"
# dir.create(out.dir, showWarnings=F, recursive=T)
# 
# min.samps <- 3
# min.transcript.prop <- 0.05 # in cpm
# min.gene.exp <- 1 # in cpm

out.dir <- "DM_0.1.2/htseq/Filtering0/"
dir.create(out.dir, showWarnings=F, recursive=T)

min.samps <- 3
min.transcript.prop <- 0 # in cpm
min.gene.exp <- 1 # in cpm


seq.depth <- colSums(expr)

dge <- dgeOrg
expr.cpm <- cpm(dge)
rownames(expr.cpm) <- dge$genes$ete_id

expr.cpm.spl <- split(data.frame(expr.cpm), dge$genes$gene_id) 


trans.to.keep <- lapply(seq(length(expr.cpm.spl)), function(g){
  # g = 1
  
  ### no genes with one transcript
  if(dim(expr.cpm.spl[[g]])[1] == 1)
    return(NULL)
  
  ### genes with min expression in all samples
  if(any(colSums(expr.cpm.spl[[g]]) < min.gene.exp))
    return(NULL)
  
  ### transcripts with min expression in samples
  # trans.to.keep = apply(expr.cpm.spl[[g]], 1, function(t) sum(t > min.transcript.exp) >= 10 )
  
  tot <- colSums(expr.cpm.spl[[g]])
  
  trans.to.keep = apply(expr.cpm.spl[[g]], 1, function(t){
    # t = expr.cpm.spl[[g]][8, ]    
    r <- t / tot
    sum(r > min.transcript.prop) >= min.samps 
    
  } )
  
  ### no genes with one transcript
  if(sum(trans.to.keep) <= 1)
    return(NULL)
  
  return(names(which(trans.to.keep == TRUE)))
  
})

trans.to.keep <- unlist(trans.to.keep)

length(trans.to.keep)


keep <- dge$genes$ete_id %in% trans.to.keep
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
rownames(dge$counts) <- dge$genes$ete_id
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)


count.method <- "htseq"
save(dge, file = paste0(out.dir, "/dge_counts_",count.method,".RData"))
load(paste0(out.dir, "/dge_counts_",count.method,".RData"))

######################################
# htseq counts: filtering from DEXSeq
######################################

out.dir <- "DM_0.1.2/htseq/Filtering_DEXSeq/"
dir.create(out.dir, showWarnings=F, recursive=T)


### DEXSeq exon results
rt <- read.table("Results_from_Katarina/dexseq_1.10.8_exon_htseq_results.txt", header = T, as.is = TRUE, sep = "\t", fill = TRUE )
head(rt)

rt <- rt[, c("groupID", "dispersion", "pvalue", "padj")]
rt <- rt[complete.cases(rt), ]

keep <- gsub(pattern = "E", replacement = "", rt$groupID)
length(keep)

dge <- dgeOrg

dge <- dge[dge$genes$ete_id %in% keep,]

dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
rownames(dge$counts) <- dge$genes$ete_id
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)


count.method <- "htseq"
save(dge, file = paste0(out.dir, "/dge_counts_",count.method,".RData"))


##############################################################################################################
# bitseq counts
##############################################################################################################

out.dir <- "DM_0.1.2/bitseq_1.10.0/"
dir.create(out.dir, showWarnings=F, recursive=T)


bitseq <- read.table("BitSeq_1.10.0/BitSeq_counts.txt", header = TRUE, as.is = TRUE)
head(bitseq)


expr <- bitseq[,3:8]
colnames(expr) <- metadata$SampleName1
gene_id <- bitseq[,"gene"]
ete_id <- paste0(bitseq[, "gene"],":", bitseq[, "transcript"])

dgeOrg <- DGEList(counts=expr, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

name1 <- "bitseq"
dge <- dgeOrg
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)


######################################
# bitseq counts: filtering
######################################

# out.dir <- "DM_0.1.2/bitseq_1.10.0/Filtering005/"
# dir.create(out.dir, showWarnings=F, recursive=T)
# 
# min.samps <- 3
# min.transcript.prop <- 0.05 # in cpm
# min.gene.exp <- 1 # in cpm


out.dir <- "DM_0.1.2/bitseq_1.10.0/Filtering0/"
dir.create(out.dir, showWarnings=F, recursive=T)

min.samps <- 3
min.transcript.prop <- 0 # in cpm
min.gene.exp <- 1 # in cpm


seq.depth <- colSums(expr)

dge <- dgeOrg
expr.cpm <- cpm(dge)
rownames(expr.cpm) <- dge$genes$ete_id

expr.cpm.spl <- split(data.frame(expr.cpm), dge$genes$gene_id) 


trans.to.keep <- lapply(seq(length(expr.cpm.spl)), function(g){
  # g = 1
  
  ### no genes with one transcript
  if(dim(expr.cpm.spl[[g]])[1] == 1)
    return(NULL)
  
  ### genes with min expression in all samples
  if(any(colSums(expr.cpm.spl[[g]]) < min.gene.exp))
    return(NULL)
  
  ### transcripts with min expression in samples
  # trans.to.keep = apply(expr.cpm.spl[[g]], 1, function(t) sum(t > min.transcript.exp) >= 10 )
  
  tot <- colSums(expr.cpm.spl[[g]])
  
  trans.to.keep = apply(expr.cpm.spl[[g]], 1, function(t){
    # t = expr.cpm.spl[[g]][8, ]    
    r <- t / tot
    sum(r > min.transcript.prop) >= min.samps 
    
  } )
  
  ### no genes with one transcript
  if(sum(trans.to.keep) <= 1)
    return(NULL)
  
  return(names(which(trans.to.keep == TRUE)))
  
})

trans.to.keep <- unlist(trans.to.keep)

length(trans.to.keep)


keep <- dge$genes$ete_id %in% trans.to.keep
dge <- dge[keep,]
dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
rownames(dge$counts) <- dge$genes$ete_id
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)


count.method <- "bitseq"
save(dge, file = paste0(out.dir, "/dge_counts_",count.method,".RData"))

######################################
# bitseq counts: filtering from DEXSeq
######################################

out.dir <- "DM_0.1.2/bitseq_1.10.0/Filtering_DEXSeq/"
dir.create(out.dir, showWarnings=F, recursive=T)


### DEXSeq exon results
rt <- read.table("DEXSeq_1.10.8/bitseq_1.10.0/dexseq_1.10.8_exon_bitseq_results.txt", header = T, as.is = TRUE, sep = "\t", fill = TRUE )
head(rt)

rt <- rt[, c("groupID", "dispersion", "pvalue", "padj")]
rt <- rt[complete.cases(rt), ]

keep <- gsub(pattern = "E", replacement = "", rt$groupID)
length(keep)

dge <- dgeOrg

dge <- dge[dge$genes$ete_id %in% keep,]

dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
rownames(dge$counts) <- dge$genes$ete_id
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

nlevels(dge$genes$gene_id)


count.method <- "bitseq"
save(dge, file = paste0(out.dir, "/dge_counts_",count.method,".RData"))



######################################
### DM_C
######################################

name1 <- paste0(count.method, "_DM_C")
mcCores <- 20


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, subset=Inf, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)


write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)


dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))


gc()

dgeDM.CommonDisp <- dgeDM


######################################
### DM_TG:optim
######################################

name1 <- paste0(count.method, "_DM_TG-optim")
mcCores <- 20


## run DM pipeline
# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

dgeDM <- dgeDM.CommonDisp

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)


dgeDM <- dmEstimateTagwiseDisp(dgeDM, group = NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[2], interval = c(0, 1e+10), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 15, gridRange = c(-7, 7), trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf = 4, span = 0.3, mcCores = mcCores, verbose = FALSE)


write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




load(paste0(out.dir, "/",name1,"_dgeDM.RData"))


pdf(paste0(out.dir, "/",name1,"_DispVsMean.pdf"))

df <- data.frame(meanExpr = log10(dgeDM$meanExpr+1), tagwiseDispersion = log10(dgeDM$tagwiseDispersion))
ggp <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion)) +
  theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold")) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)
print(ggp)

dev.off()


pdf(paste0(out.dir, "/",name1,"_DispVsMean_df.pdf"), width = 7, height = 5)

rownames(dgeDM$table) <- dgeDM$table$GeneID
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = dgeDM$table[names(dgeDM$meanExpr), "df"] )) +
  theme_bw() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)+
  # theme(legend.position="none")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
  #scale_colour_gradientn(colours = myPalette(100), limits=c(1, 15), name = "df")
  scale_colour_gradient(limits=c(1, 10), low = "blue", high="red", name = "df", na.value = "red")
print(ggp2)

dev.off()


gc()


######################################
### DM_TG:constrOptim
######################################

name1 <- paste0(count.method, "_DM_TG-constrOptim")
mcCores <- 20


## run DM pipeline
# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

dgeDM <- dgeDM.CommonDisp

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)


dgeDM <- dmEstimateTagwiseDisp(dgeDM, group = NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[3], interval = c(0, 1e+10), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 15, gridRange = c(-7, 7), trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf = 4, span = 0.3, mcCores = mcCores, verbose = FALSE)


write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))



load(paste0(out.dir, "/",name1,"_dgeDM.RData"))

pdf(paste0(out.dir, "/",name1,"_DispVsMean.pdf"))

df <- data.frame(meanExpr = log10(dgeDM$meanExpr+1), tagwiseDispersion = log10(dgeDM$tagwiseDispersion))
ggp <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion)) +
  theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold")) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)
print(ggp)

dev.off()


pdf(paste0(out.dir, "/",name1,"_DispVsMean_df.pdf"), width = 7, height = 5)

rownames(dgeDM$table) <- dgeDM$table$GeneID
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = dgeDM$table[names(dgeDM$meanExpr), "df"] )) +
  theme_bw() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)+
  # theme(legend.position="none")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
  #scale_colour_gradientn(colours = myPalette(100), limits=c(1, 15), name = "df")
  scale_colour_gradient(limits=c(1, 10), low = "blue", high="red", name = "df", na.value = "red")
print(ggp2)

dev.off()


######################################
### DM_TG:optimize
######################################

name1 <- paste0(count.method, "_DM_TG-optimize")
mcCores <- 20


## run DM pipeline
# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

dgeDM <- dgeDM.CommonDisp

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)


dgeDM <- dmEstimateTagwiseDisp(dgeDM, group = NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[1], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 15, gridRange = c(-7, 7), trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf = 4, span = 0.3, mcCores = mcCores, verbose = FALSE)


write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))



load(paste0(out.dir, "/",name1,"_dgeDM.RData"))

pdf(paste0(out.dir, "/",name1,"_DispVsMean.pdf"))

df <- data.frame(meanExpr = log10(dgeDM$meanExpr+1), tagwiseDispersion = log10(dgeDM$tagwiseDispersion))
ggp <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion)) +
  theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold")) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)
print(ggp)

dev.off()


pdf(paste0(out.dir, "/",name1,"_DispVsMean_df.pdf"), width = 7, height = 5)

rownames(dgeDM$table) <- dgeDM$table$GeneID
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = dgeDM$table[names(dgeDM$meanExpr), "df"] )) +
  theme_bw() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)+
  # theme(legend.position="none")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
  #scale_colour_gradientn(colours = myPalette(100), limits=c(1, 15), name = "df")
  scale_colour_gradient(limits=c(1, 10), low = "blue", high="red", name = "df", na.value = "red")
print(ggp2)

dev.off()


######################################
### DM_TG:gridCommon
######################################

name1 <- paste0(count.method, "_DM_TG-gridCommon")
mcCores <- 20


## run DM pipeline
# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

dgeDM <- dgeDM.CommonDisp

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)


dgeDM <- dmEstimateTagwiseDisp(dgeDM, group = NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 15, gridRange = c(-7, 7), trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf = 4, span = 0.3, mcCores = mcCores, verbose = FALSE, plot = FALSE)


write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)


dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))

gc()


load(paste0(out.dir, "/",name1,"_dgeDM.RData"))


pdf(paste0(out.dir, "/",name1,"_DispVsMean.pdf"))

df <- data.frame(meanExpr = log10(dgeDM$meanExpr+1), tagwiseDispersion = log10(dgeDM$tagwiseDispersion))
ggp <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion)) +
  theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold")) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)
print(ggp)

dev.off()


pdf(paste0(out.dir, "/",name1,"_DispVsMean_df.pdf"), width = 7, height = 5)

rownames(dgeDM$table) <- dgeDM$table$GeneID
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = dgeDM$table[names(dgeDM$meanExpr), "df"] )) +
  theme_bw() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)+
  # theme(legend.position="none")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
  #scale_colour_gradientn(colours = myPalette(100), limits=c(1, 15), name = "df")
  scale_colour_gradient(limits=c(1, 10), low = "blue", high="red", name = "df", na.value = "red")
print(ggp2)

dev.off()



######################################
### DM_TG:gridNone
######################################

name1 <- paste0(count.method, "_DM_TG-gridNone20")
mcCores <- 20


## run DM pipeline
# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

dgeDM <- dgeDM.CommonDisp


write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)


dgeDM <- dmEstimateTagwiseDisp(dgeDM, group = NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 31, gridRange = c(-10, 20), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 4, span = 0.3, mcCores = mcCores, verbose = FALSE)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)


dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))



load(paste0(out.dir, "/",name1,"_dgeDM.RData"))


pdf(paste0(out.dir, "/",name1,"_DispVsMean.pdf"))

df <- data.frame(meanExpr = log10(dgeDM$meanExpr+1), tagwiseDispersion = log10(dgeDM$tagwiseDispersion))
ggp <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion)) +
  theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold")) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)
print(ggp)

dev.off()


pdf(paste0(out.dir, "/",name1,"_DispVsMean_df.pdf"), width = 7, height = 5)

rownames(dgeDM$table) <- dgeDM$table$GeneID
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = dgeDM$table[names(dgeDM$meanExpr), "df"] )) +
  theme_bw() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)+
  # theme(legend.position="none")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
  #scale_colour_gradientn(colours = myPalette(100), limits=c(1, 15), name = "df")
  scale_colour_gradient(limits=c(1, 10), low = "blue", high="red", name = "df", na.value = "red")
print(ggp2)

dev.off()



######################################
### DM_TG:gridTrend
######################################

name1 <- paste0(count.method, "_DM_TG-gridTrend20-span03")
mcCores <- 20


## run DM pipeline
# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)

dgeDM <- dgeDM.CommonDisp


write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)


dgeDM <- dmEstimateTagwiseDisp(dgeDM, group = NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 31, gridRange = c(-10, 20), trend = c("none", "commonDispersion", "trendedDispersion")[3], priorDf = 4, span = 0.3, mcCores = mcCores, verbose = FALSE)

write.table(dgeDM$tagwiseDispersion, paste0(out.dir, "/",name1,"_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "tagwiseDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




load(paste0(out.dir, "/",name1,"_dgeDM.RData"))


pdf(paste0(out.dir, "/",name1,"_DispVsMean.pdf"))

df <- data.frame(meanExpr = log10(dgeDM$meanExpr+1), tagwiseDispersion = log10(dgeDM$tagwiseDispersion))
ggp <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion)) +
  theme_bw() +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold")) +
  geom_point(size = 1) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)
print(ggp)

dev.off()


pdf(paste0(out.dir, "/",name1,"_DispVsMean_df.pdf"), width = 7, height = 5)

rownames(dgeDM$table) <- dgeDM$table$GeneID
#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
myPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = dgeDM$table[names(dgeDM$meanExpr), "df"] )) +
  theme_bw() +
  geom_point(size = 2) +
  geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)+
  # theme(legend.position="none")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
  #scale_colour_gradientn(colours = myPalette(100), limits=c(1, 15), name = "df")
  scale_colour_gradient(limits=c(1, 10), low = "blue", high="red", name = "df", na.value = "red")
print(ggp2)

dev.off()


gc()


























































