######################################################
# BioC 3.0

# Created 11 June 2015 
# Analysis of Kim_adenocarcinoma data with DM_0.1.4 using htseq counts

######################################################


setwd("/home/Shared/data/seq/Kim_adenocarcinoma/")

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM_0.1.4/R/", full.names=TRUE)
for(i in Rfiles) source(i)


##########################################################################
# load metadata
##########################################################################


metadata <- read.table("metadata/Malgorzata Nowicka2014-11-04GSE37764.csv", stringsAsFactors=F, sep=",", header=T) 

metadata <- metadata[metadata$X == "RNA-seq",]

metadata

##############################################################################################################
# htseq counts
##############################################################################################################

count.method <- "htseq"

out.dir <- paste0("DM_0_1_4/",count.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


# ### load htseq counts
# htseq <- list()
# 
# for(i in metadata$ids){
#   # i = 1
#   htseq[[i]] <- read.table(paste0("DEXSeq_1.10.8/Exon_counts/", i, ".counts"), header = FALSE, as.is = TRUE)
#   colnames(htseq[[i]]) <- c("exon", i)  
# }
# 
# htseq.m <- htseq[[1]]
# 
# for(i in metadata$ids[-1]){  
#   htseq.m <- merge(htseq.m, htseq[[i]], by = "exon", all = TRUE, sort = FALSE)
# }
# 
# write.table(htseq.m, paste0("DEXSeq_1.10.8/Exon_counts/htseq_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


htseq <- read.table("DEXSeq_1.10.8/Exon_counts/htseq_counts.txt", header = TRUE, as.is = TRUE)
htseq <- htseq[grepl(pattern = "ENSG", htseq$exon),]

expr <- htseq[,metadata$ids]
gene_id <- strsplit2(htseq[,1], ":")[,1]
ete_id <- htseq[,1]

dgeOrg <- DGEList(counts=expr, group = metadata$Tissue.Type, genes=data.frame(gene_id=gene_id, ete_id=ete_id))

dge <- dgeOrg

write.table(data.frame(dge$genes, dge$counts), paste0(out.dir, "/dge_counts_", count.method ,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)





######################################
# htseq counts: filtering from DEXSeq
######################################

filter.method <- "Filtering_DEXSeq"

out.dir <- paste0("DM_0_1_4/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


### DEXSeq exon results
rt <- read.table("DEXSeq_1.10.8/diff_out_Full_tumorVSnormal/DEXSeq_exon_results.txt", header = T, as.is = TRUE, sep = "\t", fill = TRUE )
head(rt)

rt <- rt[, c("groupID", "featureID","dispersion", "pvalue", "padj")]
rt <- rt[complete.cases(rt), ]

### make the ete id sane as in dge
keep <- paste0(rt$groupID, ":", gsub(pattern = "E", replacement = "", rt$featureID))
head(keep)
length(keep)

dge <- dgeOrg

dge <- dge[dge$genes$ete_id %in% keep,]

dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting
rownames(dge$counts) <- dge$genes$ete_id
dge$counts <- split.data.frame(dge$counts, factor(dge$genes$gene_id, levels = unique(dge$genes$gene_id)))

nlevels(dge$genes$gene_id)

save(dge, file = paste0(out.dir, "/dge_counts_",count.method,".RData"))

pdf(paste0(out.dir, "/Hist_numberOfTranscripts.pdf"))
tt <- table(dge$genes$gene_id)
hist(tt, breaks = max(tt), col = "darkseagreen2", main = paste0(length(tt), " genes \n ", sum(tt) , " transcripts "), xlab = "Number of transcripts per gene")
dev.off()




##############################################

### run DM with different modes 

##############################################


count.methodList <- c("htseq", "bitseq")
count.method <- count.methodList[1]

filter.methodList <- c("Filtering_DEXSeq", "Filtering0", "Filtering001") 
filter.method <- filter.methodList[1]

out.dir <- paste0("DM_0_1_4/",count.method,"/",filter.method,"/")

load(paste0(out.dir, "/dge_counts_",count.method,".RData"))

modeList = c("constrOptim2G")
mode <- modeList[1]


out.name <- paste0(count.method, "_DM_", mode , "_commonDispersion")

mcCores <- 7


## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, subset=Inf, mode = mode, epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-02, mcCores = mcCores, verbose = FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, out.name ,".txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "commonDispersion", mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


dgeDM <- dmTest(dgeDM, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)


write.table(dgeDM$table, paste0(out.dir, out.name ,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir , out.name, "_dgeDM.RData"))


gc()



modeDispList = c("optim", "constrOptim", "grid", "grid")
trendList <- c("none", "none", "none", "commonDispersion")


for(i in length(modeDispList):2){
  # i = 1
  
  extraName <- ""
  
  modeDisp <- modeDispList[i]
  trend <- trendList[i]
  
  out.name <- paste0(count.method, "_DM_", mode , "_", modeDisp, extraName)
  
  if(modeDisp == "grid")
    out.name <- paste0(count.method, "_DM_", mode , "_", modeDisp, "-", trend, extraName)
  
  
  dgeDM <- dmEstimateTagwiseDisp(dgeDM, group = NULL, adjust = TRUE, mode = mode, epsilon = 1e-05, maxIte = 1000, modeDisp = modeDisp, interval = c(0, 1e+10), tol = sqrt(.Machine$double.eps) , initDisp = 10, initWeirMoM = TRUE, gridLength = 21, gridRange = c(-10, 10), trend = trend, priorDf = 4, span = 0.3, mcCores = mcCores, verbose = FALSE)
  
  write.table(dgeDM$tagwiseDispersion, paste0(out.dir, out.name, "_tagwiseDispersion.txt"), quote=F, sep="\t", row.names=T, col.names=F)
  
  pdf(paste0(out.dir, out.name ,"_DispVsMean.pdf"), 7, 5)
  
  df <- data.frame(meanExpr = log10(dgeDM$meanExpr+1), tagwiseDispersion = log10(dgeDM$tagwiseDispersion))
  
  rownames(dgeDM$table) <- dgeDM$table$GeneID
  myPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
  ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = dgeDM$table[names(dgeDM$meanExpr), "df"] )) +
    theme_bw() +
    xlab("Log10 meanExpr") +
    ylab("Log10 tagwiseDispersion") +
    geom_point(size = 2, alpha = 0.5) +
    geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
    scale_colour_gradient(limits=c(1, 10), low = "blue", high="red", name = "df", na.value = "red")
  print(ggp2)
  
  dev.off()
  
  
  dgeDM <- dmFit(dgeDM, group=NULL, dispersion= "tagwiseDispersion", mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
  
  
  dgeDM <- dmTest(dgeDM, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
  
  write.table(dgeDM$table, paste0(out.dir, out.name ,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
  save(dgeDM, file=paste0(out.dir, out.name ,"_dgeDM.RData"))
  
  gc()
  
  
}




































