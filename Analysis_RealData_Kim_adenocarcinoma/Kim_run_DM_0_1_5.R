######################################################
# BioC 3.1

# Created 19 June 2015 
# Analysis of Kim_adenocarcinoma data with DM_0.1.5 using htseq counts
# Updated 22 June 2015
# -  using kallisto counts

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

Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)

DM_out <- "DM_0_1_5/"


##########################################################################
# load metadata
##########################################################################


metadata <- read.table("3_metadata/Malgorzata Nowicka2014-11-04GSE37764.csv", stringsAsFactors=F, sep=",", header=T) 
metadata <- metadata[metadata$X == "RNA-seq",]

metadata$sampleName <- metadata$ids
metadata$condition <- metadata$Tissue.Type

metadata <- metadata[order(metadata$condition), ]

metadata


##############################################################################################################
# htseq counts
##############################################################################################################


### load htseq counts
# htseqList <- lapply(metadata$sampleName, function(i){
#   # i = 1
#   htseq <- read.table(paste0("2_counts/htseq/", i, ".counts"), header = FALSE, as.is = TRUE)
#   colnames(htseq) <- c("group_id", paste0(i))  
#   return(htseq)
# })
# 
# htseq <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), htseqList)
# 
# write.table(htseq, paste0("2_counts/htseq/htseq_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


htseq <- read.table(paste0("2_counts/htseq/htseq_counts.txt"), header = TRUE, as.is = TRUE)

htseq <- htseq[!grepl(pattern = "_", htseq$group_id), ]

expr <- htseq[,-1]
colnames(expr) <- metadata$sampleName
gene_id <- strsplit2(htseq[,1], ":")[,1]
ete_id <- htseq[,1]

dgeOrg <- DGEList(counts=expr, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id, stringsAsFactors = FALSE))

dgeOrg


######################################
# htseq counts: filtering
######################################


############### filtering function

dmFilteringExons <- function(dgeOrg, min_samps_gene_expr = 3, min_gene_expr = 1, min_samps_transcript_prop =3, min_transcript_prop = 0.01, max_transcripts = NA, out.dir){
  
  expr_cpm <- cpm(dgeOrg)
  rownames(expr_cpm) <- dgeOrg$genes$ete_id
  expr_cpm_spl <- split.data.frame(expr_cpm, dgeOrg$genes$gene_id) 
  
  expr <- dgeOrg$counts
  rownames(expr) <- dgeOrg$genes$ete_id
  expr_spl <- split.data.frame(expr, dgeOrg$genes$gene_id) 
  
  
  counts <- lapply(names(expr_spl), function(g){
    # g = "FBgn0002528"
    # print(g)
    expr_cpm_gene <- expr_cpm_spl[[g]]
    expr_gene <- expr_spl[[g]]
    
    ### no genes with one transcript
    if(dim(expr_gene)[1] == 1)
      return(NULL)
    
    ### genes with min expression in all samples
    # if(!all(colSums(expr_cpm_gene) > min_gene_expr))
    if(! sum(colSums(expr_cpm_gene) > min_gene_expr) >= min_samps_gene_expr )
      return(NULL)
    
    samps2keep <- colSums(expr_cpm_gene) != 0
    
    prop <- prop.table(expr_gene[, samps2keep], 2) 
    trans2keep <- rowSums(prop > min_transcript_prop) >= min_samps_transcript_prop
    
    ### no genes with one transcript
    if(sum(trans2keep) <= 1)
      return(NULL)
    
    
    if(!is.na(max_transcripts)){
      if(sum(trans2keep) > max_transcripts){
        
        tr_order <- order(apply(aggregate(t(prop[trans2keep, ]), by = list(Condition = dgeOrg$samples$group[samps2keep]), median)[, -1], 2, max), decreasing = TRUE)
        
        trans2keep <- trans2keep[trans2keep]
        
        trans2keep <- names(trans2keep[tr_order[1:max_transcripts]])
        
      }
    }
    
    expr <- expr_gene[trans2keep, ] 
    
    return(expr)
    
  })
  
  names(counts) <- names(expr_cpm_spl)
  counts2keep <- !sapply(counts, is.null)
  counts <- counts[counts2keep]
  
  dge <- DGEList()
  dge$counts <- counts
  dge$samples <- data.frame(group = dgeOrg$samples$group)
  
  save(dge, file = paste0(out.dir, "/dge.RData"))
  
  
  tt <- sapply(dge$counts, nrow)
  
  pdf(paste0(out.dir, "/Hist_numberOfExons.pdf"))
  hist(tt, breaks = seq(0, max(tt), by = 1), col = "darkseagreen2", main = paste0(length(tt), " genes \n ", sum(tt) , " exons/bins "), xlab = "Number of exons/bins per gene")
  dev.off()
  
  
}



############### run different filterings


model <- "diff_out_Full_tumorVSnormal"
count.method <- "htseq"


min_samps_transcript_prop <- c(6, 6)
min_transcript_prop <- c(0.01, 0.01) # in cpm
min_samps_gene_expr <- c(12, 9)
min_gene_expr <- c(1, 1) # in cpm
max_transcripts <- c(NA, NA)


for(i in 2){
  # i = 1
  print(i)
  
  filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop[i], "prop",min_transcript_prop[i], "min", min_samps_gene_expr[i], "cpm", min_gene_expr[i], "max", max_transcripts[i] )) 
  
  out.dir <- paste0(DM_out, model, "/", count.method,"/",filter.method,"/")
  dir.create(out.dir, showWarnings=F, recursive=T)
  
  dmFilteringExons(dgeOrg, min_samps_gene_expr[i], min_gene_expr[i], min_samps_transcript_prop[i], min_transcript_prop[i], max_transcripts[i], out.dir)
  
  
}





######################################
# htseq counts: filtering from DEXSeq
######################################

model <- "diff_out_Full_tumorVSnormal"
count.method <- "htseq"
filter.method <- "Filtering_DEXSeq"


out.dir <- paste0(DM_out, model, "/", count.method, "/", filter.method, "/")
dir.create(out.dir, showWarnings=F, recursive=T)


### DEXSeq exon results
rt <- read.table(paste0("4_results/DEXSeq_1.10.8/", model, "/DEXSeq_", count.method, "_exon_results.txt"), header = T, as.is = TRUE, sep = "\t", fill = TRUE )
head(rt)


rt <- rt[, c("groupID", "featureID","dispersion", "pvalue", "padj")]
rt <- rt[complete.cases(rt), ]

keep <- paste0(rt$groupID, ":" , substring(rt$featureID, 2))
head(keep)


dge_tmp <- dgeOrg[dgeOrg$genes$ete_id %in% keep,]
dge_tmp$genes$gene_id <- as.character(dge_tmp$genes$gene_id)
rownames(dge_tmp$counts) <- dge_tmp$genes$ete_id
colnames(dge_tmp$counts) <- dge_tmp$samples$group

dge <- DGEList()
dge$counts <- split.data.frame(dge_tmp$counts, dge_tmp$genes$gene_id)
dge$samples <- data.frame(group = dgeOrg$samples$group)
dge

save(dge, file = paste0(out.dir, "/dge.RData"))



tt <- sapply(dge$counts, nrow)

pdf(paste0(out.dir, "/Hist_numberOfExons.pdf"))
hist(tt, breaks = seq(0, max(tt), by = 1), col = "darkseagreen2", main = paste0(length(tt), " genes \n ", sum(tt) , " exons/bins "), xlab = "Number of exons/bins per gene")
dev.off()


##############################################################################################################
# kallisto counts
##############################################################################################################


### load kallisto counts
kallisto <- read.table(paste0("2_counts/kallisto/kallisto_counts.txt"), header = TRUE, as.is = TRUE)

kallisto <- kallisto[!grepl(pattern = "_", kallisto$group_id), ]

expr <- kallisto[,-1]
colnames(expr) <- metadata$sampleName
gene_id <- strsplit2(kallisto[,1], ":")[,1]
ete_id <- kallisto[,1]

dgeOrg <- DGEList(counts=expr, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id, stringsAsFactors = FALSE))

dgeOrg


######################################
# kallisto counts: filtering
######################################


############### filtering function

dmFilteringTranscripts <- function(dgeOrg, min_samps_gene_expr = 3, min_gene_expr = 1, min_samps_transcript_prop =3, min_transcript_prop = 0.01, max_transcripts = NA, out.dir){
  
  expr_cpm <- cpm(dgeOrg)
  rownames(expr_cpm) <- dgeOrg$genes$ete_id
  expr_cpm_spl <- split.data.frame(expr_cpm, dgeOrg$genes$gene_id) 
  
  expr <- dgeOrg$counts
  rownames(expr) <- dgeOrg$genes$ete_id
  expr_spl <- split.data.frame(expr, dgeOrg$genes$gene_id) 
  
  
  counts <- lapply(names(expr_spl), function(g){
    # g = "FBgn0002528"
    # print(g)
    expr_cpm_gene <- expr_cpm_spl[[g]]
    expr_gene <- expr_spl[[g]]
    
    ### no genes with one transcript
    if(dim(expr_gene)[1] == 1)
      return(NULL)
    
    ### genes with min expression in all samples
    # if(!all(colSums(expr_cpm_gene) > min_gene_expr))
    if(! sum(colSums(expr_cpm_gene) > min_gene_expr) >= min_samps_gene_expr )
      return(NULL)
    
    samps2keep <- colSums(expr_cpm_gene) != 0
    
    prop <- prop.table(expr_gene[, samps2keep], 2) 
    trans2keep <- rowSums(prop > min_transcript_prop) >= min_samps_transcript_prop
    
    ### no genes with one transcript
    if(sum(trans2keep) <= 1)
      return(NULL)
    
    
    if(!is.na(max_transcripts)){
      if(sum(trans2keep) > max_transcripts){
        
        tr_order <- order(apply(aggregate(t(prop[trans2keep, ]), by = list(Condition = dgeOrg$samples$group[samps2keep]), median)[, -1], 2, max), decreasing = TRUE)
        
        trans2keep <- trans2keep[trans2keep]
        
        trans2keep <- names(trans2keep[tr_order[1:max_transcripts]])
        
      }
    }
    
    expr <- expr_gene[trans2keep, ] 
    
    return(expr)
    
  })
  
  names(counts) <- names(expr_cpm_spl)
  counts2keep <- !sapply(counts, is.null)
  counts <- counts[counts2keep]
  
  dge <- DGEList()
  dge$counts <- counts
  dge$samples <- data.frame(group = dgeOrg$samples$group)
  
  save(dge, file = paste0(out.dir, "/dge.RData"))
  
  
  tt <- sapply(dge$counts, nrow)
  
  pdf(paste0(out.dir, "/Hist_numberOfTranscripts.pdf"))
  hist(tt, breaks = seq(0, max(tt), by = 1), col = "paleturquoise2", main = paste0(length(tt), " genes \n ", sum(tt) , " transcripts "), xlab = "Number of transcripts per gene")
  dev.off()
  
  
}



############### run different filterings


model <- "diff_out_Full_tumorVSnormal"
count.method <- "kallisto"


min_samps_transcript_prop <- c(6, 6)
min_transcript_prop <- c(0.01, 0.01) # in cpm
min_samps_gene_expr <- c(12, 9)
min_gene_expr <- c(1, 1) # in cpm
max_transcripts <- c(NA, NA)


for(i in 1){
  # i = 1
  print(i)
  
  filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop[i], "prop",min_transcript_prop[i], "min", min_samps_gene_expr[i], "cpm", min_gene_expr[i], "max", max_transcripts[i] )) 
  
  out.dir <- paste0(DM_out, model, "/", count.method,"/",filter.method,"/")
  dir.create(out.dir, showWarnings=F, recursive=T)
  
  dmFilteringTranscripts(dgeOrg, min_samps_gene_expr[i], min_gene_expr[i], min_samps_transcript_prop[i], min_transcript_prop[i], max_transcripts[i], out.dir)
  
  
}





######################################
# kallisto counts: filtering from DEXSeq
######################################

model <- "diff_out_Full_tumorVSnormal"
count.method <- "kallisto"
filter.method <- "Filtering_DEXSeq"


out.dir <- paste0(DM_out, model, "/", count.method, "/", filter.method, "/")
dir.create(out.dir, showWarnings=F, recursive=T)


### DEXSeq exon results
rt <- read.table(paste0("4_results/DEXSeq_1.10.8/", model, "/DEXSeq_", count.method, "_exon_results.txt"), header = T, as.is = TRUE, sep = "\t", fill = TRUE )
head(rt)


rt <- rt[, c("groupID", "featureID","dispersion", "pvalue", "padj")]
rt <- rt[complete.cases(rt), ]

keep <- paste0(rt$groupID, ":" , substring(rt$featureID, 2))
head(keep)


dge_tmp <- dgeOrg[dgeOrg$genes$ete_id %in% keep,]
dge_tmp$genes$gene_id <- as.character(dge_tmp$genes$gene_id)
rownames(dge_tmp$counts) <- dge_tmp$genes$ete_id
colnames(dge_tmp$counts) <- dge_tmp$samples$group

dge <- DGEList()
dge$counts <- split.data.frame(dge_tmp$counts, dge_tmp$genes$gene_id)
dge$samples <- data.frame(group = dgeOrg$samples$group)
dge

save(dge, file = paste0(out.dir, "/dge.RData"))



tt <- sapply(dge$counts, nrow)

pdf(paste0(out.dir, "/Hist_numberOfTranscripts.pdf"))
hist(tt, breaks = seq(0, max(tt), by = 1), col = "darkseagreen2", main = paste0(length(tt), " genes \n ", sum(tt) , " transcripts "), xlab = "Number of transcripts per gene")
dev.off()



##############################################

### run DM with different modes 

##############################################


model <- "diff_out_Full_tumorVSnormal"

count.methodList <- c("htseq", "kallisto")
count.method <- count.methodList[1]

filter.methodList <- c("Filtering_DEXSeq", "Filtering_min6prop0_01min9cpm1maxNA") 
filter.method <- filter.methodList[1]

out.dir <- paste0(DM_out, model, "/", count.method, "/", filter.method, "/")

load(paste0(out.dir, "/dge.RData"))

modePropList = c("constrOptim2G")
modeProp <- modePropList[1]


BPPARAM <- MulticoreParam(workers = 5)



##### run DM pipeline with commonDispersion 


out.name <- paste0(count.method, "_DM_", modeProp , "_", "commonDispersion")


dgeDM <- dmEstimateCommonDisp(dge, adjustDisp = TRUE, intervalDisp = c(0, 1e+5), tolDisp = 1e-00, modeProp = modeProp, tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)

write.table(dgeDM$commonDispersion, paste0(out.dir, out.name ,".txt"), quote=F, sep="\t", row.names=F, col.names=F)


dgeDM <- dmFit(dgeDM, model = c("full", "null")[1], dispersion = "commonDispersion", modeProp = modeProp, tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)


dgeDM <- dmTest(dgeDM, dispersion = "commonDispersion" , modeProp = modeProp, tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)


extraName <- ""
out.name <- paste0(count.method, "_DM_", modeProp , "_", "commonDispersion", extraName)
pdf(paste0(out.dir, out.name ,"_Histogram_pValues.pdf"))
hist(dgeDM$table$pValue, breaks = 100, main = "", col = "hotpink", xlab = "p-values")
dev.off()


write.table(dgeDM$table, paste0(out.dir, out.name ,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)

save(dgeDM, file=paste0(out.dir , out.name, "_dgeDM.RData"))




##### run DM pipeline with tagwiseDispersion

modeDispList = c("optim", "constrOptim", "grid", "grid")
trendDispList <- c("none", "none", "none", "commonDispersion")


for(i in length(modeDispList):3){
  # i = 2
  print(i)
  
  load(paste0(out.dir, "/dge.RData"))
  out.name <- paste0(count.method, "_DM_", modeProp , "_", "commonDispersion")
  dge$commonDispersion <- as.numeric(read.table(paste0(out.dir, out.name ,".txt")))
  
  extraName <- ""
  
  modeDisp <- modeDispList[i]
  trendDisp <- trendDispList[i]
  
  out.name <- paste0(count.method, "_DM_", modeProp , "_", modeDisp, extraName)
  
  if(modeDisp == "grid")
    out.name <- paste0(count.method, "_DM_", modeProp , "_", modeDisp, "-", trendDisp, extraName)
  
  
  dgeDM <- dmEstimateTagwiseDisp(dge, adjustDisp = TRUE, modeDisp = modeDisp, intervalDisp = c(0, 1e+5), tolDisp = 1e-08,  initDisp = 10, initWeirMoMDisp = TRUE, gridLengthDisp = 15, gridRangeDisp = c(-7, 7), trendDisp = trendDisp, priorDfDisp = 5, spanDisp = 0.3, modeProp = modeProp, tolProp = 1e-12, verbose = FALSE, plot = FALSE, BPPARAM = BPPARAM)
  
  
  dgeDM <- dmFit(dgeDM, model = c("full", "null")[1], dispersion = "tagwiseDispersion", modeProp = modeProp, tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)
  
  
  dgeDM <- dmTest(dgeDM, dispersion = "tagwiseDispersion" , modeProp = modeProp, tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)
  
  write.table(dgeDM$table, paste0(out.dir, out.name ,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
  save(dgeDM, file=paste0(out.dir, out.name ,"_dgeDM.RData"))
  
  
  ### plot histogram of p-values
  pdf(paste0(out.dir, out.name ,"_Histogram_pValues.pdf"))
  hist(dgeDM$table$pValue, breaks = 100, main = "", col = "hotpink", xlab = "p-values")
  dev.off()
  
  ### plot dispersion versus mean
  rownames(dgeDM$table) <- dgeDM$table$geneID
  df <- data.frame(meanExpr = log10(dgeDM$meanExpr+1), tagwiseDispersion = log10(dgeDM$tagwiseDispersion), df = dgeDM$table[names(dgeDM$meanExpr), "df"])
  commonDispersion <- log10(dgeDM$commonDispersion)  
  df_quant <- quantile(na.omit(df$df), probs = 0.75)
  
  ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = df )) +
    theme_bw() +
    xlab("Log10 meanExpr") +
    ylab("Log10 tagwiseDispersion") +
    geom_point(size = 1, alpha = 0.3) +
    geom_hline(aes(yintercept = commonDispersion), colour = "deeppink", linetype="dashed", size = 1)+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
    scale_colour_gradient(limits=c(1, df_quant), breaks = seq(2, df_quant, 2), low = "blue", high="red", name = "df", na.value = "red")
  
  pdf(paste0(out.dir, out.name ,"_DispersionVersusMean.pdf"), 7, 5)
  print(ggp2)
  dev.off()
  
  
}



































