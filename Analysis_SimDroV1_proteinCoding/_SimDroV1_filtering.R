##############################################################################

# BioC 3.0
# Created 28 Jan 2015:

# Compare filtering approaches for DM with DEXSeq filtering


##############################################################################


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding")

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)


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








































