##############################################################################

# BioC 3.1
# Created 14 June 2015:

# Run DM_0.1.5 on 
# - htseq couts
# - BitSeq counts

# Updated 16 June 2015
# - filtering counts

##############################################################################


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding")

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)


### Source all R files in DM package / DM_0.1.5
Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])


##############################################################################################################
# load metadata
##############################################################################################################

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


##############################################################################################################
# htseq counts
##############################################################################################################

count.method <- "htseq"

out.dir <- paste0("DM_0_1_5/",count.method,"/")
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



######################################
# htseq counts: filtering
######################################


############### filtering function

dmFilteringExons <- function(dgeOrg, min_samps_gene_expr = 3, min_gene_expr = 1, min_samps_transcript_prop =3, min_transcript_prop = 0.01, max_transcripts = NULL, simu_info_exons, out.dir){
  
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
    
    
    if(!is.null(max_transcripts)){
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
  hist(tt, breaks = max(tt), col = "darkseagreen2", main = paste0(length(tt), " genes \n ", sum(tt) , " exons/bins "), xlab = "Number of exons/bins per gene")
  dev.off()
  
  
  
  #### info about 1000 AS genes
  
  simu_info_exons <- simu_info_exons[, c("gene_id", "transcript_id", "exon_id", "status", "status_trans" ,"status_exon")]
  
  simu_info_exons_spl <- split(simu_info_exons, simu_info_exons$gene_id)
  
  tas <- table(unlist(lapply(names(simu_info_exons_spl), function(g){
    # g = "FBgn0052820"
    as_exons_simu <- simu_info_exons_spl[[g]][simu_info_exons_spl[[g]]$status_exon == 1, "exon_id"]
    
    if(is.null(dge$counts[[g]]))
      return(NA)
    
    as_exons <- sum(rownames(dge$counts[[g]]) %in% as_exons_simu)
    
    return(as_exons)
    
  })), useNA = "always")
  
  
  tas <- tas[c(length(tas), 1:(length(tas) -1))]
  
  names(tas)[1] <- "NG"
  
  colors <- c("darkred", "darkred", rep("grey", 98))
  names(colors) <- c("NG", "0", 1:98)
  colors <- colors[names(tas)]
  
  pdf(paste0(out.dir, "/Hist_filtering.pdf"), width = 7, height = 7)
  xx <- barplot(tas, xlab = "Number of AS exons left within AS gene", ylab = "Number of AS genes", col = colors)
  text(x = xx, y = as.numeric(tas), label = as.numeric(tas), pos = 3)
  dev.off()
  
  
}



############### run different filterings


simu_info_exons <- read.table("Simu_info/true_exons_simulation.txt", header = TRUE, as.is = TRUE)



min_samps_transcript_prop <- 3
min_transcript_prop <- 0.05 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 1 # in cpm
max_transcripts <- NULL

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


dmFilteringExons(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts = NULL, simu_info_exons, out.dir)




min_samps_transcript_prop <- 3
min_transcript_prop <- 0.05 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 0 # in cpm
max_transcripts <- NULL

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)

dmFilteringExons(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts = NULL, simu_info_exons, out.dir)




min_samps_transcript_prop <- 3
min_transcript_prop <- 0.01 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 1 # in cpm
max_transcripts <- NULL

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)

dmFilteringExons(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts = NULL, simu_info_exons, out.dir)




min_samps_transcript_prop <- 3
min_transcript_prop <- 0.01 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 0 # in cpm
max_transcripts <- NULL

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


dmFilteringExons(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts, simu_info_exons, out.dir)




min_samps_transcript_prop <- 3
min_transcript_prop <- 0.01 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 1 # in cpm
max_transcripts <- 10

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


dmFilteringExons(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts, simu_info_exons, out.dir)






min_samps_transcript_prop <- 3
min_transcript_prop <- 0.01 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 0 # in cpm
max_transcripts <- 10

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


dmFilteringExons(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts, simu_info_exons, out.dir)




min_samps_transcript_prop <- 3
min_transcript_prop <- 0.01 # in cpm
min_samps_gene_expr <- 3
min_gene_expr <- 1 # in cpm
max_transcripts <- 10


filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


dmFilteringExons(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts, simu_info_exons, out.dir)


######################################
# htseq counts: filtering from DEXSeq
######################################

filter.method <- "Filtering_DEXSeq"

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


### DEXSeq exon results
rt <- read.table("Results_from_Katarina/dexseq_1.10.8_exon_htseq_results.txt", header = T, as.is = TRUE, sep = "\t", fill = TRUE )
head(rt)

rt <- rt[, c("groupID", "dispersion", "pvalue", "padj")]
rt <- rt[complete.cases(rt), ]

keep <- gsub(pattern = "E", replacement = "", rt$groupID)
length(keep)


dge_tmp <- dgeOrg[dgeOrg$genes$ete_id %in% keep,]
dge_tmp$genes$gene_id <- as.character(dge_tmp$genes$gene_id)
rownames(dge_tmp$counts) <- dge_tmp$genes$ete_id
colnames(dge_tmp$counts) <- dge_tmp$samples$group

dge <- DGEList()
dge$counts <- split.data.frame(dge_tmp$counts, dge_tmp$genes$gene_id)
dge$samples <- data.frame(group = dgeOrg$samples$group)


save(dge, file = paste0(out.dir, "/dge.RData"))



tt <- sapply(dge$counts, nrow)

pdf(paste0(out.dir, "/Hist_numberOfExons.pdf"))
hist(tt, breaks = max(tt), col = "darkseagreen2", main = paste0(length(tt), " genes \n ", sum(tt) , " exons/bins "), xlab = "Number of exons/bins per gene")
dev.off()



#### info about 1000 AS genes

simu_info_exons <- simu_info_exons[, c("gene_id", "transcript_id", "exon_id", "status", "status_trans" ,"status_exon")]

simu_info_exons_spl <- split(simu_info_exons, simu_info_exons$gene_id)

tas <- table(unlist(lapply(names(simu_info_exons_spl), function(g){
  # g = "FBgn0052820"
  as_exons_simu <- simu_info_exons_spl[[g]][simu_info_exons_spl[[g]]$status_exon == 1, "exon_id"]
  
  if(is.null(dge$counts[[g]]))
    return(NA)
  
  as_exons <- sum(rownames(dge$counts[[g]]) %in% as_exons_simu)
  
  return(as_exons)
  
})), useNA = "always")


tas <- tas[c(length(tas), 1:(length(tas) -1))]

names(tas)[1] <- "NG"

colors <- c("darkred", "darkred", rep("grey", 98))
names(colors) <- c("NG", "0", 1:98)
colors <- colors[names(tas)]

pdf(paste0(out.dir, "/Hist_filtering.pdf"), width = 7, height = 7)
xx <- barplot(tas, xlab = "Number of AS exons left within AS gene", ylab = "Number of AS genes", col = colors)
text(x = xx, y = as.numeric(tas), label = as.numeric(tas), pos = 3)
dev.off()



##############################################################################################################
# bitseq counts
##############################################################################################################


count.method <- "bitseq"

out.dir <- paste0("DM_0_1_5/",count.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


bitseq <- read.table("BitSeq_1.10.0/BitSeq_counts.txt", header = TRUE, as.is = TRUE)
head(bitseq)

expr <- bitseq[,3:8]
colnames(expr) <- metadata$SampleName1
gene_id <- bitseq[,"gene"]
ete_id <- paste0(bitseq[, "gene"],":", bitseq[, "transcript"])

dgeOrg <- DGEList(counts=expr, group = metadata$condition, genes=data.frame(gene_id=gene_id, ete_id=ete_id))
dgeOrg$genes$ete_id <- strsplit2(dgeOrg$genes$ete_id, ":")[, 2]
dgeOrg

######################################
# bitseq counts: filtering
######################################



############### filtering function

dmFilteringTranscripts <- function(dgeOrg, min_samps_gene_expr = 3, min_gene_expr = 1, min_samps_transcript_prop =3, min_transcript_prop = 0.01, max_transcripts = NULL, simu_info_exons, out.dir){
  
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
    
    
    if(!is.null(max_transcripts)){
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
  hist(tt, breaks = max(tt), col = "paleturquoise2", main = paste0(length(tt), " genes \n ", sum(tt) , " transcripts "), xlab = "Number of transcripts per gene")
  dev.off()
  
  
  
  #### info about 1000 AS genes
  
  simu_info_exons <- simu_info_exons[, c("gene_id", "transcript_id", "exon_id", "status", "status_trans" ,"status_exon")]
  
  simu_info_exons_spl <- split(simu_info_exons, simu_info_exons$gene_id)
  
  tas <- table(unlist(lapply(names(simu_info_exons_spl), function(g){
    # g = "FBgn0052820"
    as_trans_simu <- unique(simu_info_exons_spl[[g]][simu_info_exons_spl[[g]]$status_trans == 1, "transcript_id"])
    
    if(is.null(dge$counts[[g]]))
      return(NA)
    
    as_trans <- sum(rownames(dge$counts[[g]]) %in% as_trans_simu)
    
    return(as_trans)
    
  })), useNA = "always")
  
  
  tas <- tas[c(length(tas), 1:(length(tas) -1))]
  
  names(tas)[1] <- "NG"
  
  colors <- c("darkred", "darkred", rep("grey", 98))
  names(colors) <- c("NG", "0", 1:98)
  colors <- colors[names(tas)]
  
  pdf(paste0(out.dir, "/Hist_filtering.pdf"), width = 7, height = 7)
  xx <- barplot(tas, xlab = "Number of AS transcripts left within AS gene", ylab = "Number of AS genes", col = colors)
  text(x = xx, y = as.numeric(tas), label = as.numeric(tas), pos = 3)
  dev.off()
  
  
}



############### run different filterings


simu_info_exons <- read.table("Simu_info/true_exons_simulation.txt", header = TRUE, as.is = TRUE)



min_samps_transcript_prop <- 3
min_transcript_prop <- 0.05 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 1 # in cpm
max_transcripts <- NULL

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


dmFilteringTranscripts(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts = NULL, simu_info_exons, out.dir)





min_samps_transcript_prop <- 3
min_transcript_prop <- 0.05 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 0 # in cpm
max_transcripts <- NULL

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)

dmFilteringTranscripts(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts = NULL, simu_info_exons, out.dir)




min_samps_transcript_prop <- 3
min_transcript_prop <- 0.01 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 1 # in cpm
max_transcripts <- NULL

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)

dmFilteringTranscripts(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts = NULL, simu_info_exons, out.dir)




min_samps_transcript_prop <- 3
min_transcript_prop <- 0.01 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 1 # in cpm
max_transcripts <- 10

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


dmFilteringTranscripts(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts, simu_info_exons, out.dir)




min_samps_transcript_prop <- 3
min_transcript_prop <- 0.01 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 0 # in cpm
max_transcripts <- NULL

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


dmFilteringTranscripts(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts, simu_info_exons, out.dir)



min_samps_transcript_prop <- 3
min_transcript_prop <- 0.01 # in cpm
min_samps_gene_expr <- 6
min_gene_expr <- 0 # in cpm
max_transcripts <- 10

filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


dmFilteringTranscripts(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts, simu_info_exons, out.dir)






min_samps_transcript_prop <- 3
min_transcript_prop <- 0.01 # in cpm
min_samps_gene_expr <- 3
min_gene_expr <- 1 # in cpm
max_transcripts <- 10


filter.method <- gsub("\\.", "_" , paste0("Filtering_", "min", min_samps_transcript_prop, "prop",min_transcript_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_transcripts )) 

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


dmFilteringTranscripts(dgeOrg, min_samps_gene_expr, min_gene_expr, min_samps_transcript_prop, min_transcript_prop, max_transcripts, simu_info_exons, out.dir)



######################################
# bitseq counts: filtering from DEXSeq
######################################


filter.method <- "Filtering_DEXSeq"

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)


### DEXSeq exon results
rt <- read.table("DEXSeq_1.10.8/bitseq_1.10.0/dexseq_1.10.8_exon_bitseq_results.txt", header = T, as.is = TRUE, sep = "\t", fill = TRUE )
head(rt)

rt <- rt[, c("groupID", "featureID", "dispersion", "pvalue", "padj")]
rt <- rt[complete.cases(rt), ]

keep <- gsub(pattern = "E", replacement = "", rt$featureID)
length(keep)

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
hist(tt, breaks = max(tt), col = "paleturquoise2", main = paste0(length(tt), " genes \n ", sum(tt) , " transcripts "), xlab = "Number of transcripts per gene")
dev.off()



#### info about 1000 AS genes

simu_info_exons <- simu_info_exons[, c("gene_id", "transcript_id", "exon_id", "status", "status_trans" ,"status_exon")]

simu_info_exons_spl <- split(simu_info_exons, simu_info_exons$gene_id)

tas <- table(unlist(lapply(names(simu_info_exons_spl), function(g){
  # g = "FBgn0052820"
  as_trans_simu <- unique(simu_info_exons_spl[[g]][simu_info_exons_spl[[g]]$status_trans == 1, "transcript_id"])
  
  if(is.null(dge$counts[[g]]))
    return(NA)
  
  as_trans <- sum(rownames(dge$counts[[g]]) %in% as_trans_simu)
  
  return(as_trans)
  
})), useNA = "always")


tas <- tas[c(length(tas), 1:(length(tas) -1))]

names(tas)[1] <- "NG"

colors <- c("darkred", "darkred", rep("grey", 98))
names(colors) <- c("NG", "0", 1:98)
colors <- colors[names(tas)]

pdf(paste0(out.dir, "/Hist_filtering.pdf"), width = 7, height = 7)
xx <- barplot(tas, xlab = "Number of AS transcripts left within AS gene", ylab = "Number of AS genes", col = colors)
text(x = xx, y = as.numeric(tas), label = as.numeric(tas), pos = 3)
dev.off()





##############################################################################################################

### run DM with different modes 

##############################################################################################################


count.methodList <- c("htseq", "bitseq")
count.method <- count.methodList[1]

filter.methodList <- c("Filtering_DEXSeq", "Filtering0", "Filtering001") 
filter.method <- filter.methodList[1]

out.dir <- paste0("DM_0_1_5/",count.method,"/",filter.method,"/")
  
load(paste0(out.dir, "/dge.RData"))

modePropList = c("constrOptim2G")
modeProp <- modePropList[1]


BPPARAM <- MulticoreParam(workers = 10)




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


for(i in length(modeDispList):2){
  # i = 2
  print(i)
  
  load(paste0(out.dir, "/dge_counts_",count.method,".RData"))
  out.name <- paste0(count.method, "_DM_", modeProp , "_", "commonDispersion")
  dge$commonDispersion <- as.numeric(read.table(paste0(out.dir, out.name ,".txt")))
    
  extraName <- ""
  
  modeDisp <- modeDispList[i]
  trendDisp <- trendDispList[i]
  
  out.name <- paste0(count.method, "_DM_", modeProp , "_", modeDisp, extraName)
  
  if(modeDisp == "grid")
    out.name <- paste0(count.method, "_DM_", modeProp , "_", modeDisp, "-", trendDisp, extraName)

  
  dgeDM <- dmEstimateTagwiseDisp(dge, adjustDisp = TRUE, modeDisp = modeDisp, intervalDisp = c(0, 1e+5), tolDisp = 1e-08,  initDisp = 10, initWeirMoMDisp = TRUE, gridLengthDisp = 15, gridRangeDisp = c(-7, 7), trendDisp = trendDisp, priorDfDisp = 10, spanDisp = 0.3, modeProp = modeProp, tolProp = 1e-12, verbose = FALSE, plot = FALSE, BPPARAM = BPPARAM)
  
  
  dgeDM <- dmFit(dgeDM, model = c("full", "null")[1], dispersion = "tagwiseDispersion", modeProp = modeProp, tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)
  
  
  dgeDM <- dmTest(dgeDM, dispersion = "tagwiseDispersion" , modeProp = modeProp, tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)
  
  write.table(dgeDM$table, paste0(out.dir, out.name ,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
  save(dgeDM, file=paste0(out.dir, out.name ,"_dgeDM.RData"))
  
  
  ### plot histogram of p-values
  pdf(paste0(out.dir, out.name ,"_Histogram_pValues.pdf"))
  hist(dgeDM$table$pValue, breaks = 100, main = "", col = "hotpink", xlab = "p-values")
  dev.off()
  
  ### plot dispersion versus mean
  df <- data.frame(meanExpr = log10(dgeDM$meanExpr+1), tagwiseDispersion = log10(dgeDM$tagwiseDispersion))
  rownames(dgeDM$table) <- dgeDM$table$geneID
  
  ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = dgeDM$table[names(dgeDM$meanExpr), "df"] )) +
    theme_bw() +
    xlab("Log10 meanExpr") +
    ylab("Log10 tagwiseDispersion") +
    geom_point(size = 2, alpha = 0.5) +
    geom_hline(aes(yintercept=log10(dgeDM$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
    scale_colour_gradient(limits=c(1, 10), breaks = c(2, 4, 6, 8, 10), low = "blue", high="red", name = "df", na.value = "red")
  
  pdf(paste0(out.dir, out.name ,"_DispersionVersusMean.pdf"), 7, 5)
  print(ggp2)
  dev.off()
  

}
























































