##############################################################################

# BioC 3.1
# Created 10 Sep 2015:

# Run DM_0.2.3 on 
# - htseq couts
# - kallisto counts

# Updated 11 Sep 2015
# - Use htseq counts from INCOMPLETE_KALLISTOEST


##############################################################################

setwd("/home/gosia/multinomial_project/simulations_sim5_drosophila_noDE_noNull/")

setwd("/home/gosia/multinomial_project/simulations_sim5_hsapiens_noDE_noNull/")

setwd("/home/gosia/multinomial_project/simulations_sim5_hsapiens_withDE_noNull/")



library(ggplot2)
library(DM)


########################################################
# load metadata
########################################################

# create metadata file
metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata


##############################################################################################################
# htseq counts
##############################################################################################################

count_method <- "htseq"

out_dir <- paste0("dm_0_2_3/", count_method, "/")
dir.create(out_dir, showWarnings=F, recursive=T)



### table with list of DS exons and DU transcripts
simu_info <- read.table("3_truth/simu_info.txt", header = TRUE, as.is = TRUE)

group_id <- na.omit(simu_info[,"exon_id"])
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]

ds_info <- data.frame(gene_id = gene_id, feature_id = feature_id, stringsAsFactors = FALSE)




### load htseq counts
htseqList <- lapply(1:6, function(i){
  # i = 1
  htseq <- read.table(paste0("2_counts/dexseq_nomerge/dexseq", i, ".txt"), header = FALSE, as.is = TRUE)
  colnames(htseq) <- c("group_id", paste0("sample_", i))  
  return(htseq)
})

htseq <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), htseqList)

# write.table(htseq, paste0("2_counts/dexseq_nomerge/htseq_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# htseq <- read.table("2_counts/dexseq_nomerge/htseq_counts.txt", header = TRUE, as.is = TRUE)

htseq <- htseq[!grepl(pattern = "_", htseq$group_id), ]



counts <- htseq[,-1]
group_id <- htseq[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]

sample_id = metadata$sample_id
group = metadata$group



### DM pipeline


do <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)


min_samps_feature_prop_list <- c(0, 3, 3)
min_feature_prop_list <- c(0, 0.05, 0.01) # in cpm
min_samps_gene_expr_list <- c(0, 6, 6)
min_gene_expr_list <- c(0, 1, 1) # in cpm
max_features_list <- c(Inf, Inf, Inf)


for(i in 1:3){
  
  # i = 1
  print(i)
  
  min_samps_gene_expr = min_samps_gene_expr_list[i]
  min_gene_expr = min_gene_expr_list[i]
  min_samps_feature_prop = min_samps_feature_prop_list[i]
  min_feature_prop = min_feature_prop_list[i]
  max_features = max_features_list[i]
  
  
  filter_method <- gsub("\\.", "_" , paste0("filtering_", "min", min_samps_feature_prop, "prop",min_feature_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_features)) 
  
  out_dir <- paste0("dm_0_2_3/", count_method, "/", filter_method,"/")
  dir.create(out_dir, showWarnings=F, recursive=T)
  
  d <- dmFilter(do, min_samps_gene_expr, min_gene_expr, min_samps_feature_prop, min_feature_prop, max_features)
  
  plotData(d, out_dir)
  
  DM:::dmDS_plotDataInfo(info = DM::counts(d)[, c("gene_id", "feature_id")], ds_info, out_dir = out_dir)
  
  save(d, file = paste0(out_dir, "d.Rdata"))
  
 
}



######################################
# htseq counts: filtering from DEXSeq
######################################

filter_method <- "filtering_dexseq"

out_dir <- paste0("dm_0_2_3/",count_method,"/",filter_method,"/")
dir.create(out_dir, showWarnings=F, recursive=T)


### DEXSeq exon results
library(DEXSeq)
load("4_results/dexseq_htseq_nomerge.Rdata")
# res


rt <- as.data.frame(res[, 1:7])
rt <- rt[, c("groupID", "featureID", "dispersion", "pvalue", "padj")]
keep <- complete.cases(rt)

counts <- htseq[keep,-1]
group_id <- htseq[keep,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]

sample_id = metadata$sample_id
group = metadata$group


d <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)

plotData(d, out_dir)

DM:::dmDS_plotDataInfo(info = DM::counts(d)[, c("gene_id", "feature_id")], ds_info, out_dir = out_dir)

save(d, file = paste0(out_dir, "d.Rdata"))


##############################################################################################################
# htseq counts INCOMPLETE_KALLISTOEST
##############################################################################################################

count_method <- "htseq_kallisto_filter"

out_dir <- paste0("dm_0_2_3/", count_method, "/")
dir.create(out_dir, showWarnings=F, recursive=T)



### table with list of DS exons and DU transcripts
simu_info <- read.table("3_truth/simu_info.txt", header = TRUE, as.is = TRUE)

group_id <- na.omit(simu_info[,"exon_id"])
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]

ds_info <- data.frame(gene_id = gene_id, feature_id = feature_id, stringsAsFactors = FALSE)




### load htseq counts
htseqList <- lapply(1:6, function(i){
  # i = 1
  htseq <- read.table(paste0("2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast5/dexseq", i, ".txt"), header = FALSE, as.is = TRUE)
  colnames(htseq) <- c("group_id", paste0("sample_", i))  
  return(htseq)
})

htseq <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), htseqList)

htseq <- htseq[!grepl(pattern = "_", htseq$group_id), ]

counts <- htseq[,-1]
group_id <- htseq[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]

sample_id = metadata$sample_id
group = metadata$group



### DM pipeline


do <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)


min_samps_feature_prop_list <- c(0, 3, 3)
min_feature_prop_list <- c(0, 0.05, 0.01) # in cpm
min_samps_gene_expr_list <- c(0, 6, 6)
min_gene_expr_list <- c(0, 1, 1) # in cpm
max_features_list <- c(Inf, Inf, Inf)


for(i in 1:3){
  
  # i = 1
  print(i)
  
  min_samps_gene_expr = min_samps_gene_expr_list[i]
  min_gene_expr = min_gene_expr_list[i]
  min_samps_feature_prop = min_samps_feature_prop_list[i]
  min_feature_prop = min_feature_prop_list[i]
  max_features = max_features_list[i]
  
  
  filter_method <- gsub("\\.", "_" , paste0("filtering_", "min", min_samps_feature_prop, "prop",min_feature_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_features)) 
  
  out_dir <- paste0("dm_0_2_3/", count_method, "/", filter_method,"/")
  dir.create(out_dir, showWarnings=F, recursive=T)
  
  d <- dmFilter(do, min_samps_gene_expr, min_gene_expr, min_samps_feature_prop, min_feature_prop, max_features)
  
  plotData(d, out_dir)
  
  DM:::dmDS_plotDataInfo(info = DM::counts(d)[, c("gene_id", "feature_id")], ds_info, out_dir = out_dir)
  
  save(d, file = paste0(out_dir, "d.Rdata"))
  
 
}



######################################
# htseq counts: filtering from DEXSeq
######################################

filter_method <- "filtering_dexseq"

out_dir <- paste0("dm_0_2_3/",count_method,"/",filter_method,"/")
dir.create(out_dir, showWarnings=F, recursive=T)


### DEXSeq exon results
library(DEXSeq)
load("4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast5.Rdata")
# res


rt <- as.data.frame(res[, 1:7])
rt <- rt[, c("groupID", "featureID", "dispersion", "pvalue", "padj")]
keep <- complete.cases(rt)

counts <- htseq[keep,-1]
group_id <- htseq[keep,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]

sample_id = metadata$sample_id
group = metadata$group


d <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)

plotData(d, out_dir)

DM:::dmDS_plotDataInfo(info = DM::counts(d)[, c("gene_id", "feature_id")], ds_info, out_dir = out_dir)

save(d, file = paste0(out_dir, "d.Rdata"))





##############################################################################################################
# kallisto counts
##############################################################################################################


count_method <- "kallisto"

out_dir <- paste0("dm_0_2_3/", count_method, "/")
dir.create(out_dir, showWarnings=F, recursive=T)


simu_info <- read.table("3_truth/simu_info.txt", header = TRUE, as.is = TRUE)

gene_id <- simu_info$gene_id
feature_id <- simu_info$transcript_id

ds_info <- data.frame(gene_id = gene_id, feature_id = feature_id, stringsAsFactors = FALSE)



### load kallisto counts
kallistoList <- lapply(1:6, function(i){
  # i = 1
  kallisto <- read.table(paste0("2_counts/kallisto/kallisto", i, ".txt"), header = FALSE, as.is = TRUE)
  colnames(kallisto) <- c("group_id", paste0("sample_", i))  
  return(kallisto)
})

kallisto <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), kallistoList)

# write.table(kallisto, paste0("2_counts/kallisto/kallisto_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# kallisto <- read.table("2_counts/kallisto/kallisto_counts.txt", header = TRUE, as.is = TRUE)
head(kallisto)





counts <- kallisto[,-1]
group_id <- kallisto[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]

### DM pipeline

do <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)


min_samps_feature_prop_list <- c(0, 3, 3)
min_feature_prop_list <- c(0, 0.05, 0.01) # in cpm
min_samps_gene_expr_list <- c(0, 6, 6)
min_gene_expr_list <- c(0, 1, 1) # in cpm
max_features_list <- c(Inf, Inf, Inf)


for(i in 1:3){
  
  # i = 1
  print(i)
  
  min_samps_gene_expr = min_samps_gene_expr_list[i]
  min_gene_expr = min_gene_expr_list[i]
  min_samps_feature_prop = min_samps_feature_prop_list[i]
  min_feature_prop = min_feature_prop_list[i]
  max_features = max_features_list[i]
  
  
  filter_method <- gsub("\\.", "_" , paste0("filtering_", "min", min_samps_feature_prop, "prop",min_feature_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_features)) 
  
  out_dir <- paste0("dm_0_2_3/", count_method, "/", filter_method,"/")
  dir.create(out_dir, showWarnings=F, recursive=T)
  
  d <- dmFilter(do, min_samps_gene_expr, min_gene_expr, min_samps_feature_prop, min_feature_prop, max_features)
  
  plotData(d, out_dir)
  
  DM:::dmDS_plotDataInfo(info = DM::counts(d)[, c("gene_id", "feature_id")], ds_info, out_dir = out_dir)

  save(d, file = paste0(out_dir, "d.Rdata"))

}


######################################
# kallisto counts: filtering from DEXSeq
######################################


filter_method <- "filtering_dexseq"

out_dir <- paste0("dm_0_2_3/",count_method,"/",filter_method,"/")
dir.create(out_dir, showWarnings=F, recursive=T)


### DEXSeq exon results
library(DEXSeq)
load("4_results/dexseq_kallisto.Rdata")
# res


rt <- as.data.frame(res[, 1:7])
rt <- rt[, c("groupID", "featureID", "dispersion", "pvalue", "padj")]
keep <- complete.cases(rt)


counts <- kallisto[keep,-1]
group_id <- kallisto[keep,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]



d <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)

plotData(d, out_dir)

DM:::dmDS_plotDataInfo(info = DM::counts(d)[, c("gene_id", "feature_id")], ds_info, out_dir = out_dir)

save(d, file = paste0(out_dir, "d.Rdata"))























