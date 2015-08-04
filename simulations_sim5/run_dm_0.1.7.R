##############################################################################

# BioC 3.1
# Created 24 July 2015:

# Run DM_0.1.7 on 
# - htseq couts
# - kallisto counts


##############################################################################


setwd("/home/gosia/multinomial_project/simulations_sim5_hsapiens_withDE_noNull/")


library(GenomicRanges)
library(BiocParallel)
library(edgeR)

library(ggplot2)
library(reshape2)

library(rtracklayer)

library(DM)



Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)


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

out_dir <- paste0("dm_0_1_7/", count_method, "/")
dir.create(out_dir, showWarnings=F, recursive=T)


### load htseq counts
htseqList <- lapply(1:6, function(i){
  # i = 1
  htseq <- read.table(paste0("2_counts/dexseq_nomerge/dexseq", i, ".txt"), header = FALSE, as.is = TRUE)
  colnames(htseq) <- c("group_id", paste0("sample_", i))  
  return(htseq)
})

htseq <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), htseqList)

write.table(htseq, paste0("2_counts/dexseq_nomerge/htseq_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


htseq <- read.table("2_counts/dexseq_nomerge/htseq_counts.txt", header = TRUE, as.is = TRUE)
htseq <- htseq[!grepl(pattern = "_", htseq$group_id), ]

counts <- as.matrix(htseq[,-1])
group_id <- htseq[,1]
group_split <- strsplit2(group_id, ":")
gene_id_counts <- group_split[, 1]
feature_id_counts <- group_split[, 2]
sample_id = metadata$sample_id
group = metadata$group


data <- data_org <- dmDSdata(counts = counts, gene_id_counts = gene_id_counts, feature_id_counts = feature_id_counts, sample_id = sample_id, group = group)



######################################
# htseq counts: filtering
######################################

### table with list of DS exons and DU transcripts
simu_info <- read.table("3_truth/simu_info.txt", header = TRUE, as.is = TRUE)

group_id <- na.omit(simu_info[,"exon_id"])
group_split <- strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]

info <- data.frame(gene_id = gene_id, feature_id = feature_id, stringsAsFactors = FALSE)



min_samps_feature_prop_list <- c(3, 3, 3, 3, 3, 3)
min_feature_prop_list <- c(0.05, 0.05, 0.01, 0.01, 0.01, 0.01) # in cpm
min_samps_gene_expr_list <- c(6, 6, 6, 6, 6, 6)
min_gene_expr_list <- c(1, 0, 1, 1, 0, 0) # in cpm
max_features_list <- c(Inf, Inf, Inf, 10, Inf, 10)


for(i in 1:6){
  # i = 1
  print(i)
  
  min_samps_gene_expr = min_samps_gene_expr_list[i]
  min_gene_expr = min_gene_expr_list[i]
  min_samps_feature_prop = min_samps_feature_prop_list[i]
  min_feature_prop = min_feature_prop_list[i]
  max_features = max_features_list[i]
  
  
  filter_method <- gsub("\\.", "_" , paste0("filtering_", "min", min_samps_feature_prop, "prop",min_feature_prop, "min", min_samps_gene_expr, "cpm", min_gene_expr, "max", max_features)) 
  
  out_dir <- paste0("dm_0_1_7/", count_method, "/", filter_method,"/")
  dir.create(out_dir, showWarnings=F, recursive=T)
  
  data <- dmDSfilter(data_org, min_samps_gene_expr = min_samps_gene_expr, min_gene_expr = min_gene_expr, min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, max_features = max_features)
  
  save(data, file = paste0(out_dir, "/data.RData"))

  dmDSplotData(data, out_dir = out_dir, info = info)
  
  
}




######################################
# htseq counts: filtering from DEXSeq
######################################

filter_method <- "filtering_dexseq"

out_dir <- paste0("dm_0_1_7/",count_method,"/",filter_method,"/")
dir.create(out_dir, showWarnings=F, recursive=T)


### DEXSeq exon results
library(DEXSeq)
load("4_results/dexseq_htseq_nomerge.Rdata")
# res


rt <- as.data.frame(res[, 1:7])
rt <- rt[, c("groupID", "featureID", "dispersion", "pvalue", "padj")]
keep <- complete.cases(rt)

counts <- as.matrix(htseq[keep,-1])
group_id <- htseq[keep,1]
group_split <- strsplit2(group_id, ":")
gene_id_counts <- group_split[, 1]
feature_id_counts <- group_split[, 2]
sample_id = metadata$sample_id
group = metadata$group


data <- dmDSdata(counts = counts, gene_id_counts = gene_id_counts, feature_id_counts = feature_id_counts, sample_id = sample_id, group = group)

save(data, file = paste0(out_dir, "/data.RData"))

dmDSplotData(data, out_dir = out_dir, info = info)


##############################################################################################################
# kallisto counts
##############################################################################################################


count_method <- "kallisto"

out_dir <- paste0("dm_0_1_7/",count_method,"/")
dir.create(out_dir, showWarnings=F, recursive=T)


### load kallisto counts
kallistoList <- lapply(1:6, function(i){
  # i = 1
  kallisto <- read.table(paste0("2_counts/kallisto/kallisto", i, ".txt"), header = FALSE, as.is = TRUE)
  colnames(kallisto) <- c("group_id", paste0("sample_", i))  
  return(kallisto)
})

kallisto <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), kallistoList)

write.table(kallisto, paste0("2_counts/kallisto/kallisto_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


kallisto <- read.table("2_counts/kallisto/kallisto_counts.txt", header = TRUE, as.is = TRUE)
head(kallisto)


counts <- kallisto[,-1]
group_id <- kallisto[,1]
group_split <- strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]

data_org <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = metadata$sample_id, group = metadata$group)


######################################
# kallisto counts: filtering
######################################


simu_info <- read.table("3_truth/simu_info.txt", header = TRUE, as.is = TRUE)

gene_id <- simu_info$gene_id
feature_id <- simu_info$transcript_id

info <- data.frame(gene_id = gene_id, feature_id = feature_id, stringsAsFactors = FALSE)



min_samps_feature_prop <- c(3, 3, 3, 3, 3, 3)
min_feature_prop <- c(0.05, 0.05, 0.01, 0.01, 0.01, 0.01) # in cpm
min_samps_gene_expr <- c(6, 6, 6, 6, 6, 6)
min_gene_expr <- c(1, 0, 1, 1, 0, 0) # in cpm
max_features <- c(Inf, Inf, Inf, 10, Inf, 10)


for(i in 1:6){
  # i = 1
  print(i)
  
  filter_method <- gsub("\\.", "_" , paste0("filtering_", "min", min_samps_feature_prop[i], "prop",min_feature_prop[i], "min", min_samps_gene_expr[i], "cpm", min_gene_expr[i], "max", max_features[i] )) 
  
  out_dir <- paste0("dm_0_1_7/", count_method, "/", filter_method,"/")
  dir.create(out_dir, showWarnings=F, recursive=T)
  
  data <- dmDS_filter(data = data_org, min_samps_gene_expr = min_samps_gene_expr[i], min_gene_expr = min_gene_expr[i], min_samps_feature_prop = min_samps_feature_prop[i], min_feature_prop = min_feature_prop[i], max_features = max_features[i])
  
  save(data, file = paste0(out_dir, "/data.RData"))
  
  dmDSplotData(data, out_dir = out_dir, info = info)
  
}






######################################
# kallisto counts: filtering from DEXSeq
######################################


filter_method <- "filtering_dexseq"

out_dir <- paste0("dm_0_1_7/",count_method,"/",filter_method,"/")
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
group_split <- strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]


data <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = metadata$sample_id, group = metadata$group)

save(data, file = paste0(out_dir, "/data.RData"))

dmDS_plotData(data, out_dir = out_dir, info = info)


##############################################################################################################

### run DM with different modes 

##############################################################################################################


########################################################
# run
########################################################


count_method_list <- c("htseq", "kallisto")
count_method <- count_method_list[1]

filter_method_list <- c("filtering_dexseq", "filtering_min3prop0_01min6cpm1maxInf", "filtering_min3prop0_05min6cpm1maxInf") 
filter_method <- filter_method_list[1]

out_dir <- paste0("dm_0_1_7/", count_method, "/", filter_method, "/")

load(paste0(out_dir, "/data.RData"))


# data@samples$group <- factor(rep(c("aa", "bb", "zz"), each = 2))



##### run DM pipeline 


data <- dmDSdispersion(data, mean_expression = TRUE, common_dispersion = TRUE, tagwise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = c("optimize", "optim", "constrOptim", "grid")[4], disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = c("none", "common", "trended")[1], disp_prior_df = 10, disp_span = 0.3, prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers = 10))


data_fit <- dmDSfit(data, dispersion = "tagwise_dispersion", prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers = 10))


table <- dmDStest(data_fit)


### PLOTS
out_name <- paste0(out_dir, "/", count_method, "_DM_common", "_")


dmDSplotTest(table, out_dir = out_name)

dmDSplotFit(data_fit, gene_id = c("FBgn0003721", "FBgn0263325"), plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[1], order = TRUE, plot_full = FALSE, plot_nunll = FALSE, out_dir = out_name)

dmDSplotFit(table, gene_id = c("FBgn0003721", "FBgn0263325"), plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[1], order = TRUE, plot_full = FALSE, plot_nunll = FALSE, out_dir = out_name)

dmDSplotDispersion(table, out_dir = out_name)



write.table(table, paste0(out_name, "results.xls"), quote=F, sep="\t", row.names=F, col.names=T)

save(dispersion, fit, table, file=paste0(out_name, "dmDS.RData"))



##### run DM pipeline with tagwiseDispersion

disp_mode_list = c("optim", "constrOptim", "grid", "grid", "grid")
disp_moderation_list <- c("none", "none", "none", "common", "trended")


for(i in length(disp_mode_list):3){
  # i = 5
  print(i)
  
  disp_mode <- disp_mode_list[i]
  disp_moderation <- disp_moderation_list[i]
  
  out_name <- paste0(out_dir, "/", count_method, "_DM_", disp_mode, "_")
  
  if(disp_mode == "grid")
    out_name <- paste0(out_dir, "/", count_method, "_DM_", disp_mode, "_", disp_moderation, "_")
  
  
  dispersion <- dmDS_estimateTagwiseDispersion(data, disp_adjust = TRUE, disp_mode = disp_mode, disp_interval = c(0, 1e+5), disp_tol = 1e-08, disp_init = common_dispersion, disp_init_weirMoM = FALSE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = disp_moderation, disp_prior_df = 10, disp_span = 0.3, prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers = 10))
  
  
  dmDS_plotDispersion(dispersion, mean_expression, common_dispersion = common_dispersion, out_dir = out_name)
  
  
  fit <- dmDS_fit(data, dispersion, prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers = 10))
  
  
  table <- dmDS_test(fit, verbose = FALSE, BPPARAM = MulticoreParam(workers = 10))
  
  dmDS_plotTable(table, out_dir = out_name)
  
  
  write.table(table, paste0(out_name, "results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
  
  save(dispersion, fit, table, file = paste0(out_name, "dmDS.RData"))
  
  dmDS_plotFit(data, gene_id, fit = fit, table = table, plot_type = c("barplot", "boxplot1", "boxplot2", "lineplot", "ribbonplot")[3], order = TRUE, plot_full = ifelse(is.null(fit), FALSE, TRUE), plot_null = ifelse(is.null(fit), FALSE, TRUE), out_dir = out_name)
  
}



























