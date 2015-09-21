# BioC 3.0

# Run sQTLseekeR 2.0 analysis 

# Created 3 Sep 2015
# Modyfied 13 Sep 2015


setwd("/home/gosia/multinomial_project/simulations_sqtl2_hsapiens_noDE_noNull")


################################################################################
### Prepare data for sQTLSeekeR
################################################################################

data_dir <- "sqtlseeker_2_0/data/"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)




##############################################
### Prepare expression file (use expected counts)
##############################################

# sim_det <- simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, sep = "\t", as.is = TRUE)
# save(simulation_details, file = "3_truth/simulation_details.Rdata")

load("3_truth/simulation_details.Rdata")
sim_det <- simulation_details


counts <- sim_det[, c(2, 1, grep("isoformCount", colnames(sim_det)))]

ns <- length(grep("isoformCount", colnames(sim_det)))

colnames(counts) <- c("trId", "geneId", paste0("s", 1:ns))

head(counts)

write.table(counts, paste0(data_dir,"counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##############################################
### Prepare BED file with gene location
##############################################

# counts <- read.table(paste0(data_dir,"counts.txt"), as.is = TRUE, header = TRUE)

genes <- unique(counts$geneId)
ng <- length(genes)

bed <- data.frame(chr = "chr", start = 1:ng, end = (1:ng)+1, geneId = genes)

write.table(bed, paste0(data_dir,"genes.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



##############################################
### Prepare genotype
##############################################


genotypes1 <- matrix(0, nrow = ng, ncol = ns/2)
genotypes2 <- matrix(1, nrow = ng, ncol = ns/2)


genotypes <- data.frame(chr = "chr", start = 1:ng, end = 1:ng, snpId = paste0("snp", 1:ng), genotypes1, genotypes2)

colnames(genotypes)[-(1:4)] <- paste0("s", 1:ns)


write.table(genotypes, paste0(data_dir,"genotypes.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




################################################################################
### Run sQTLSeekeR
################################################################################

library(sQTLseekeR)
library(ggplot2)
library(GenomicRanges)


# results_dir <- "sqtlseeker_2_0/results/"
# dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)


results_dir <- "sqtlseeker_2_0/results_min_dispersion001/"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)



## 1) Index the genotype file (if not done externally before)

genotypes_file <- paste0(data_dir,"genotypes.tsv")
genotypes_index_file <- index.genotype(genotypes_file)

genotypes_index_file <- paste0(data_dir, "genotypes.tsv.bgz")


## 2) Prepare transcript expression

### paly with min.dispersion
tre_df <- prepare.trans.exp(te.df = counts, min.transcript.exp = 1, min.gene.exp = 1, min.dispersion = 0.01)

write.table(tre_df, paste0(results_dir, "tre_df_seeker.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


### read in the results to have right levels for geneId and trId
# tre_df <- read.table(paste0(results_dir, "tre_df_seeker.tsv"), as.is = TRUE, header = TRUE) 




tt <- table(tre_df$geneId)

df <- data.frame(tt = as.numeric(tt))
main <- "sQTLseekeR"

ggp <- ggplot(df, aes(x = tt)) +
  theme_bw() +
  ggtitle(main) +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  geom_histogram(binwidth = 1, fill = "seagreen4") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
  coord_cartesian(xlim = c(0, max(tt) + 2)) +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(tt), " genes   \n ", sum(tt) , " features   ")), aes(x = x, y = y, label = label), hjust = 1, vjust = 2, size = 6)

pdf(paste0(results_dir, "seeker_hist_features.pdf"))
print(ggp)
dev.off()



## 3) Test gene/SNP associations

bed_file <- paste0(data_dir,"genes.bed")
bed <- read.table(bed_file, as.is = TRUE, header = TRUE)


results <- sqtl.seeker(tre.df = tre_df, genotype.f = genotypes_index_file, gene.loc = bed, genic.window = 0, min.nb.ext.scores = 1000, nb.perm.max = 1e+06, nb.perm.max.svQTL = 10000, svQTL = FALSE, approx = TRUE, verbose = TRUE)

### add adj p-values (as in sqtls.R)
results$qv = qvalue::qvalue(results$pv)$qvalues


write.table(results, paste0(results_dir, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


# ### test if the right genotypes are readed in
# genic.window <- 0

# grgene <- GRanges(Rle(bed[, 1]), IRanges(bed[, 2], bed[, 3]))
# grgene <- resize(grgene, width(grgene) + 2*genic.window, fix="center")

# genotype_gene = sQTLseekeR:::read.bedix(genotypes_index_file, grgene[23])
# genotype_gene



df <- data.frame(pvalues = results$pv)

ggp <- ggplot(df, aes(x = pvalues)) +
  theme_bw() +
  xlab("p-values") +
  ylab("Frequency") +
  geom_histogram(binwidth = 0.01, fill = "deeppink4") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=16, face="bold")) +
  coord_cartesian(xlim = c(-0.02, 1.02))


pdf(paste0(results_dir, "seeker_hist_pvalue.pdf"))
print(ggp)
dev.off()




## 4) Get significant sQTLs
sqtls <- sqtls(res.df = results, FDR = 0.05, md.min = 0.01, out.pdf = NULL, svQTL.removal = TRUE, FDR.svQTL = 0.01)

write.table(sqtls, paste0(results_dir, "sqtls_fdr05.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




################################################################################
### Check the truth 
################################################################################


truth <- read.table("3_truth/truth_human_non_null.txt", as.is = TRUE, header = TRUE, sep = "\t")

genes_ds <- truth[truth$ds_status == 1, "gene"]

table(results$geneId %in% genes_ds)





truth2 <- sim_det[, c("gene_id", "transcript_id", "gene_ds_status", "transcript_ds_status")]

ds_info <- truth2[truth2$transcript_ds_status == 1, c(1, 2)]
colnames(ds_info) <- c("gene_id", "feature_id")


table(duplicated(results$geneId))

info <- tre_df[, c(2, 1)]
colnames(info) <- c("gene_id", "feature_id")


table(info$feature_id %in% ds_info$feature_id)
table(unique(info$gene_id) %in% unique(ds_info$gene_id))


source("/home/gosia/R/multinomial_project/package_devel/DM/R/dmDS_plotDataInfo.R")

dmDS_plotDataInfo(info, ds_info, out_dir = results_dir)













