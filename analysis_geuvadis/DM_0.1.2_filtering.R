# BioC 3.0

# Created 14 Jan 2014
# Modyfied 15 Jan 2014

# Perform filtering similar to the one in sQTLseekeR analysis of GEUVADIS based on  sQTLseeker results

# Create dgeSQTL object 

# Update 28 May 2015 
# Remove the genotypes with less than 5 replicates

##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

library(GenomicRanges)

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

# library(sQTLseekeR)
# sQTLseekeR.files <- list.files("/home/gosia/R/R_Multinomial_project/Analysis_GEUVADIS/sQTLseekeR-master/R/", full.names = TRUE)
# for(i in sQTLseekeR.files)
# source(i)

out.dir <- "DM_0_1_2_Data/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


##########################################################################################
# Filtering based on sQTLseeker results
##########################################################################################


data.dir <- "sQTLseekeR20_analysis/Data/"

## Input files: transcript expression, gene location and genotype information
trans.exp.f = paste0(data.dir, "trExpRPKM.tsv")
genotype.f = paste0(data.dir, "snps_CEU_full.tsv")


## Getting the IDs of samples in CEU population
groups.f <- paste0(data.dir, "sample-groups.tsv")
groups = read.table(groups.f, header=TRUE, as.is=TRUE)
groupsCEU = subset(groups,group=="CEU")



######################### Filter transcripts and genes

te.df.all = read.table(trans.exp.f, as.is=TRUE, header=TRUE, sep="\t")

te.df = te.df.all[,c("trId", "geneId", groupsCEU$sample)]
colnames(te.df) <- c("trId", "geneId", groupsCEU$sampleShort)

## If I want to have exactly the same filtering I can not use prepare.trans.exp() because there is a sampling step
## I load the relativized results for sQTLseeker and filter transcripts based on that 
## tre.df.rel = prepare.trans.exp(te.df)

tre.df.rel <- read.table(paste0("sQTLseekeR20_analysis/Data/trExpRPKM_CEU_clean.tsv"), header = TRUE, as.is=TRUE)

tre.df <- te.df[te.df$trId %in% tre.df.rel$trId, colnames(tre.df.rel)]
all(tre.df$trId == tre.df.rel$trId)

# x <- c("A", "B", "C")
# y <- c("C", "A", "B")
# y[match(x, y)]


tre.df[is.na(tre.df.rel)] <- NA
all(is.na(tre.df) == is.na(tre.df.rel))

# write.table(tre.df, paste0(out.dir, "tre.df.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(tre.df, file = paste0(out.dir, "tre.df.RData"))

load(paste0(out.dir, "tre.df.RData"))

pdf(paste0(out.dir, "/Hist_numberOfTranscripts.pdf"))
tt <- table(tre.df$geneId)
hist(tt, breaks = max(tt), col = "orangered", main = paste0(length(tt), " genes \n ", sum(tt) , " transcripts "), xlab = "number of transcripts per gene")
dev.off()





######################### SNPs and genotype

### get the SNP - gene match from seeker results 
res.seeker.all <- read.table("sQTLseekeR20_analysis/Results/CEU_results_all.txt", header = TRUE, as.is = TRUE)

snp2gene <- unique(res.seeker.all[, c("geneId", "snpId")])

### load genotypes
# snps_CEU_full <- read.table(genotype.f, header = TRUE, as.is = TRUE)
# save(snps_CEU_full, file = paste0(data.dir, "snps_CEU_full.RData"))

load(paste0(data.dir, "snps_CEU_full.RData"))

snps_CEU_filtered <- snps_CEU_full[snps_CEU_full$snpId %in% unique(res.seeker.all$snpId), ]

##### This does not give the guarantie that there are no groups with < 5 samples; additional filtering is needed.
genotypes <- merge(snp2gene, snps_CEU_filtered, by = "snpId", all.x = TRUE, sort = FALSE)

genotypes[genotypes == -1]  <- NA




#### Check if there are NAs
for(i in 1:ncol(genotypes))
  cat(sum(is.na(genotypes[,i])), " ")

for(i in 1:ncol(genotypes))
  cat(sum(genotypes[,i] == -1), " ")


# write.table(genotypes, paste0(out.dir, "genotypes_basedOnSQTLseekeRresults.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


genotypes$geneId <- as.character(genotypes$geneId)
genotypes$snpId <- as.character(genotypes$snpId)

save(genotypes, file = paste0(out.dir, "genotypes_basedOnSQTLseekeRresults.RData"))

##########################################################################################
# Put NA for the genotypes that have less than 5 replicates
##########################################################################################

load(paste0(out.dir, "/genotypes_basedOnSQTLseekeRresults.RData"))

geno <- as.matrix(genotypes[, -c(1:5)])



for(i in 1:nrow(geno)){
  # i = 85883
  tab <- table(geno[i, ])
  
  if(sum(tab < 5) == 1){
    # print(tab)
    variant <- as.numeric(names(which(tab < 5)))
    geno[i, geno[i, ] == variant ] <- NA
    cat("Less than 5 replicates in ", i, "\n")
  }
  
  
}


genotypes_clean <- data.frame(genotypes[, c(1:5)], geno)

save(genotypes_clean, file = paste0(out.dir, "/genotypes_clean_basedOnSQTLseekeRresults.RData"))


load(paste0(out.dir, "/genotypes_clean_basedOnSQTLseekeRresults.RData"))

tt <- table(genotypes_clean$geneId)
pdf(paste0(out.dir, "/Hist_numberOfSnps.pdf"))
hist(tt, breaks = 100, col = "chartreuse2", main = paste0(length(tt), " genes \n ", sum(tt) , " snps "), xlab = "number of snps per gene")
dev.off()


##########################################################################################
# Create dgeSQTL
##########################################################################################

load(paste0(out.dir, "genotypes_basedOnSQTLseekeRresults.RData"))

### load clean genotypes

load(paste0(out.dir, "/genotypes_clean_basedOnSQTLseekeRresults.RData"))
genotypes <- genotypes_clean

### load transcript expression in RPKM
load(paste0(out.dir, "/tre.df.RData"))











library(edgeR)

dgeSQTL <- DGEList(counts = as.matrix(tre.df[, -c(1, 2)]), genes = tre.df[, c(1, 2)])
rownames(dgeSQTL$counts) <- tre.df[,"trId"]

dgeSQTL$counts <- round(dgeSQTL$counts * 100) ### RPKM -> counts

dgeSQTL$samples <- colnames(tre.df[, -c(1, 2)])
colnames(dgeSQTL$genes) <- c("ete_id", "gene_id")

dgeSQTL$genotypes <- as.matrix(genotypes[, -c(1:5)])
dgeSQTL$SNPs <- genotypes[, c("snpId", "geneId", "chr", "start")]
colnames(dgeSQTL$SNPs)[1:2] <- c("SNP_id", "gene_id")

### check
all(colnames(dgeSQTL$counts) == colnames(dgeSQTL$genotypes))


dgeSQTL$counts <- split(data.frame(dgeSQTL$counts), factor(dgeSQTL$genes$gene_id, levels = unique(dgeSQTL$genes$gene_id)))
dgeSQTL$counts <- lapply(dgeSQTL$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


save(dgeSQTL, file=paste0(out.dir, "/dgeSQTL_clean.RData"))


load(paste0(out.dir, "/dgeSQTL_clean.RData"))


tt <- table(dgeSQTL$SNPs$gene_id)

pdf(paste0(out.dir, "/dgeSQTL_Hist_numberOfSnps.pdf"))
hist(tt, breaks = 100, col = "chartreuse2", main = paste0(length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "Number of SNPs per gene")
dev.off()




tt <- sapply(dgeSQTL$counts[unique(dgeSQTL$SNPs$gene_id)], nrow)

pdf(paste0(out.dir, "/dgeSQTL_Hist_numberOfTranscripts.pdf"))
hist(tt, breaks = max(tt), col = "orangered", main = paste0(length(tt), " genes \n ", sum(tt) , " transcripts "), xlab = "Number of transcripts per gene")
dev.off()








### subset the dgeSQTL object 
keep.genes <- intersect(unique(dgeSQTL$genes$gene_id), unique(dgeSQTL$SNPs$gene_id))[1:20]

dgeSQTL$counts <- dgeSQTL$counts[keep.genes]
dgeSQTL$genes <- dgeSQTL$genes[dgeSQTL$genes$gene_id %in% names(dgeSQTL$counts), ]
dgeSQTL$genotypes <- dgeSQTL$genotypes[dgeSQTL$SNPs$gene_id %in% names(dgeSQTL$counts), ]
dgeSQTL$SNPs <- dgeSQTL$SNPs[dgeSQTL$SNPs$gene_id %in% names(dgeSQTL$counts), ]


dgeSQTL$counts$ENSG00000159733.9


























