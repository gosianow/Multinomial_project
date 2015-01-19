# BioC 3.0

# Created 15 Jan 2014
# Modyfied 15 Jan 2014

# Compare the results of sQTLseekeR 2.0 & DM_0.1.2
# Make the plots of transcript usage for interesting SNPs

############################################################################################


setwd("/home/Shared/data/seq/GEUVADIS/")

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(sQTLseekeR)

library(ggplot2)
library(reshape2)

##############################################################
# Results from SQTLseekeR 2.0
##############################################################


out.dir <- "sQTLseekeR20_analysis/Plots/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


################### Load seekeR & DM results

res.seeker <- read.table("sQTLseekeR20_analysis/Results/CEU_results_FDR05.txt", header = TRUE, as.is = TRUE)
res.seeker.all <- read.table("sQTLseekeR20_analysis/Results/CEU_results_all.txt", header = TRUE, as.is = TRUE)

tre.df.rel <- read.table(paste0("sQTLseekeR20_analysis/Data/trExpRPKM_CEU_clean.tsv"), header = TRUE, as.is=TRUE)

## load genotypes
load(paste0("DM_0.1.1_sQTL_analysis/Data/", "genotypes_basedOnSQTLseekeRresults.RData"))

## load tre.df
load(paste0("DM_0.1.1_sQTL_analysis/Data/", "tre.df.RData"))


################### Import Phase1.Geuvadis_dbSnp137_idconvert_snpOnly.txt

snpIDconvert <- read.table(paste0("GEUVADIS_genotypes/Phase1.Geuvadis_dbSnp137_idconvert_snpOnly.txt"), header = FALSE, as.is=TRUE)

save(snpIDconvert, file = paste0("sQTLseekeR20_analysis/Data/snpIDconvert.RData"))

snpIDconvert.filtered <- snpIDconvert[snpIDconvert[, 2] %in% res.seeker.all$snpId, ]
colnames(snpIDconvert.filtered) <- c("snpId_GWAS", "snpId")


save(snpIDconvert.filtered, file = paste0("sQTLseekeR20_analysis/Data/snpIDconvert.filtered.RData"))


res.seeker <- merge(x = res.seeker, y = snpIDconvert.filtered, by="snpId", all.x = TRUE, sort = FALSE)



################## Plot transcirpt expression for diff. genotypes

plot.snps <- res.seeker[1:10, c("geneId", "snpId","snpId_GWAS")]

# plot.snps <- res.seeker[res.seeker$snpId_GWAS == "rs2549794", c("geneId", "snpId", "snpId_GWAS")]
# plot.snps <- res.seeker[res.seeker$snpId_GWAS == "rs17318596", c("geneId", "snpId", "snpId_GWAS")]



for(i in 1:nrow(plot.snps)){
  # i = 1
  gene <- plot.snps[i, 1]
  snp <- plot.snps[i, 2]
  
  expr <- tre.df.rel[tre.df.rel$geneId == gene, -c(1,2)]
  
  rownames(expr) <- tre.df.rel[tre.df.rel$geneId == gene, "trId"]
  
  geno <- genotypes.org[genotypes.org$snpId == snp & genotypes.org$geneId == gene , -c(1:5)]
  
  samps.keep <- !is.na(geno) & !is.na(expr[1,])
  
  expr <- expr[,samps.keep]
  geno <- geno[samps.keep]
  names(geno) <- colnames(expr)
  
  
  tom <- cbind(expr, Transcript=rownames(expr))  
  m <- melt(tom, id.vars = "Transcript" )
  m$genotype <- NA
  
  geno.val <- sort(unique(geno))  
  var.counts <- rep(0, length(geno.val))
  names(var.counts) <- paste0("variant ", geno.val)
  
 for(j in 1:length(geno.val)){
   
   m$genotype[m$variable %in% names(geno[geno == geno.val[j]] )] <- paste0("variant ", geno.val[j])
   var.counts[j] <- length(geno[geno == geno.val[j]])
  
 }
  
  ggplot(data = m, aes(x = genotype, y = value)) + ggtitle(paste0(snp," / ", plot.snps[i, "snpId_GWAS"])) +
    geom_boxplot(aes(fill = Transcript), width = 1)  + scale_x_discrete(labels = paste0(names(var.counts), " (",var.counts, ")" ), name="") + ylab("Splicing ratio")
  
ggsave(paste0(out.dir, "Expression_", gene, "_", snp, ".pdf"), width = 10, height = 7, units = "in")
  

}















































