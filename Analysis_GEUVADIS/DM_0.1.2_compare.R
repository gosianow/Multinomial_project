# BioC 3.0

# Created 15 Jan 2014
# Modyfied 19 Jan 2014

# Compare the results of sQTLseekeR 2.0 & DM_0.1.2
# Make the plots of transcript usage for interesting SNPs

# Update 22 Apr 2015
# - check why the DM p-values histogram is so strange, why there are so many pvs = 1


############################################################################################


setwd("/home/Shared/data/seq/GEUVADIS/")

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(sQTLseekeR)

library(ggplot2)
library(reshape2)

##############################################################
# Load results
##############################################################



################### Load seekeR & DM results

res.seeker <- read.table("sQTLseekeR20_analysis/Results/CEU_results_FDR05.txt", header = TRUE, as.is = TRUE)
res.seeker.all <- read.table("sQTLseekeR20_analysis/Results/CEU_results_all.txt", header = TRUE, as.is = TRUE)

tre.df.rel <- read.table(paste0("sQTLseekeR20_analysis/Data/trExpRPKM_CEU_clean.tsv"), header = TRUE, as.is=TRUE)

res.dm <- read.table("DM_0_1_2_sQTL_analysis/Results_Data_DM_0_1_2_TagwiseDisp_gridCommonDispersion/CEU_results_FDR05.txt", header = TRUE, as.is = TRUE)
res.dm.all <- read.table("DM_0_1_2_sQTL_analysis/Results_Data_DM_0_1_2_TagwiseDisp_gridCommonDispersion/CEU_results_all.txt", header = TRUE, as.is = TRUE)

## load tre.df
load(paste0("DM_0.1.2_sQTL_analysis/Data/", "tre.df.RData"))
head(tre.df)

## load genotypes
load(paste0("DM_0.1.2_sQTL_analysis/Data/", "genotypes_basedOnSQTLseekeRresults.RData"))
head(genotypes)

tagwDisp <- read.table("DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/tagwiseDispersion.txt", header = TRUE, as.is = TRUE)


################### Import Phase1.Geuvadis_dbSnp137_idconvert_snpOnly.txt
### Add the GWAS snp ID

# snpIDconvert <- read.table(paste0("GEUVADIS_genotypes/Phase1.Geuvadis_dbSnp137_idconvert_snpOnly.txt"), header = FALSE, as.is=TRUE)
# save(snpIDconvert, file = paste0("sQTLseekeR20_analysis/Data/snpIDconvert.RData"))
# snpIDconvert.filtered <- snpIDconvert[snpIDconvert[, 2] %in% res.seeker.all$snpId, ]
# colnames(snpIDconvert.filtered) <- c("snpId_GWAS", "snpId")
# save(snpIDconvert.filtered, file = paste0("sQTLseekeR20_analysis/Data/snpIDconvert.filtered.RData"))

load(paste0("sQTLseekeR20_analysis/Data/snpIDconvert.filtered.RData"))
head(snpIDconvert.filtered)

res.seeker <- merge(x = res.seeker, y = snpIDconvert.filtered, by="snpId", all.x = TRUE, sort = FALSE)
head(res.seeker)

#######################################################
# generate hist of p-values
#######################################################

out.dir <- "Comparison_DM_0.1.2_sQTLseekeR_2.0/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


pdf(paste0(out.dir, "hist_pvalues.pdf"))

hist(res.dm.all[, "PValue"], col = "#FF7F00", breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "DM", cex.main = 3)

hist(res.seeker.all[, "pv"], col = "#1874CD" , breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "sQTLseekeR", cex.main = 3)

dev.off()


#######################################################
# generate venn diagrams 
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/plotVenn.R")

## TP as a separate circle
venne.list <- list()

venne.list[["sQTLseeker"]] <- paste0(res.seeker[,"snpId"], "-", res.seeker[,"geneId"])
# venne.list[["Seeker"]] <- intersect(res.seeker[,"snpId"], res.dm[, "SNP_id"])
venne.list[["DM"]] <- paste0(res.dm[,"SNP_id"], "-", res.dm[, "gene_id"])


plotVenn(venne.list, colors=c(sQTLseeker= "#1874CD", DM= "#FF7F00" ), venn.methods = names(venne.list), margin=0.1, cat.cex=0, cex=1.7, out.dir=out.dir, name2="snp")


## TP as a separate circle
venne.list <- list()

venne.list[["sQTLseeker"]] <- unique(res.seeker[,"geneId"])
venne.list[["DM"]] <- unique(res.dm[,"gene_id"])


plotVenn(venne.list, colors=c(sQTLseeker= "#1874CD", DM= "#FF7F00" ), venn.methods = names(venne.list), margin=0.1, cat.cex=0, cex=1.7, out.dir=out.dir, name2="genes")


#######################################################
# generate sets of snps for plotting - FP for DM 
#######################################################

res.seeker$SNPgene <- paste0(res.seeker[,"snpId"], "-", res.seeker[,"geneId"])
res.dm$SNPgene <- paste0(res.dm[,"SNP_id"], "-", res.dm[, "gene_id"])


fp.dm <- setdiff(res.dm$SNPgene, res.seeker$SNPgene)

res.dm.o <- res.dm[order(res.dm$PValue, decreasing = FALSE), ]

### number of genes with x number of snps

res.dm.t <- table(res.dm.o$gene_id)
res.dm.tt <- table(res.dm.t)
res.dm.tt <- data.frame(Number_of_SNPs = names(res.dm.tt), Number_of_genes = as.numeric(res.dm.tt))

write.table(res.dm.tt, paste0(out.dir, "Number_of_SNPs_per_gene_DM.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

gene <- names(which(res.dm.t == 4738))

head(res.dm.o[res.dm.o$gene_id == gene, ])


res.seeker.t <- table(table(res.seeker$geneId))
res.seeker.t <- data.frame(Number_of_SNPs = names(res.seeker.t), Number_of_genes = as.numeric(res.seeker.t))

write.table(res.seeker.t, paste0(out.dir, "Number_of_SNPs_per_gene_sQTLseekeR.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


res.dm.fp <- res.dm.o[res.dm.o$SNPgene %in% fp.dm, ]

res.dm.t <- table(res.dm.fp$gene_id)
res.dm.tt <- table(res.dm.t)
res.dm.tt <- data.frame(Number_of_SNPs = names(res.dm.tt), Number_of_genes = as.numeric(res.dm.tt))


#### check if this gene is significant / No
# gene <- "ENSG00000164308.12"
# res.dm.fp[res.dm.fp$gene_id == gene, ]

res.dm.fp <- split(res.dm.fp, f = factor(res.dm.fp$gene_id, levels = unique(res.dm.fp$gene_id) ))

### keep genes with large number of snps > 100
gene <- names(which(res.dm.t > 100))

res.dm.fp <- res.dm.fp[gene] 



### keep top snp from each gene
res.dm.fp.top <-  lapply(seq(length(res.dm.fp)), function(g){ res.dm.fp[[g]][1, , drop = FALSE] })

res.dm.fp.top <- do.call(rbind, res.dm.fp.top)
res.dm.fp.top <- res.dm.fp.top[order(res.dm.fp.top$PValue, decreasing = FALSE), ]

res.dm.fp.top$rank <- 1:nrow(res.dm.fp.top)
res.dm.fp.top$nrSNPs <- as.numeric(res.dm.t[res.dm.fp.top$gene_id])


res.dm.fp.top <- merge(res.dm.fp.top, tagwDisp, by.x = "SNPgene", by.y = "SNP_id", all.x = TRUE, sort = FALSE)

#######################################################
# generate sets of snps for plotting - FN for DM 
#######################################################

res.seeker$SNPgene <- paste0(res.seeker[,"snpId"], "-", res.seeker[,"geneId"])
res.dm$SNPgene <- paste0(res.dm[,"SNP_id"], "-", res.dm[, "gene_id"])
res.dm.all$SNPgene <- paste0(res.dm.all[,"SNP_id"], "-", res.dm.all[, "gene_id"])

fn.dm <- setdiff(res.seeker$SNPgene, res.dm$SNPgene)

res.dm.fn <- res.dm.all[res.dm.all$SNPgene %in% fn.dm, ]

res.dm.fn <- merge(res.dm.fn, res.seeker[, c("geneId", "qv")], by.x = "gene_id", by.y = "geneId", all.x = TRUE, sort = FALSE)

res.dm.fn <- res.dm.fn[order(res.dm.fn$qv, decreasing = FALSE), ]


res.dm.fn <- merge(res.dm.fn, tagwDisp, by.x = "SNPgene", by.y = "SNP_id", all.x = TRUE, sort = FALSE)

res.dm.fn <- split(res.dm.fn, f = factor(res.dm.fn$gene_id, levels = unique(res.dm.fn$gene_id)))

### keep top snp from each gene
res.dm.fn.top <-  lapply(seq(length(res.dm.fn)), function(g){ res.dm.fn[[g]][1, , drop = FALSE] })

res.dm.fn.top <- do.call(rbind, res.dm.fn.top)
res.dm.fn.top <- res.dm.fn.top[order(res.dm.fn.top$qv, decreasing = FALSE), ]

res.dm.fn.top$rank <- 1:nrow(res.dm.fn.top)


##############################################################
# Plot raw expression for FP called by DM
##############################################################

# out.dir <- "DM_0.1.2_sQTL_analysis/Plots_TagwiseDisp_gridCommonDispersion/DM_FP/"
# dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
# 
# 
# plot.snps <- res.dm.fp.top[, c("gene_id", "SNP_id", colnames(res.dm.fp.top)[-c(1, 2)])]
# plot.names <- paste0("_rank", res.dm.fp.top$rank, "_nrSNPs", res.dm.fp.top$nrSNPs, "_df", res.dm.fp.top$df)
# plot.main <- paste0(plot.names, "\n", "tagwise dispersion ", res.dm.fp.top$tagwiseDispersion)
# 
# #### gene with over 4'000 significant snps 
# which(plot.snps$gene_id == "ENSG00000189283.5")



out.dir <- "DM_0.1.2_sQTL_analysis/Plots_TagwiseDisp_gridCommonDispersion/DM_FN/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


plot.snps <- res.dm.fn.top[, c("gene_id", "SNP_id", "SNPgene", colnames(res.dm.fn.top)[-c(1, 2, 3)])]
plot.names <- paste0("_rank", res.dm.fn.top$rank)
plot.main <- paste0(plot.names, "\n", "tagwise dispersion ", res.dm.fn.top$tagwiseDispersion, "\n FDR " , res.dm.fn.top$FDR, " // qv ", res.dm.fn.top$qv)



for(i in 1:nrow(plot.snps)){
  # i = 10
  cat(paste0("SNPgene ", i, "\n"))
  
  gene <- plot.snps[i, 1]
  snp <- plot.snps[i, 2]
  
  expr.rel <- tre.df.rel[tre.df.rel$geneId == gene, -c(1,2)]
  expr <- tre.df[tre.df$geneId == gene, -c(1,2)]
  
  rownames(expr.rel) <- tre.df.rel[tre.df.rel$geneId == gene, "trId"]
  rownames(expr) <- tre.df[tre.df$geneId == gene, "trId"]
  
  geno <- genotypes[genotypes$snpId == snp & genotypes$geneId == gene , -c(1:5)]
  
  samps.keep <- !is.na(geno) & !is.na(expr.rel[1,])
  
  expr.rel <- expr.rel[,samps.keep]
  expr <- expr[,samps.keep]
  geno <- geno[samps.keep]
  names(geno) <- colnames(expr.rel)
  
  
  
  ######## Plot Splicing ratio
  
  tom <- cbind(expr.rel, Transcript=rownames(expr.rel))  
  m <- melt(tom, id.vars = "Transcript" )
  m$genotype <- NA
  
  geno.val <- sort(unique(geno))  
  var.counts <- rep(0, length(geno.val))
  names(var.counts) <- paste0("variant ", geno.val)
  
  for(j in 1:length(geno.val)){    
   m$genotype[m$variable %in% names(geno[geno == geno.val[j]] )] <- paste0("variant ", geno.val[j])
    var.counts[j] <- length(geno[geno == geno.val[j]])   
  }
  
  ggplot(data = m, aes(x = genotype, y = value)) + ggtitle(paste0(snp, "-", gene, "\n", plot.main[i])) +
    geom_boxplot(aes(fill = Transcript), width = 1)  + scale_x_discrete(labels = paste0(names(var.counts), " (",var.counts, ")" ), name="") + ylab("Splicing ratio") + theme(panel.grid.major=element_blank(), axis.text=element_text(size=13),axis.title=element_text(size=14,face="bold")) + geom_vline(xintercept=c(1.5,2.5),color="white")
  
  ggsave(paste0(out.dir, "sQTLseekeR_expr",plot.names[i],"_", snp, "-", gene, ".pdf"), width = 15, height = 7, units = "in")
  
  
  
  ######## Plot Expression (RPKM)
  
  # expr <- round(expr * 100)
  
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
  
  ggplot(data = m, aes(x = genotype, y = value)) + ggtitle(paste0(snp, "-", gene, "\n", plot.main[i])) +
    geom_boxplot(aes(fill = Transcript), width = 1)  + scale_x_discrete(labels = paste0(names(var.counts), " (",var.counts, ")" ), name="") + ylab("Expression (RPKM)") + theme(panel.grid.major=element_blank(), axis.text=element_text(size=13),axis.title=element_text(size=14,face="bold")) + geom_vline(xintercept=c(1.5,2.5),color="white")
  
  ggsave(paste0(out.dir, "DM_expr",plot.names[i],"_", snp, "-", gene, ".pdf"), width = 15, height = 7, units = "in")
  
  
}




##############################################################
# Plot transcirpt ratios for sQTLseekeR & raw expression for DM
# Plot validated snp-gene pairs from sQTLseekeR paper
##############################################################

# out.dir <- "sQTLseekeR20_analysis/Plots/"
# dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
# 
#### top significant for seekeR
# plot.snps <- res.seeker[1:10, c("geneId", "snpId")]

#### nps validated from seekeR paper
# plot.snps <- res.seeker[res.seeker$snpId == "snp_19_41937095", c("geneId", "snpId")]
# plot.snps <- res.seeker[res.seeker$snpId == "snp_5_96244549", c("geneId", "snpId")]



for(i in 1:nrow(plot.snps)){
  # i = 1
  gene <- plot.snps[i, 1]
  snp <- plot.snps[i, 2]
  
  expr.rel <- tre.df.rel[tre.df.rel$geneId == gene, -c(1,2)]
  expr <- tre.df[tre.df$geneId == gene, -c(1,2)]
  
  rownames(expr.rel) <- tre.df.rel[tre.df.rel$geneId == gene, "trId"]
  rownames(expr) <- tre.df[tre.df$geneId == gene, "trId"]
  
  geno <- genotypes[genotypes$snpId == snp & genotypes$geneId == gene , -c(1:5)]
  
  samps.keep <- !is.na(geno) & !is.na(expr.rel[1,])
  
  expr.rel <- expr.rel[,samps.keep]
  expr <- expr[,samps.keep]
  geno <- geno[samps.keep]
  names(geno) <- colnames(expr.rel)
  
  
  
  ######## Plot Splicing ratio
  
  tom <- cbind(expr.rel, Transcript=rownames(expr.rel))  
  m <- melt(tom, id.vars = "Transcript" )
  m$genotype <- NA
  
  geno.val <- sort(unique(geno))  
  var.counts <- rep(0, length(geno.val))
  names(var.counts) <- paste0("variant ", geno.val)
  
  for(j in 1:length(geno.val)){    
    m$genotype[m$variable %in% names(geno[geno == geno.val[j]] )] <- paste0("variant ", geno.val[j])
    var.counts[j] <- length(geno[geno == geno.val[j]])   
  }
  
  ggplot(data = m, aes(x = genotype, y = value)) + ggtitle(paste0(snp)) +
    geom_boxplot(aes(fill = Transcript), width = 1)  + scale_x_discrete(labels = paste0(names(var.counts), " (",var.counts, ")" ), name="") + ylab("Splicing ratio") + theme(panel.grid.major=element_blank(), axis.text=element_text(size=13),axis.title=element_text(size=14,face="bold")) + geom_vline(xintercept=c(1.5,2.5),color="white")
  
  ggsave(paste0(out.dir, "sQTLseekeR_expr_",i,"_", gene, "-", snp, ".pdf"), width = 10, height = 7, units = "in")
  
  
  
  ######## Plot Expression (RPKM)
  
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
  
  ggplot(data = m, aes(x = genotype, y = value)) + ggtitle(paste0(snp)) +
    geom_boxplot(aes(fill = Transcript), width = 1)  + scale_x_discrete(labels = paste0(names(var.counts), " (",var.counts, ")" ), name="") + ylab("Expression (RPKM)") + theme(panel.grid.major=element_blank(), axis.text=element_text(size=13),axis.title=element_text(size=14,face="bold")) + geom_vline(xintercept=c(1.5,2.5),color="white")
  
  ggsave(paste0(out.dir, "DM_expr_",i,"_", gene, "-", snp, ".pdf"), width = 10, height = 7, units = "in")
  
  
}
























































