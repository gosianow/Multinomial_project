# BioC 3.0

# Created 12 May 2015
# Analyse data from DM_0.1.2_filtering / Compare with DM_0.1.2_analysis for only one chromosome chr5 

# Updated 18 May 2015
# Analyse data from DM_0.1.3_filtering
# Use different rel.tol in constrOptim for estimating  pi


##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

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




##########################################################################################
# Run the DMsQTL pipeline by chromosome 
##########################################################################################

######### run on DM_0_1_2 data 
out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_2_TagwiseDisp_gridNone/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

load(paste0("DM_0_1_2_Data/dgeSQTL.RData"))



######### run on chr5 DM_0_1_3 data / default tol = sqrt(.Machine$double.eps)
out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

load(paste0("DM_0_1_3_Data/dgeSQTL_chr5.RData"))



######### run on chr5 DM_0_1_3 data / version with tol = 1e-10
out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tol10/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

load(paste0("DM_0_1_3_Data/dgeSQTL_chr5.RData"))




######### run on chr5 DM_0_1_3 data / version with tol = .Machine$double.eps
out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tolEps/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

load(paste0("DM_0_1_3_Data/dgeSQTL_chr5.RData"))



######### run on chr5 DM_0_1_3 data / version with tol = 1e-10
out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tol10_constrOptim2/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

load(paste0("DM_0_1_3_Data/dgeSQTL_chr5.RData"))





dgeSQTL.org <- dgeSQTL

head(dgeSQTL.org$SNPs$SNP_id)
head(dgeSQTL.org$genes$gene_id)

keep.genes <- data.frame(gene_id = intersect(unique(dgeSQTL.org$genes$gene_id), unique(dgeSQTL.org$SNPs$gene_id)), stringsAsFactors = FALSE)
dim(keep.genes)

keep.genes <- merge(keep.genes, unique(dgeSQTL.org$SNPs[, c("gene_id", "chr")]), by = "gene_id", all.x = TRUE, sort = FALSE)
dim(keep.genes)


mcCores <- 20



for(chr in 5){
  # chr = 5
  
  dgeSQTL <- dgeSQTL.org

  ### subset the dgeSQTL object 
  keep.genes.chr <- keep.genes[keep.genes$chr == chr, "gene_id"]

  dgeSQTL$counts <- dgeSQTL$counts[keep.genes.chr]
  dgeSQTL$genes <- dgeSQTL$genes[dgeSQTL$genes$gene_id %in% keep.genes.chr, ]
  dgeSQTL$genotypes <- dgeSQTL$genotypes[dgeSQTL$SNPs$gene_id %in% keep.genes.chr, ]
  dgeSQTL$SNPs <- dgeSQTL$SNPs[dgeSQTL$SNPs$gene_id %in% keep.genes.chr, ]

  save(dgeSQTL, file = paste0(out.dir, "dgeSQTL_chr",chr,".RData"))
#   load(paste0(out.dir, "dgeSQTL_chr",chr,".RData"))
  
  ######### commonDispersion
  cat(paste0("chr", chr," estimate common dispersion \n"))
  
#   dgeSQTL <- dmSQTLEstimateCommonDisp(dgeSQTL, adjust = TRUE, subset = round(nrow(dgeSQTL$genotypes)/170) , mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+2), tol = 1e-02, mcCores=mcCores, verbose=FALSE)
# 
#   write.table(dgeSQTL$commonDispersion, paste0(out.dir, "commonDispersion_chr", chr, ".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  dgeSQTL$commonDispersion <- as.numeric(read.table(paste0(out.dir, "commonDispersion_chr", chr, ".txt" )))
  
  
  ######### tagwiseDispersion
  cat(paste0("chr", chr," estimate tagwise dispersion \n"))
  
  dgeSQTL <- dmSQTLEstimateTagwiseDisp(dgeSQTL, adjust = TRUE, mode = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[2], epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 50), tol = 1e-08,  initDisp = 2, initWeirMoM = TRUE, gridLength = 10, gridRange = c(-6, 6), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 10, span = 0.3, mcCores = mcCores, verbose = FALSE, plot = FALSE)
  
  tagwiseDispersion <- data.frame(SNP_id = names( dgeSQTL$tagwiseDispersion), tagwiseDispersion = dgeSQTL$tagwiseDispersion, stringsAsFactors = FALSE, row.names = NULL)
  
  write.table(tagwiseDispersion, paste0(out.dir, "tagwiseDispersion_chr", chr, ".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
#   tagwiseDispersion <- read.table(paste0(out.dir, "tagwiseDispersion_chr", chr, ".txt" ), header = TRUE)
#   dgeSQTL$tagwiseDispersion <- tagwiseDispersion[, "tagwiseDispersion"]
#   names(dgeSQTL$tagwiseDispersion) <- tagwiseDispersion[, "SNP_id"]
  
  
  ######### DM fitting 
  cat(paste0("chr", chr," fitting DM \n"))
  
  dgeSQTL <- dmSQTLFit(dgeSQTL, model="full", dispersion=c("commonDispersion", "tagwiseDispersion")[2], mode=c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[2], epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores )
  
  
  ######### LR testing
  cat(paste0("chr", chr," testing \n"))
  
  dgeSQTL <- dmSQTLTest(dgeSQTL, dispersion=c("commonDispersion", "tagwiseDispersion")[2], mode=c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[2], epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
  
  save(dgeSQTL, file = paste0(out.dir, "dgeSQTL_chr",chr,".RData"))
  
  write.table(dgeSQTL$table, paste0(out.dir, "results_chr", chr, ".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat(paste0("chr", chr," done \n \n"))
  
  gc()
  
}



### check how many snps are called significant 
table(dgeSQTL$table$FDR < 0.05, useNA = "always")

table(dgeSQTL$table$LR < 0, useNA = "always")

### check how many genes are called significant 
length(unique((dgeSQTL$table$gene_id)))
length(unique((dgeSQTL$table$gene_id[dgeSQTL$table$FDR < 0.05])))



##########################################################################################
#### Choose the results to work with 
##########################################################################################


out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tol10_constrOptim2/"
out.dir.plots <- "DM_0_1_3_sQTL_analysis/Plots_Data_DM_0_1_3_TagwiseDisp_gridNone_tol10_constrOptim2/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)

load(paste0(out.dir,"/dgeSQTL_chr5.RData"))



out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tolEps/"
out.dir.plots <- "DM_0_1_3_sQTL_analysis/Plots_Data_DM_0_1_3_TagwiseDisp_gridNone_tolEps/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)

load(paste0(out.dir,"/dgeSQTL_chr5.RData"))



out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tol10/"
out.dir.plots <- "DM_0_1_3_sQTL_analysis/Plots_Data_DM_0_1_3_TagwiseDisp_gridNone_tol10/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)

load(paste0(out.dir,"/dgeSQTL_chr5.RData"))



out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone/"
out.dir.plots <- "DM_0_1_3_sQTL_analysis/Plots_Data_DM_0_1_3_TagwiseDisp_gridNone/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)

load(paste0(out.dir,"/dgeSQTL_chr5.RData"))



out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_2_TagwiseDisp_gridNone/"
out.dir.plots <- "DM_0_1_3_sQTL_analysis/Plots_Data_DM_0_1_2_TagwiseDisp_gridNone/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)

load(paste0(out.dir,"/dgeSQTL_chr5.RData"))




out.dir <- "DM_0_1_2_sQTL_analysis/Results_Data_DM_0_1_2_TagwiseDisp_gridCommonDispersion/"
out.dir.plots <- "DM_0_1_2_sQTL_analysis/Plots_Data_DM_0_1_2_TagwiseDisp_gridCommonDispersion/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)




##########################################################################################
# Merge DM tables with results and recalculate FDR 
# Merge tagwise dispersion 
# Merge mean gene expression & nr of transcripts
##########################################################################################



##### merge tables with results

res <- lapply(5, function(chr){ 
  # chr=5
  results <- read.table(paste0(out.dir, "results_chr", chr, ".txt" ), header = TRUE, as.is = TRUE)
  
  results <- results[complete.cases(results), ]
  
  pdf(paste0(out.dir, "hist_pvalues_chr",chr,".pdf"))
  
  hist(results[, "PValue"], breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "", col = "#1E90FF")
  hist(results[, "LR"], breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="LR", main = "", col = "#1E90FF")
  
  dev.off()
  
  return(results)
  })
  

res <- do.call(rbind, res)


FDR <- p.adjust(res$PValue, method="BH")

res$FDR <- FDR


res.table <- res

res.table$SNPgene <- paste0(res.table$SNP_id, "-", res.table$gene_id)


write.table(res, paste0(out.dir, "CEU_results_all.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(res[res$FDR < 0.05 & !is.na(res$FDR), ], paste0(out.dir, "CEU_results_FDR05.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




#### check how many genes have LR < 0 

table(res$LR < 0, useNA = "always")




##### merge tagwise dispersion 

res <- lapply(5, function(chr){ results <- read.table(paste0(out.dir, "tagwiseDispersion_chr", chr, ".txt" ), header = TRUE, as.is = TRUE)})

res <- do.call(rbind, res)

write.table(res, paste0(out.dir, "tagwiseDispersion.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



gene <- "ENSG00000169045.13"
snp <- "snp_5_179056159"

res[res$SNP_id == paste(snp, gene, sep="-"), ] # 0.1092517 - this is the values estimated with DM_0.1.2.





##### merge mean gene expression & nr of transcripts

res <- lapply(5, function(chr){ 
  # chr = 5
  
  cat(paste0("chr", chr," \n"))
  
  load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))
  
  meanExpr <- merge(dgeSQTL$SNPs[, c("SNP_id", "gene_id")], data.frame(gene_id = names(dgeSQTL$meanExpr), meanExpr = dgeSQTL$meanExpr, stringsAsFactors = FALSE), by = "gene_id", all.x = TRUE, sort = FALSE)
  
  meanExpr$SNPgene <- paste0(meanExpr$SNP_id, "-", meanExpr$gene_id)
    
  tagwiseDispersion <- data.frame(SNPgene = names(dgeSQTL$tagwiseDispersion), tagwiseDispersion = dgeSQTL$tagwiseDispersion, stringsAsFactors = FALSE)
  meanExpr <- merge(meanExpr, tagwiseDispersion, by = "SNPgene", sort = FALSE)
  
  nrTrans <- t(sapply(dgeSQTL$counts, dim, USE.NAMES = TRUE))  
  nrTrans <- data.frame(gene_id = rownames(nrTrans), nrTrans = nrTrans[, 1], stringsAsFactors = FALSE, row.names = NULL )  
  meanExpr <- merge(meanExpr, nrTrans, by = "gene_id", all.x = TRUE, sort = FALSE)
  
  minSampSize <- apply(dgeSQTL$genotypes, 1, function(gt){ min(table(gt)) })
  minSampSize <- data.frame(SNPgene = paste0(dgeSQTL$SNPs$SNP_id, "-", dgeSQTL$SNPs$gene_id) , minSampSize = minSampSize, stringsAsFactors = FALSE)             
  meanExpr <- merge(meanExpr, minSampSize, by = "SNPgene", sort = FALSE)
  
  return(meanExpr)
  
})

res <- do.call(rbind, res)

res.info <- res

write.table(res, paste0(out.dir, "meanExpr.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




#### merge res.table and res.info
res.all <- merge(res.table, res.info[, c("SNPgene", "meanExpr", "tagwiseDispersion", "nrTrans", "minSampSize")], by = "SNPgene", all.x = TRUE, sort = FALSE)

write.table(res.all, paste0(out.dir, "CEU_results_all_info.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)





### SNPs with negative LR and 2 transcripts
table(res.all$LR < 0, useNA = "always")

res.all[res.all$LR < 0 & res.all$nrTrans == 2, ]

dim(res.all[res.all$LR < 0 & res.all$nrTrans == 3, ])


### check how many snps are called significant 
dim(res)
sum(res$FDR < 0.05, na.rm = TRUE)

### check how many genes are called significant 
length(unique((res$gene_id)))
length(unique((res$gene_id[res$FDR < 0.05])))

##########################################################################################
####### Plot number of SNPs per gene
##########################################################################################




pdf(paste0(out.dir.plots, "/Hist_numberOfSNPsPerGene.pdf"))
tt <- table(res.all$gene_id)

hist(tt, breaks = 100, col = "darkseagreen2", main = paste0("All ",length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "number of SNPs per gene")

tt <- table(res.all[res.all$FDR < 0.05, "gene_id"])

hist(tt, breaks = 100, col = "darkturquoise", main = paste0( "Significant ",length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "number of SNPs per gene")

dev.off()








##########################################################################################
####### Compare LR in res.all with different tol 
##########################################################################################


# tol = 1e-10

out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone/"
res.all.tol8 <- read.table(paste0(out.dir, "CEU_results_all_info.txt" ), header = TRUE, sep = "\t", as.is = TRUE)


out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tol10/"
res.all.tol10 <- read.table(paste0(out.dir, "CEU_results_all_info.txt" ), header = TRUE, sep = "\t", as.is = TRUE)


res <- merge(res.alltol8, res.all.tol10, by = "SNPgene", sort = FALSE)

out.dir.plots.c <- "DM_0_1_3_sQTL_analysis/Plots_Data_DM_0_1_3_Compare/"
dir.create(out.dir.plots.c, showWarnings = FALSE, recursive = TRUE)


png(paste0(out.dir.plots.c, "LRvsLRtol10.png"), 700, 700)
smoothScatter(res[, "LR.x"], res[, "LR.y"], nrpoints = Inf, xlab = "tol = 1e-8", ylab = "tol = 1e-10")
dev.off()

png(paste0(out.dir.plots.c, "LRvsLRtol10_b.png"), 700, 700)
smoothScatter(res[, "LR.x"], res[, "LR.y"], nrpoints = Inf, xlab = "tol = 1e-8", ylab = "tol = 1e-10", xlim = c(-100, 100), ylim = c(-100, 100))
dev.off()




### tol = eps

out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone/"
res.all.tol8 <- read.table(paste0(out.dir, "CEU_results_all_info.txt" ), header = TRUE, sep = "\t", as.is = TRUE)


out.dir <- "DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tolEps/"
res.all.tol10 <- read.table(paste0(out.dir, "CEU_results_all_info.txt" ), header = TRUE, sep = "\t", as.is = TRUE)


res <- merge(res.all.tol8, res.all.tol10, by = "SNPgene", sort = FALSE)

out.dir.plots.c <- "DM_0_1_3_sQTL_analysis/Plots_Data_DM_0_1_3_Compare/"
dir.create(out.dir.plots.c, showWarnings = FALSE, recursive = TRUE)


png(paste0(out.dir.plots.c, "LRvsLRtolEps.png"), 700, 700)
smoothScatter(res[, "LR.x"], res[, "LR.y"], nrpoints = Inf, xlab = "tol = 1e-8", ylab = "tol = eps")
dev.off()

png(paste0(out.dir.plots.c, "LRvsLRtolEps_b.png"), 700, 700)
smoothScatter(res[, "LR.x"], res[, "LR.y"], nrpoints = Inf, xlab = "tol = 1e-8", ylab = "tol = eps", xlim = c(-100, 100), ylim = c(-100, 100))
dev.off()


##########################################################################################
####### Plot LR versus nrTrans and min number of samples with one genotype
##########################################################################################


png(paste0(out.dir.plots, "LRvsNrTrans1.png"), 700, 700)
ggp <- ggplot(res.all, aes(nrTrans, LR)) + geom_point(colour = "blue", alpha = 0.2, position = position_jitter(width = .4))
print(ggp)
dev.off()



png(paste0(out.dir.plots, "LRvsNrTrans2.png"), 800, 700)
ggp <- ggplot(res.all, aes(nrTrans, LR, colour = minSampSize)) + geom_point(alpha = 0.2, position = position_jitter(width = .4)) + scale_color_gradient(low = "red", high = "blue", limits=c(1, 30))
print(ggp)
dev.off()


pdf(paste0(out.dir.plots, "Hist_minSampSize.pdf"))
hist(res.all$minSampSize, breaks = 45, col = "mediumpurple")
hist(res.all$minSampSize[res.all$LR < 0], breaks = 45, col = "mediumpurple")
ggp <- ggplot(res.all, aes(minSampSize, fill = factor(LR < 0))) + geom_density(alpha = 0.5)
print(ggp)
dev.off()


pdf(paste0(out.dir.plots, "Hist_nrTrans.pdf"))
hist(res.all$nrTrans, breaks = 45, col = "mediumpurple")
hist(res.all$nrTrans[res.all$LR < 0], breaks = 45, col = "mediumpurple")
ggp <- ggplot(res.all, aes(nrTrans, fill = factor(LR < 0))) + geom_density(alpha = 0.5)
print(ggp)
dev.off()




##########################################################################################
# Merge proportions piH from dgeSQTL
##########################################################################################


ngroups <- 3
lgroups <- c(0, 1, 2)

res <- mclapply(5, function(chr){ 
  # chr = 22
  
  cat(paste0("chr", chr," \n"))
  
  load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))
  
  piH <- lapply(seq(length(dgeSQTL$fit.null)), function(snp){
    # snp = 1
    if(is.null(dgeSQTL$fit.null[[snp]])) return(NULL)
    
    piHsnp <- matrix(NA, nrow(dgeSQTL$fit[[snp]]$piH), ngroups + 1)
    colnames(piHsnp) <- c(lgroups, "null")
    rownames(piHsnp) <- rownames(dgeSQTL$fit.null[[snp]]$piH)
    
    piHsnp[, colnames(dgeSQTL$fit[[snp]]$piH)] <- dgeSQTL$fit[[snp]]$piH
    piHsnp[, "null"] <- dgeSQTL$fit.null[[snp]]$piH
    
    data.frame(SNPgene = names(dgeSQTL$fit[[snp]]$gamma0), ete_id = rownames(piHsnp), piHsnp, row.names = NULL)    
    
  })
  
  names(piH) <- names(dgeSQTL$tagwiseDispersion)
  
  #   piH <- do.call(rbind, piH) 
  #   write.table(piH, paste0(out.dir, "proportions_chr",chr,".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  save(piH, file = paste0(out.dir, "piH_chr", chr, ".RData" ))
  
  return(NULL)
  
}, mc.cores = 3)







##########################################################################################
# check if snps from sQTLseekeR paper are significant 
##########################################################################################

####### check the significans of snps from sQTLseekeR paper


snp <- "snp_19_41937095"
gene <- "ENSG00000105341.11"

res[res$SNP_id == snp & res$gene_id == gene, ]


snp <- "snp_5_96244549"
gene <- "ENSG00000164308.12"

res[res$SNP_id == snp & res$gene_id == gene, ]

chr = 5

load(paste0(out.dir, "dgeSQTL_chr",chr,".RData"))


snpn <- which(dgeSQTL$SNPs$SNP_id == snp & dgeSQTL$SNPs$gene_id == gene)

dgeSQTL$SNPs[snpn, ]


dgeSQTL$fit[[snpn]]

expr <- dgeSQTL$counts[[gene]]

piH <- dgeSQTL$fit[[snpn]]$piH

rownames(piH) <- rownames(expr)




##########################################################################################
####### Plot p-values distribution
##########################################################################################



pdf(paste0(out.dir, "hist_pvalues.pdf"))

hist(res.table[, "PValue"], col = 1, breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "")

dev.off()




##########################################################################################
# Plot dispersion vs mean expression 
##########################################################################################


res <- mclapply(5, function(chr){ 
  # chr = 22
  
  cat(paste0("chr", chr," \n"))
  
  load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))
  

  m <- match(dgeSQTL$SNPs$gene_id, names(dgeSQTL$meanExpr))
   m.df <- match(names(dgeSQTL$tagwiseDispersion), paste0(dgeSQTL$table$SNP_id, "-",dgeSQTL$table$gene_id) )
  
  
  df <- data.frame(meanExpr = log10(dgeSQTL$meanExpr[m]+1), tagwiseDispersion = log10(dgeSQTL$tagwiseDispersion), df = dgeSQTL$table[m.df, "df"])
  

  table(df$df)
  sum(is.na(df$tagwiseDispersion))
  
  
  pdf(paste0(out.dir, "DispVSmeanExpression_chr", chr, ".pdf" ), 7, 5)
  
#   ggp <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion)) +
#     theme_bw() +
#     geom_point(size = 1) +
#     geom_hline(aes(yintercept=log10(dgeSQTL$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)
#   print(ggp)
  
  #myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  myPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
  ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = df )) +
    theme_bw() +
    geom_point(size = 2) +
  xlab("Log10 mean expression") +
  ylab("Log10 tagwise dispersion") +
    geom_hline(aes(yintercept=log10(dgeSQTL$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)+
    # theme(legend.position="none")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
    #scale_colour_gradientn(colours = myPalette(100), limits=c(1, 15), name = "df")
    scale_colour_gradient(limits=c(1, 30), low = "blue", high="red", name = "df", na.value = "red")
  print(ggp2)
  
  dev.off()
  
gc()

  return(NULL)
  
}, mc.cores = 3)







##########################################################################################
# Plot proportions of some interesting SNPs
##########################################################################################


library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


plotProportions <- function(dgeSQTL, plot.snps, plotPath){
  
  
  pdf(plotPath, width = 12, height = 7)
  
  
  for(i in 1:nrow(plot.snps)){
    # i = 1
    cat(paste0("SNPgene ", i, "\n"))
    
    gene <- as.character(plot.snps[i, 1])
    snp <- as.character(plot.snps[i, 2])
    snpgene <- paste0(snp, "-", gene)
    
    
    index <- which(dgeSQTL$SNPs$SNP_id == snp & dgeSQTL$SNPs$gene_id == gene)
    
    Condition <- factor(dgeSQTL$genotypes[index, ])
    expr <- dgeSQTL$counts[[gene]]
    
    nas <- is.na(Condition) | is.na(expr[1, ])
    
    
    Condition <- Condition[!nas]
    Condition.counts <- table(Condition)
    
    expr <- expr[, !nas]
    
    
    labels <- rownames(expr)
    prop.smp <- prop.table(expr, 2)
    
    #### order transcipts by decreasing proportions 
    order.tr <- labels[order(apply(aggregate(t(prop.smp), by = list(Condition = Condition), median)[, -1], 2, max), decreasing = TRUE)]   
    
    prop.est <- dgeSQTL$fit[[index]]$piH
    rownames(prop.est) <- labels
    
    prop.est.null <- data.frame(ete_id = labels, Proportions = dgeSQTL$fit.null[[index]]$piH)
    prop.est.null$ete_id <- factor(prop.est.null$ete_id, levels = order.tr)
    
    
    prop.smp.m <- melt(prop.smp, varnames = c("ete_id", "Samples"), value.name = "Proportions") 
    prop.smp.m$ete_id <- factor(prop.smp.m$ete_id, levels = order.tr)
    prop.smp.m$Samples <- factor(prop.smp.m$Samples)
    prop.smp.m$Condition <- rep(Condition, each = nrow(prop.smp))
    
    prop.est.m <- melt(prop.est, varnames = c("ete_id", "Condition"), value.name = "Proportions")
    prop.est.m$ete_id <- factor(prop.est.m$ete_id, levels = order.tr)
    prop.est.m$Condition <- factor(prop.est.m$Condition)
    
    
    
    index.table <- which(dgeSQTL$table$SNP_id == snp & dgeSQTL$table$gene_id == gene)
    
    main <- paste0(gene, " - ", snp, "\n meanExpression = ", round(dgeSQTL$meanExpr[gene]), " / TagwiseDispersion = ", round(dgeSQTL$tagwiseDispersion[snpgene], 2), "\n LR = ", round(dgeSQTL$table[index.table, "LR"], 4) , " / PValue = ", sprintf("%.02e", dgeSQTL$table[index.table, "PValue"]))
    
    
    ### box plots with points
    ggb <- ggplot(prop.smp.m, aes(x = ete_id, y = Proportions)) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), legend.position="none", plot.title = element_text(size=10)) +
      ggtitle(main) +     
      geom_jitter(aes(fill = Condition, colour = factor(Condition)), position = position_jitterdodge(dodge.width = 0.75), alpha = 0.5) +
      geom_boxplot(aes(colour = Condition), fill = "white", outlier.size = NA, alpha = 0) + 
      geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, fill = Condition), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 2, shape = 19, colour = "black") +
      geom_point(data = prop.est.null, aes(x = ete_id, y = Proportions), size = 3, shape = 18, colour = "orange") +
      coord_cartesian(ylim = c(-0.1, 1.1)) 
    
    
    print(ggb)
    
    #     ### box plots per condition
    #     ggb <- ggplot(prop.smp.m, aes(x = Condition, y = Proportions)) +
    #       theme_bw() + 
    #       theme(axis.text.x = element_text(angle = 0, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), plot.title = element_text(size=10)) +
    #       ggtitle(main) +     
    #       scale_x_discrete(labels = paste0(names(Condition.counts), " (",Condition.counts, ")" ), name="Variant") +
    #       geom_boxplot(aes(fill = ete_id), width = 1) + 
    #       # geom_point(data = prop.est.m, aes(x = Condition, y = Proportions, fill = ete_id), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 2, shape = 19, colour = "red") +
    #       coord_cartesian(ylim = c(-0.1, 1.1)) +
    #       geom_vline(xintercept=c(1.5,2.5), color="white")
    #     
    #     
    #     print(ggb)
    
    
  }
  
  dev.off()
  
}





#### SNP with most genative LR in DM_0.1.2 analysis
gene <- "ENSG00000169045.13"
snp <- "snp_5_179056159"


plot.snps <- data.frame(gene_id = gene, SNP_id = snp, stringsAsFactors = FALSE)

plotPath <- paste0(out.dir.plots, "Proportions.pdf")


plotProportions(dgeSQTL, plot.snps, plotPath)





### most significant SNPs

plot.snps <- res.all[1:6, c("gene_id", "SNP_id")]


plotPath <- paste0(out.dir.plots, "Proportions_most_significant.pdf")

plotProportions(dgeSQTL, plot.snps, plotPath)




### SNPs with negative LR and 2 transcripts
plot.snps <- res.all[res.all$LR < 0 & res.all$nrTrans == 2, c("gene_id", "SNP_id")]

plotPath <- paste0(out.dir.plots, "Proportions_negativeLR_2trans.pdf")

plotProportions(dgeSQTL, plot.snps, plotPath)



### SNPs with negative LR and 3 transcripts

res.plot <- res.all[res.all$LR < 0 & res.all$nrTrans == 3, ]
res.plot <- res.plot[order(res.plot$LR, decreasing = FALSE), ]
res.plot <- res.plot[!duplicated(res.plot$gene_id), ]


plot.snps <- head(res.plot[, c("gene_id", "SNP_id")])

plotPath <- paste0(out.dir.plots, "Proportions_negativeLR_3trans.pdf")

plotProportions(dgeSQTL, plot.snps, plotPath)




### SNPs with the most negative LR 

res.plot <- res.all[res.all$LR < 0, ]
res.plot <- res.plot[order(res.plot$LR, decreasing = FALSE), ]
res.plot <- res.plot[!duplicated(res.plot$gene_id), ]


plot.snps <- head(res.plot[, c("gene_id", "SNP_id")])

plotPath <- paste0(out.dir.plots, "Proportions_mostNegativeLR.pdf")

plotProportions(dgeSQTL, plot.snps, plotPath)







plot.snps <- head(res.all[res.all$minSampSize == 1, c("gene_id", "SNP_id")])

plotPath <- paste0(out.dir.plots, "Proportions_minSampSize1.pdf")

plotProportions(dgeSQTL, plot.snps, plotPath)




plot.snps <- head(res.all[res.all$minSampSize == 2, c("gene_id", "SNP_id")])

plotPath <- paste0(out.dir.plots, "Proportions_minSampSize2.pdf")

plotProportions(dgeSQTL, plot.snps, plotPath)




























