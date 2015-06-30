# BioC 3.0

# Created 26 May 2015
# Analyse data from DM_0.1.3_filtering (only chr5) with DM_0.1.4

# Updated 28 May 2015
# Run on DM_0_1_2 data but the clean (NA for variants with less than 5 replicates) genotypes


# Updated 11 June 2015
# Run on DM_0_1_5 data

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


Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM_0.1.4/R/", full.names=TRUE)
for(i in Rfiles) source(i)

library(pryr)


##########################################################################################
# Run the DMsQTL pipeline by chromosome 
##########################################################################################


######### run on chr5 DM_0_1_3 data / version with tol = 1e-12
out.dir <- "DM_0_1_4_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tol12_constrOptim2G/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

load(paste0("DM_0_1_3_Data/dgeSQTL_chr5.RData"))




######### run on DM_0_1_2 data / clean version
out.dir <- "DM_0_1_4_sQTL_analysis/Results_Data_DM_0_1_2_clean_TagwiseDisp_gridNone_tol12_constrOptim2G/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

load(paste0("DM_0_1_2_Data/dgeSQTL_clean.RData"))



######### run on DM_0_1_5 data
out.dir <- "DM_0_1_4_sQTL_analysis/Results_Data_DM_0_1_5_TagwiseDisp_gridNone_tol12_constrOptim2G/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

# load(paste0("DM_0_1_5_Data/dgeSQTL.RData"))
# 
# dgeSQTL$samples <- colnames(dgeSQTL$counts[[1]])
# 
# genes <- lapply(names(dgeSQTL$counts), function(g){
#   data.frame(ete_id = rownames(dgeSQTL$counts[[g]]), gene_id = g, stringsAsFactors = FALSE)
# })
# dgeSQTL$genes <- do.call(rbind, genes)
# 
# 
# SNPs <- lapply(names(dgeSQTL$genotypes), function(g){
#   data.frame(SNP_id = rownames(dgeSQTL$genotypes[[g]]), gene_id = g, stringsAsFactors = FALSE)
#   # matrix(c(rownames(dgeSQTL$genotypes[[g]]), rep(g, nrow(dgeSQTL$genotypes[[g]]))), byrow = FALSE, dimnames = c("SNP_id", "gene_id"))
# })
# dgeSQTL$SNPs <- do.call(rbind, SNPs)
# dgeSQTL$SNPs$chr <- strsplit2(dgeSQTL$SNPs$SNP_id, "_")[, 2]
# 
# dgeSQTL$genotypes <- do.call(rbind, dgeSQTL$genotypes)
# rownames(dgeSQTL$genotypes) <- NULL
# 
# save(dgeSQTL, file = paste0("DM_0_1_5_Data/dgeSQTL_DM_0_1_4.RData"))


load(paste0("DM_0_1_5_Data/dgeSQTL_DM_0_1_4.RData"))





#################################
######## run the DM pipeline
#################################

head(dgeSQTL$SNPs$SNP_id)
head(dgeSQTL$genes$gene_id)

keep.genes <- data.frame(gene_id = intersect(unique(dgeSQTL$genes$gene_id), unique(dgeSQTL$SNPs$gene_id)), stringsAsFactors = FALSE)
dim(keep.genes)

keep.genes <- merge(keep.genes, unique(dgeSQTL$SNPs[, c("gene_id", "chr")]), by = "gene_id", all.x = TRUE, sort = FALSE)
dim(keep.genes)


mcCores <- 5


dgeSQTL_org <- dgeSQTL


DM_pipeline_per_chr <- function(dgeSQTL_org, keep.genes, chr, out.dir, mcCores){
  
  keep.genes.chr <- keep.genes[keep.genes$chr == chr, "gene_id"]
  
  dgeSQTL <- DGEList()
  dgeSQTL$counts <- dgeSQTL_org$counts[keep.genes.chr]
  dgeSQTL$genes <- dgeSQTL_org$genes[dgeSQTL_org$genes$gene_id %in% keep.genes.chr, ]
  dgeSQTL$genotypes <- dgeSQTL_org$genotypes[dgeSQTL_org$SNPs$gene_id %in% keep.genes.chr, ]
  dgeSQTL$SNPs <- dgeSQTL_org$SNPs[dgeSQTL_org$SNPs$gene_id %in% keep.genes.chr, ]
  
  # save(dgeSQTL, file = paste0(out.dir, "dgeSQTL_chr",chr,".RData"))
  #   load(paste0(out.dir, "dgeSQTL_chr",chr,".RData"))
  
  
  ######### commonDispersion
  cat(paste0("chr", chr," estimate common dispersion \n"))
  
  dgeSQTL <- dmSQTLEstimateCommonDisp(dgeSQTL, adjust = TRUE, subset = length(dgeSQTL$counts) , mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+2), tol = 1e-02, mcCores=mcCores, verbose=FALSE)
  
  write.table(dgeSQTL$commonDispersion, paste0(out.dir, "commonDispersion_chr", chr, ".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  #   dgeSQTL$commonDispersion <- as.numeric(read.table(paste0(out.dir, "commonDispersion_chr", chr, ".txt" )))
  #   dgeSQTL$commonDispersion <- 4
  
  
  
  ######### tagwiseDispersion
  cat(paste0("chr", chr," estimate tagwise dispersion \n"))
  
  dgeSQTL <- dmSQTLEstimateTagwiseDisp(dgeSQTL, adjust = TRUE, mode = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 50), tol = 1e-08,  initDisp = 2, initWeirMoM = TRUE, gridLength = 10, gridRange = c(-6, 6), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 10, span = 0.3, mcCores = mcCores, verbose = FALSE, plot = FALSE)
  
  tagwiseDispersion <- data.frame(SNP_id = names( dgeSQTL$tagwiseDispersion), tagwiseDispersion = dgeSQTL$tagwiseDispersion, stringsAsFactors = FALSE, row.names = NULL)
  
  write.table(tagwiseDispersion, paste0(out.dir, "tagwiseDispersion_chr", chr, ".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  #   tagwiseDispersion <- read.table(paste0(out.dir, "tagwiseDispersion_chr", chr, ".txt" ), header = TRUE)
  #   dgeSQTL$tagwiseDispersion <- tagwiseDispersion[, "tagwiseDispersion"]
  #   names(dgeSQTL$tagwiseDispersion) <- tagwiseDispersion[, "SNP_id"]
  
  
  ######### DM fitting 
  cat(paste0("chr", chr," fitting DM \n"))
  
  dgeSQTL <- dmSQTLFit(dgeSQTL, model="full", dispersion=c("commonDispersion", "tagwiseDispersion")[2], mode=c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores )
  
  
  ######### LR testing
  cat(paste0("chr", chr," testing \n"))
  
  dgeSQTL <- dmSQTLTest(dgeSQTL, dispersion=c("commonDispersion", "tagwiseDispersion")[2], mode=c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
  
  save(dgeSQTL, file = paste0(out.dir, "dgeSQTL_chr",chr,".RData"))
  
  write.table(dgeSQTL$table, paste0(out.dir, "results_chr", chr, ".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat(paste0("chr", chr," done \n \n"))
  
  return(NULL)
  
}



mem_change_chr <- rep(NA, 22)
mem_used_chr <- rep(NA, 22)


for(chr in 1:22){
  # chr = 22
  
  mem_change_chr[chr] <- mem_change(DM_pipeline_per_chr(dgeSQTL_org, keep.genes, chr, out.dir, mcCores))
  mem_used_chr[chr] <- mem_used()
  
  
}









##########################################################################################
#### Choose the results to work with 
##########################################################################################


out.dir <- "DM_0_1_4_sQTL_analysis/Results_Data_DM_0_1_2_clean_TagwiseDisp_gridNone_tol12_constrOptim2G/"
out.dir.plots <- "DM_0_1_4_sQTL_analysis/Plots_Data_DM_0_1_2_clean_TagwiseDisp_gridNone_tol12_constrOptim2G/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)




out.dir <- "DM_0_1_4_sQTL_analysis/Results_Data_DM_0_1_5_TagwiseDisp_gridNone_tol12_constrOptim2G/"
out.dir.plots <- "DM_0_1_4_sQTL_analysis/Plots_Data_DM_0_1_5_TagwiseDisp_gridNone_tol12_constrOptim2G/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)


##########################################################################################
# Merge DM tables with results and recalculate FDR 
# Merge tagwise dispersion 
# Merge mean gene expression & nr of transcripts
##########################################################################################



##### merge tables with results

res <- lapply(1:22, function(chr){ 
  # chr = 1
  results <- read.table(paste0(out.dir, "results_chr", chr, ".txt" ), header = TRUE, as.is = TRUE)
  
  results <- results[complete.cases(results), ]
  
  results[!complete.cases(results), ]
  
#   pdf(paste0(out.dir, "hist_pvalues_chr",chr,".pdf"))
#   
#   hist(results[, "PValue"], breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "", col = "#1E90FF")
#   hist(results[, "LR"], breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="LR", main = "", col = "#1E90FF")
#   
#   dev.off()
  
  return(results)
  })
  

res <- do.call(rbind, res)


FDR <- p.adjust(res$PValue, method="BH")

res$FDR <- FDR


res.table <- res

res.table$SNPgene <- paste0(res.table$SNP_id, "-", res.table$gene_id)

# write.table(res, paste0(out.dir, "CEU_results_all.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



#### check how many genes have LR < 0 

table(res$LR < 0, useNA = "always")

res[res$LR < 0, ]



##### merge mean gene expression & nr of transcripts

res <- lapply(1:22, function(chr){ 
  # chr = 1
  
  cat(paste0("chr", chr," \n"))
  
  load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))
  
  meanExpr <- merge(dgeSQTL$SNPs[, c("SNP_id", "gene_id")], data.frame(gene_id = names(dgeSQTL$meanExpr), meanExpr = dgeSQTL$meanExpr, stringsAsFactors = FALSE), by = "gene_id", all.x = TRUE, sort = FALSE)
  
  meanExpr$SNPgene <- paste0(meanExpr$SNP_id, "-", meanExpr$gene_id)
    
  tagwiseDispersion <- data.frame(SNPgene = names(dgeSQTL$tagwiseDispersion), tagwiseDispersion = dgeSQTL$tagwiseDispersion, stringsAsFactors = FALSE)
  meanExpr <- merge(meanExpr, tagwiseDispersion, by = "SNPgene", sort = FALSE)
  
#   nrTrans <- t(sapply(dgeSQTL$counts, dim, USE.NAMES = TRUE))  
#   nrTrans <- data.frame(gene_id = rownames(nrTrans), nrTrans = nrTrans[, 1], stringsAsFactors = FALSE, row.names = NULL )  
#   meanExpr <- merge(meanExpr, nrTrans, by = "gene_id", all.x = TRUE, sort = FALSE)
#   
#   minSampSize <- apply(dgeSQTL$genotypes, 1, function(gt){ min(table(gt)) })
#   minSampSize <- data.frame(SNPgene = paste0(dgeSQTL$SNPs$SNP_id, "-", dgeSQTL$SNPs$gene_id) , minSampSize = minSampSize, stringsAsFactors = FALSE)             
#   meanExpr <- merge(meanExpr, minSampSize, by = "SNPgene", sort = FALSE)
  
  return(meanExpr)
  
})

res <- do.call(rbind, res)

res.info <- res


#### merge res.table and res.info
# res.all <- merge(res.table, res.info[, c("SNPgene", "meanExpr", "tagwiseDispersion", "nrTrans", "minSampSize")], by = "SNPgene", all.x = TRUE, sort = FALSE)

res.all <- merge(res.table, res.info[, c("SNPgene", "meanExpr", "tagwiseDispersion")], by = "SNPgene", all.x = TRUE, sort = FALSE)

write.table(res.all, paste0(out.dir, "CEU_results_all_info.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##########################################################################################
# Plot dispersion vs mean expression 
##########################################################################################


res <- mclapply(1:22, function(chr){ 
  # chr = 22
  
  cat(paste0("chr", chr," \n"))
  
  load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))
  
  m <- match(dgeSQTL$SNPs$gene_id, names(dgeSQTL$meanExpr))
  m.df <- match(names(dgeSQTL$tagwiseDispersion), paste0(dgeSQTL$table$SNP_id, "-",dgeSQTL$table$gene_id) )
  
  df <- data.frame(meanExpr = log10(dgeSQTL$meanExpr[m]+1), tagwiseDispersion = log10(dgeSQTL$tagwiseDispersion), df = dgeSQTL$table[m.df, "df"])
  
  pdf(paste0(out.dir, "DispVSmeanExpression_chr", chr, ".pdf" ), 7, 5)
  
  myPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
  ggp2 <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion, colour = df )) +
    theme_bw() +
    geom_point(size = 2) +
    xlab("Log10 mean expression") +
    ylab("Log10 tagwise dispersion") +
    geom_hline(aes(yintercept=log10(dgeSQTL$commonDispersion)), colour = "deeppink", linetype="dashed", size = 1)+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14)) +
    scale_colour_gradient(limits=c(1, 30), low = "blue", high="red", name = "df", na.value = "red")
  print(ggp2)
  
  dev.off()
  
  return(NULL)
  
}, mc.cores = 5)




##########################################################################################
####### Plot p-values distribution
##########################################################################################


res.table <- read.table(paste0(out.dir, "CEU_results_all_info.txt" ), header = TRUE, sep = "\t", as.is = TRUE)

pdf(paste0(out.dir.plots, "hist_pvalues.pdf"))

hist(res.table[, "PValue"], col = "orange", breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "")

dev.off()




##########################################################################################
####### Plot number of SNPs per gene
##########################################################################################


res.all<- read.table(paste0(out.dir, "CEU_results_all_info.txt" ), header = TRUE, sep = "\t", as.is = TRUE)


pdf(paste0(out.dir.plots, "/Hist_numberOfSNPsPerGene.pdf"))
tt <- table(res.all$gene_id)

hist(tt, breaks = 100, col = "darkseagreen2", main = paste0("All ",length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "Number of SNPs per gene")

tt <- table(res.all[res.all$FDR < 0.05, "gene_id"])

hist(tt, breaks = 100, col = "darkturquoise", main = paste0( "Significant ",length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "Number of SNPs per gene")

dev.off()


##########################################################################################
####### Check which gene has the most of the significant snps
##########################################################################################


tt <- table(res.all$gene_id)
ttSign <- table(res.all[res.all$FDR < 0.05, "gene_id"])

ttSign <- sort(ttSign, decreasing = TRUE)

topNr <- length(ttSign)

df <- data.frame(gene_id = names(ttSign[1:topNr]), NumberSignSNPs = ttSign[1:topNr], NumberAllSNPs = tt[names(ttSign[1:topNr])])

df$Percent <- df$NumberSignSNPs / df$NumberAllSNPs

pdf(paste0(out.dir.plots, "/Percent_numberOfSignificantSNPsPerGene.pdf"))
ggp <- ggplot(df, aes(x = log10(NumberAllSNPs), y = Percent)) + geom_point(alpha = 0.2)
print(ggp)
ggp <- ggplot(df, aes(x = log10(NumberAllSNPs), y = log10(NumberSignSNPs))) + 
  geom_point(alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1, colour = "red")
print(ggp)
dev.off()








##########################################################################################
# check if snps from sQTLseekeR paper are significant 
##########################################################################################

####### check the significans of snps from sQTLseekeR paper


snp <- "snp_19_41937095"
gene <- "ENSG00000105341.11"

res.all[res.all$SNP_id == snp & res.all$gene_id == gene, ]


snp <- "snp_5_96244549"
gene <- "ENSG00000164308.12"

res.all[res.all$SNP_id == snp & res.all$gene_id == gene, ]


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
    
    main <- paste0(gene, " - ", snp, "\n meanExpression = ", round(dgeSQTL$meanExpr[gene]), " / TagwiseDispersion = ", round(dgeSQTL$tagwiseDispersion[snpgene], 2), "\n LR = ", round(dgeSQTL$table[index.table, "LR"], 4) , " / PValue = ", sprintf("%.02e", dgeSQTL$table[index.table, "PValue"]), " / FDR = ", sprintf("%.02e", dgeSQTL$table[index.table, "FDR"]))
    
    
    ### box plots with points
    ggb <- ggplot(prop.smp.m, aes(x = ete_id, y = Proportions)) +
      theme_bw() + 
      # theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), legend.position="none", plot.title = element_text(size=10)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), plot.title = element_text(size=10)) +
      ggtitle(main) +     
      geom_jitter(aes(fill = Condition, colour = Condition), position = position_jitterdodge(dodge.width = 0.75), alpha = 0.5, show_guide = FALSE) +
      geom_boxplot(aes(colour = Condition), fill = "white", outlier.size = NA, alpha = 0) + 
      geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, fill = Condition), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 3, shape = 23, show_guide = FALSE) +
      geom_point(data = prop.est.null, aes(x = ete_id, y = Proportions), size = 3, shape = 23, fill = "orange") +
      coord_cartesian(ylim = c(-0.1, 1.1)) +
    scale_x_discrete(name="Transcripts") +
      scale_colour_discrete(name = "Variants")
    
    print(ggb)
    
    ### box plots per condition
    ggb <- ggplot(prop.smp.m, aes(x = Condition, y = Proportions)) +
      theme(axis.text.x = element_text(angle = 0, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), plot.title = element_text(size=10), panel.grid.major=element_blank()) +
      geom_vline(xintercept=c(1.5,2.5),color="white") +
      ggtitle(main) +     
      scale_x_discrete(labels = paste0( "Variant ", names(Condition.counts), " (",Condition.counts, ")" ), name="") +
      geom_boxplot(aes(fill = ete_id), width = 1) + 
      geom_point(data = prop.est.m, aes(x = Condition, y = Proportions, fill = ete_id), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 2, shape = 23, colour = "black") +
      coord_cartesian(ylim = c(-0.1, 1.1)) +
      geom_vline(xintercept=c(1.5,2.5), color="white")+
      scale_fill_discrete(name = "Transcripts")
    
    print(ggb)
    
    
  }
  
  dev.off()
  
}





#### SNP with most negative LR in DM_0.1.2 analysis
gene <- "ENSG00000169045.13"
snp <- "snp_5_179056159"

chr = 5
load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))

plot.snps <- data.frame(gene_id = gene, SNP_id = snp, stringsAsFactors = FALSE)

plotPath <- paste0(out.dir.plots, "Proportions.pdf")


plotProportions(dgeSQTL, plot.snps, plotPath)



#### SNPs from sQTLseekeR paper

snp <- "snp_5_96244549"
gene <- "ENSG00000164308.12"


plot.snps <- data.frame(gene_id = gene, SNP_id = snp, stringsAsFactors = FALSE)

plotPath <- paste0(out.dir.plots, "Proportions_SNPs_from_Seeker_paper",snp,".pdf")


plotProportions(dgeSQTL, plot.snps, plotPath)



snp <- "snp_19_41937095"
gene <- "ENSG00000105341.11"

chr = 19
load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))


plot.snps <- data.frame(gene_id = gene, SNP_id = snp, stringsAsFactors = FALSE)

plotPath <- paste0(out.dir.plots, "Proportions_SNPs_from_Seeker_paper",snp,".pdf")


plotProportions(dgeSQTL, plot.snps, plotPath)






#### most significant SNPs that belong to a gene with highest number of siginificant SNPs

gene <- "ENSG00000196735.6"

snps <- res.all[res.all$gene_id == gene & res.all$FDR < 0.05, "SNP_id"]

chr = 6
load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))



plot.snps <- data.frame(gene_id = gene, SNP_id = snps, stringsAsFactors = FALSE)

plotPath <- paste0(out.dir.plots, "Proportions_Gene_highest_nr_sign_SNPs_head.pdf")


plotProportions(dgeSQTL, head(plot.snps, 10), plotPath)




plotPath <- paste0(out.dir.plots, "Proportions_Gene_highest_nr_sign_SNPs_tail.pdf")


plotProportions(dgeSQTL, tail(plot.snps, 10), plotPath)























