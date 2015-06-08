# BioC 3.1

# Created 2 June 2015
# Analyse data from DM_0.1.5_filtering with DM_0.1.5 (New dgeSQTL object and use of BiocParallel)


##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

library(DM)
library(BiocParallel)

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
# Run the DMsQTL pipeline
##########################################################################################


######### run on DM_0_1_5 data 
out.dir <- "DM_0_1_5_sQTL_analysis/Results_Data_DM_0_1_5_TagwiseDisp_gridNone/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

load(paste0("DM_0_1_5_Data/dgeSQTL.RData"))


dgeSQTL_org <- dgeSQTL

### subset 100 genes
geneList <- names(dgeSQTL_org$counts)[1:20]
dgeSQTL$counts <- dgeSQTL_org$counts[geneList]
dgeSQTL$genotypes <- dgeSQTL_org$genotypes[geneList]

dgeSQTL$commonDispersion <- 4

BPPARAM <- MulticoreParam(workers = 5)

# BPPARAM <- SnowParam(workers = 5)


######### commonDispersion
cat(paste0("Estimate common dispersion \n"))

dgeSQTL <- dmSQTLEstimateCommonDisp(dgeSQTL, adjustDisp = TRUE, intervalDisp = c(0, 1e+3), tolDisp = 1e-01, modeProp = "constrOptim2G", tolProp = 1e-12, verbose=FALSE, BPPARAM = BPPARAM)


dgeSQTL$commonDispersion <- 4


######### tagwiseDispersion
cat(paste0("Estimate tagwise dispersion \n"))

dgeSQTL <- dmSQTLEstimateTagwiseDisp(dgeSQTL, adjustDisp = TRUE, modeDisp = c("optimize", "optim", "constrOptim", "grid")[2], intervalDisp = c(0, 1e+5), tolDisp = 1e-08,  initDisp = 10, initWeirMoMDisp = TRUE, gridLengthDisp = 15, gridRangeDisp = c(-7, 7), trendDisp = c("none", "commonDispersion", "trendedDispersion")[1], priorDfDisp = 10, spanDisp = 0.3, modeProp = c( "constrOptim2", "constrOptim2G", "FisherScoring")[2], tolProp = 1e-12, verbose = FALSE, plot = FALSE, BPPARAM = BPPARAM)



######### DM fitting 
cat(paste0("chr", chr," fitting DM \n"))

dgeSQTL <- dmSQTLFit(dgeSQTL, model="full", dispersion=c("commonDispersion", "tagwiseDispersion")[2], mode=c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores )


######### LR testing
cat(paste0("chr", chr," testing \n"))

dgeSQTL <- dmSQTLTest(dgeSQTL, dispersion=c("commonDispersion", "tagwiseDispersion")[2], mode=c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

save(dgeSQTL, file = paste0(out.dir, "dgeSQTL_chr",chr,".RData"))

write.table(dgeSQTL$table, paste0(out.dir, "results_chr", chr, ".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




##########################################################################################
#### Choose the results to work with 
##########################################################################################


out.dir <- "DM_0_1_4_sQTL_analysis/Results_Data_DM_0_1_2_clean_TagwiseDisp_gridNone_tol12_constrOptim2G/"
out.dir.plots <- "DM_0_1_4_sQTL_analysis/Plots_Data_DM_0_1_2_clean_TagwiseDisp_gridNone_tol12_constrOptim2G/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)




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
  
  gc()
  
  return(NULL)
  
}, mc.cores = 5)




##########################################################################################
####### Plot p-values distribution
##########################################################################################



pdf(paste0(out.dir, "hist_pvalues.pdf"))

hist(res.table[, "PValue"], col = "orange", breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "")

dev.off()




##########################################################################################
####### Plot number of SNPs per gene
##########################################################################################




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




gene <- "ENSG00000196735.6"

snps <- res.all[res.all$gene_id == gene & res.all$FDR < 0.05, "SNP_id"]

chr = 6
load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))






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





#### SNP with most genative LR in DM_0.1.2 analysis
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

chr = 6
load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))

plot.snps <- data.frame(gene_id = gene, SNP_id = snps, stringsAsFactors = FALSE)

plotPath <- paste0(out.dir.plots, "Proportions_Gene_highest_nr_sign_SNPs_head.pdf")


plotProportions(dgeSQTL, head(plot.snps, 10), plotPath)




plotPath <- paste0(out.dir.plots, "Proportions_Gene_highest_nr_sign_SNPs_tail.pdf")


plotProportions(dgeSQTL, tail(plot.snps, 10), plotPath)























