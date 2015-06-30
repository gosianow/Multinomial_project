# BioC 3.0

# Created 2 June 2015
# Updated 10 June 2015
# Update 16 June 2015
# - compare the Data_DM_0_1_5 results

############################################################################################


setwd("/home/Shared/data/seq/GEUVADIS/")

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(sQTLseekeR)

library(ggplot2)
library(reshape2)



out.dir.plots <- "DM_0_1_4_sQTL_analysis/Comparison_Data_DM_0_1_2_clean_TagwiseDisp_gridNone_tol12_constrOptim2G/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)


out.dir <- "DM_0_1_4_sQTL_analysis/Results_Data_DM_0_1_5_TagwiseDisp_gridNone_tol12_constrOptim2G/"
out.dir.plots <- "DM_0_1_4_sQTL_analysis/Comparison_Data_DM_0_1_5_TagwiseDisp_gridNone_tol12_constrOptim2G/"
dir.create(out.dir.plots, showWarnings = FALSE, recursive = TRUE)



##############################################################
# Load results
##############################################################

################### Load seekeR


res.seeker.all <- read.table("sQTLseekeR20_analysis/Results/CEU_results_all.txt", header = TRUE, as.is = TRUE)
res.seeker <- read.table("sQTLseekeR20_analysis/Results/CEU_results_FDR05.txt", header = TRUE, as.is = TRUE)


################### Load DM results


res.dm.all <- read.table("DM_0_1_4_sQTL_analysis/Results_Data_DM_0_1_2_clean_TagwiseDisp_gridNone_tol12_constrOptim2G/CEU_results_all_info.txt", header = TRUE, as.is = TRUE)
res.dm <- res.dm.all[res.dm.all$FDR < 0.05, ]



res.dm.all <- read.table("DM_0_1_4_sQTL_analysis/Results_Data_DM_0_1_5_TagwiseDisp_gridNone_tol12_constrOptim2G/CEU_results_all_info.txt", header = TRUE, as.is = TRUE)
res.dm <- res.dm.all[res.dm.all$FDR < 0.05, ]



#######################################################
# generate hist of p-values
#######################################################


pdf(paste0(out.dir.plots, "hist_pvalues_sQTLseekeR.pdf"))

hist(res.seeker.all[, "pv"], col = "#1874CD" , breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "sQTLseekeR", cex.main = 3)

dev.off()

pdf(paste0(out.dir.plots, "hist_pvalues_DM.pdf"))

hist(res.dm.all[, "PValue"], col = "#FF7F00", breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "DM", cex.main = 3)

dev.off()


#######################################################
# generate venn diagrams 
#######################################################


source("/home/gosia/R/R_Multinomial_project/Plot_Functions/plotVenn_0_1_4.R")

#### gene-snp pairs
venn_list <- list()

venn_list[["sQTLseeker"]] <- paste0(res.seeker[,"snpId"], "-", res.seeker[,"geneId"])
venn_list[["DM"]] <- paste0(res.dm[,"SNP_id"], "-", res.dm[, "gene_id"])

plotPath <- paste0(out.dir.plots, "Venn_Diagr_", "snps",".pdf")

plotVenn2(venn_list, venn.colors = c(sQTLseeker= "#1874CD", DM= "#FF7F00" ), venn.subset = names(venn_list), margin = 0, cat.cex=1.8, cex=1.7, plotPath = plotPath)


#### genes
venn_list <- list()

venn_list[["sQTLseeker"]] <- unique(res.seeker[,"geneId"])
venn_list[["DM"]] <- unique(res.dm[,"gene_id"])

plotPath <- paste0(out.dir.plots, "Venn_Diagr_", "genes",".pdf")

plotVenn2(venn_list, venn.colors = c(sQTLseeker= "#1874CD", DM= "#FF7F00" ), venn.subset = names(venn_list), margin = 0, cat.cex=1.8, cex=1.7, plotPath = plotPath)




#######################################################
# generate sets of snps for plotting - FP for DM 
#######################################################

res.seeker$SNPgene <- paste0(res.seeker[,"snpId"], "-", res.seeker[,"geneId"])
res.dm$SNPgene <- paste0(res.dm[,"SNP_id"], "-", res.dm[, "gene_id"])


fp.dm <- setdiff(res.dm$SNPgene, res.seeker$SNPgene)
length(fp.dm)

res.dm.fp <- res.dm[res.dm$SNPgene %in% fp.dm, ]


pdf(paste0(out.dir.plots, "/Hist_numberOfSNPsPerGene_FP.pdf"))
tt <- table(res.dm.fp$gene_id)
hist(tt, breaks = 100, col = "brown", main = paste0("All ",length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "Number of SNPs per gene")
dev.off()


res.dm.fp <- res.dm.fp[order(res.dm.fp$PValue, decreasing = FALSE), ]
res.dm.fp.top <- split(res.dm.fp, f = factor(res.dm.fp$gene_id, levels = unique(res.dm.fp$gene_id)))
res.dm.fp.top <-  lapply(res.dm.fp.top, function(g){ 
  nrSNPS <- nrow(g)
  g <- g[1, , drop = FALSE]
  g$nrSNPS <- nrSNPS
  return(g)
  })
res.dm.fp.top <- do.call(rbind, res.dm.fp.top)

dim(res.dm.fp.top)

max(res.dm.fp.top$nrSNPS)

res.dm.fp.top$chr <- strsplit2(res.dm.fp.top$SNP_id, "_")[, 2]


#######################################################
# generate sets of snps for plotting - FN for DM 
#######################################################

res.seeker$SNPgene <- paste0(res.seeker[,"snpId"], "-", res.seeker[,"geneId"])
res.dm$SNPgene <- paste0(res.dm[,"SNP_id"], "-", res.dm[, "gene_id"])
res.dm.all$SNPgene <- paste0(res.dm.all[,"SNP_id"], "-", res.dm.all[, "gene_id"])

fn.dm <- setdiff(res.seeker$SNPgene, res.dm$SNPgene)

res.dm.fn <- res.dm.all[res.dm.all$SNPgene %in% fn.dm, ]

res.dm.fn <- merge(res.dm.fn, res.seeker[, c("geneId", "pv", "qv")], by.x = "gene_id", by.y = "geneId", all.x = TRUE, sort = FALSE)


### keep top snp from each gene
res.dm.fn <- res.dm.fn[order(res.dm.fn$qv, decreasing = FALSE), ]
res.dm.fn <- split(res.dm.fn, f = factor(res.dm.fn$gene_id, levels = unique(res.dm.fn$gene_id)))
res.dm.fn.top <-  lapply(res.dm.fn, function(g){ 
  g <- g[1, , drop = FALSE]
  return(g)
  })

res.dm.fn.top <- do.call(rbind, res.dm.fn.top)


res.dm.fn.top <- res.dm.fn.top[order(res.dm.fn.top$PValue, decreasing = TRUE), ]
res.dm.fn.top <- res.dm.fn.top[order(res.dm.fn.top$qv, decreasing = FALSE), ]


dim(res.dm.fn.top)


res.dm.fn.top$chr <- strsplit2(res.dm.fn.top$SNP_id, "_")[, 2]


##############################################################
# Plot proportions of FN and FP
##############################################################

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
    
#     ### box plots per condition
#     ggb <- ggplot(prop.smp.m, aes(x = Condition, y = Proportions)) +
#       theme(axis.text.x = element_text(angle = 0, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), plot.title = element_text(size=10), panel.grid.major=element_blank()) +
#       geom_vline(xintercept=c(1.5,2.5),color="white") +
#       ggtitle(main) +     
#       scale_x_discrete(labels = paste0( "Variant ", names(Condition.counts), " (",Condition.counts, ")" ), name="") +
#       geom_boxplot(aes(fill = ete_id), width = 1) + 
#       geom_point(data = prop.est.m, aes(x = Condition, y = Proportions, fill = ete_id), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 2, shape = 23, colour = "black") +
#       coord_cartesian(ylim = c(-0.1, 1.1)) +
#       geom_vline(xintercept=c(1.5,2.5), color="white")+
#       scale_fill_discrete(name = "Transcripts")
#     
#     print(ggb)
    
    
  }
  
  dev.off()
  
}







########## plot FN
snps2plot <- res.dm.fn.top


for(i in 1:10){
  
  snp <- snps2plot$SNP_id[i]
  gene <- snps2plot$gene_id[i]
  
  chr = snps2plot$chr[i]
  
  load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))
  
  
  plot.snps <- data.frame(gene_id = gene, SNP_id = snp, stringsAsFactors = FALSE)
  
  plotPath <- paste0(out.dir.plots, "Propportions_FN_",i,".pdf")
  
  plotProportions(dgeSQTL, plot.snps, plotPath)
  
}




########## plot FP
snps2plot <- res.dm.fp.top


for(i in 1:10){
  
  snp <- snps2plot$SNP_id[i]
  gene <- snps2plot$gene_id[i]
  nrSNPS <- snps2plot$nrSNPS[i]
  chr = snps2plot$chr[i]
  
  load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))
  
  
  plot.snps <- data.frame(gene_id = gene, SNP_id = snp, stringsAsFactors = FALSE)
  
  plotPath <- paste0(out.dir.plots, "Propportions_FP_",i,"_nrSNPS_",nrSNPS,".pdf")
  
  plotProportions(dgeSQTL, plot.snps, plotPath)
  
}






##############################################################
# Plot proportions of genes that are filtered out in DM_0_1_5_Data
##############################################################


load("DM_0_1_5_Data/counts.RData")

head(counts)

genes5 <- names(counts)

res.dm.filtr <- res.dm.all[!res.dm.all$gene_id %in% genes5, ]


res.dm.filtr <- res.dm.filtr[order(res.dm.filtr$PValue, decreasing = FALSE), ]
res.dm.filtr.top <- split(res.dm.filtr, f = factor(res.dm.filtr$gene_id, levels = unique(res.dm.filtr$gene_id)))
res.dm.filtr.top <-  lapply(res.dm.filtr.top, function(g){ 
  nrSNPS <- nrow(g)
  g <- g[1, , drop = FALSE]
  g$nrSNPS <- nrSNPS
  return(g)
})
res.dm.filtr.top <- do.call(rbind, res.dm.filtr.top)

dim(res.dm.filtr.top)

max(res.dm.filtr.top$nrSNPS)

res.dm.filtr.top$chr <- strsplit2(res.dm.filtr.top$SNP_id, "_")[, 2]




########## plot filtered snps
snps2plot <- res.dm.filtr.top


for(i in 1:10){
  
  snp <- snps2plot$SNP_id[i]
  gene <- snps2plot$gene_id[i]
  nrSNPS <- snps2plot$nrSNPS[i]
  chr = snps2plot$chr[i]
  
  load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))
  
  
  plot.snps <- data.frame(gene_id = gene, SNP_id = snp, stringsAsFactors = FALSE)
  
  plotPath <- paste0(out.dir.plots, "Propportions_Filter_",i,"_nrSNPS_",nrSNPS,".pdf")
  
  plotProportions(dgeSQTL, plot.snps, plotPath)
  
}


########## some checks

i = 1

snp <- snps2plot$SNP_id[i]
gene <- snps2plot$gene_id[i]
nrSNPS <- snps2plot$nrSNPS[i]
chr = snps2plot$chr[i]

load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))

dgeSQTL$counts[[gene]]

colSums(dgeSQTL$counts[[gene]])




gene <- "ENSG00000184674.7"
chr = 22

load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))

dgeSQTL$counts[[gene]]

colSums(dgeSQTL$counts[[gene]])



























