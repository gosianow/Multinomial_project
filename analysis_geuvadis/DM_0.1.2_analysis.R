# BioC 3.0

# Created 16 Jan 2014
# Analyse data from DM_0.1.2_filtering

# Modyfied 9 Feb 2014
# Add plots of mean expression vs dispersion 


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


# Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
# for(i in Rfiles) source(i)



##########################################################################################
# Run the DMsQTL pipeline by chromosome // df <- DFnull * (length(fit.full[[snp]]$df) - 1)
##########################################################################################

out.dir <- "DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)



load(paste0("DM_0.1.2_sQTL_analysis/Data/dgeSQTL.RData"))

dgeSQTL.org <- dgeSQTL

head(dgeSQTL.org$SNPs$SNP_id)
head(dgeSQTL.org$genes$gene_id)

keep.genes <- data.frame(gene_id = intersect(unique(dgeSQTL.org$genes$gene_id), unique(dgeSQTL.org$SNPs$gene_id)), stringsAsFactors = FALSE)
dim(keep.genes)

keep.genes <- merge(keep.genes, unique(dgeSQTL.org$SNPs[, c("gene_id", "chr")]), by = "gene_id", all.x = TRUE, sort = FALSE)
dim(keep.genes)


mcCores <- 25




for(chr in 22:1){
  # chr = 19
  
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
  
  dgeSQTL <- dmSQTLEstimateCommonDisp(dgeSQTL, adjust = TRUE, subset = round(nrow(dgeSQTL$genotypes)/170) , mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+2), tol = 1e-02, mcCores=mcCores, verbose=FALSE)

  write.table(dgeSQTL$commonDispersion, paste0(out.dir, "commonDispersion_chr", chr, ".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#   dgeSQTL$commonDispersion <- read.table(paste0(out.dir, "commonDispersion_chr", chr, ".txt" ))
  
  
  ######### tagwiseDispersion
  cat(paste0("chr", chr," estimate tagwise dispersion \n"))
  
  dgeSQTL <- dmSQTLEstimateTagwiseDisp(dgeSQTL, adjust = TRUE, mode = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[3], epsilon = 1e-05, maxIte = 1000, modeDisp = c("optimize", "optim", "constrOptim", "grid")[4], interval = c(0, 50), tol = 1e-01,  initDisp = 2, initWeirMoM = TRUE, gridLength = 10, gridRange = c(-6, 3), trend = c("none", "commonDispersion", "trendedDispersion")[2], priorDf = 10, span = 0.3, mcCores = mcCores, verbose = FALSE, plot = FALSE)
  
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
  
  dgeSQTL <- dmSQTLTest(dgeSQTL, dispersion=c("commonDispersion", "tagwiseDispersion")[2], mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
  
  save(dgeSQTL, file = paste0(out.dir, "dgeSQTL_chr",chr,".RData"))
  
  write.table(dgeSQTL$table, paste0(out.dir, "results_chr", chr, ".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat(paste0("chr", chr," done \n \n"))
  
  gc()
  
}



### check how many snps are called significant 
dim(dgeSQTL$table)
sum(dgeSQTL$table$FDR < 0.05, na.rm = TRUE)


### check how many genes are called significant 
length(unique((dgeSQTL$table$gene_id)))
length(unique((dgeSQTL$table$gene_id[dgeSQTL$table$FDR < 0.05])))



##########################################################################################
# Run the DMsQTL pipeline by chromosome 
# with degrees of freedom df <- DFfull - DFnull
##########################################################################################


setwd("/home/Shared/data/seq/GEUVADIS/")

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)


Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)



out.dir <- "DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion_df/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


mcCores <- 10


for(chr in 16:1){
  # chr = 19
  
  load(paste0("DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/", "dgeSQTL_chr",chr,".RData"))
  
  ######### LR testing
  cat(paste0("chr", chr," testing \n"))
  
  fit.full <- dgeSQTL$fit
  fit.null <- dgeSQTL$fit.null
  
  cat("Calculating LR \n")
  LRT <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){
    # snp = 1
    
    #     if(verbose) 
    # cat("testing SNP: ", paste0(dgeSQTL$SNPs[snp, c("SNP_id", "gene_id")], collapse = "-"), fill = TRUE)
    
    if(is.null(fit.null[[snp]]) || is.null(fit.full[[snp]])) 
      return(data.frame(LR=NA, df=NA, PValue=NA, LLfull=NA, LLnull=NA))
    
    LLnull <- fit.null[[snp]]$logLik
    
    LLfull <- sum(fit.full[[snp]]$logLik)
    
    LR <-  2*(LLfull - LLnull)
    
    DFnull <- fit.null[[snp]]$df
    DFfull <- sum(fit.full[[snp]]$df)
    
    df <- DFfull - DFnull
    # df <- DFnull * (length(fit.full[[snp]]$df) - 1)
    
    pValue <- pchisq(LR, df = df , lower.tail = FALSE)
    
    #     gc()
    
    return(data.frame(LR=LR, df=df, PValue=pValue, LLfull=LLfull, LLnull=LLnull))
    
  }, mc.cores=mcCores)
  
  # save(fit.null, LRT, file="LRT.RData")
  
  LRT <- do.call(rbind, LRT)
  FDR <- p.adjust(LRT$PValue, method="BH")
  o <- order(LRT$PValue)
  
  table <- data.frame(SNP_id = dgeSQTL$SNPs$SNP_id, gene_id = dgeSQTL$SNPs$gene_id, LR=LRT$LR, df=LRT$df, LLfull=LRT$LLfull, LLnull=LRT$LLnull , PValue=LRT$PValue, FDR=FDR, stringsAsFactors = FALSE)[o,]
  
  write.table(table, paste0(out.dir, "results_chr", chr, ".txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat(paste0("chr", chr," done \n \n"))
  
  gc()
  
}




##########################################################################################
# Merge DM tables with results and recalculate FDR 
# Merge tagwise dispersion 
# check if snps from sQTLseekeR paper are significant 
##########################################################################################

out.dir <- "DM_0_1_2_sQTL_analysis/Results_Data_DM_0_1_2_TagwiseDisp_gridCommonDispersion/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

##### merge tables with results

res <- lapply(22:1, function(chr){ results <- read.table(paste0(out.dir, "results_chr", chr, ".txt" ), header = TRUE, as.is = TRUE)})
  

res <- do.call(rbind, res)


FDR <- p.adjust(res$PValue, method="BH")

res$FDR <- FDR

write.table(res, paste0(out.dir, "CEU_results_all.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(res[res$FDR < 0.05 & !is.na(res$FDR), ], paste0(out.dir, "CEU_results_FDR05.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




##### merge tagwise dispersion 

res <- lapply(22:1, function(chr){ results <- read.table(paste0(out.dir, "tagwiseDispersion_chr", chr, ".txt" ), header = TRUE, as.is = TRUE)})

res <- do.call(rbind, res)

write.table(res, paste0(out.dir, "tagwiseDispersion.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##### merge mean gene expression & nr of transcripts

res <- lapply(22:1, function(chr){ 
  # chr = 22
  
  cat(paste0("chr", chr," \n"))
  
  load(paste0(out.dir, "dgeSQTL_chr", chr, ".RData" ))
  
  meanExpr <- merge(dgeSQTL$SNPs[, c("SNP_id", "gene_id")], data.frame(gene_id = names(dgeSQTL$meanExpr), meanExpr = dgeSQTL$meanExpr, stringsAsFactors = FALSE), by = "gene_id", all.x = TRUE, sort = FALSE)
  
  meanExpr$tagwiseDispersion <- dgeSQTL$tagwiseDispersion
  
  meanExpr$SNPgene <- paste0(meanExpr$SNP_id, "-", meanExpr$gene_id)

  nrTrans <- t(sapply(dgeSQTL$counts, dim, USE.NAMES = TRUE))
  
  nrTrans <- data.frame(gene_id = rownames(nrTrans), nrTrans = nrTrans[, 1], stringsAsFactors = FALSE, row.names = NULL )
  
  meanExpr <- merge(meanExpr, nrTrans, by = "gene_id", all.x = TRUE, sort = FALSE)
  
  return(meanExpr)
  
})

res <- do.call(rbind, res)

write.table(res, paste0(out.dir, "meanExpr.txt" ), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




### check how many snps are called significant 
dim(res)
sum(res$FDR < 0.05, na.rm = TRUE)

### check how many genes are called significant 
length(unique((res$gene_id)))
length(unique((res$gene_id[res$FDR < 0.05])))



##########################################################################################
####### plot p-values distribution
##########################################################################################
out.dir.plots <- "DM_0_1_2_sQTL_analysis/Plots_Data_DM_0_1_2_TagwiseDisp_gridCommonDispersion/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


pdf(paste0(out.dir.plots, "Hist_pvalues.pdf"))

hist(res[, "PValue"], col = 1, breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "")
  
dev.off()


##########################################################################################
####### Plot number of SNPs per gene
##########################################################################################


pdf(paste0(out.dir.plots, "/Hist_numberOfSNPsPerGene.pdf"))
tt <- table(res$gene_id)

hist(tt, breaks = 100, col = "darkseagreen2", main = paste0("All ",length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "Number of SNPs per gene")

tt <- table(res[res$FDR < 0.05, "gene_id"])

hist(tt, breaks = 100, col = "darkturquoise", main = paste0( "Significant ",length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "Number of SNPs per gene")

dev.off()


##########################################################################################
####### check the significans of snps from sQTLseekeR paper
##########################################################################################

snp <- "snp_19_41937095"
gene <- "ENSG00000105341.11"

res[res$SNP_id == snp & res$gene_id == gene, ]


snp <- "snp_5_96244549"
gene <- "ENSG00000164308.12"

res[res$SNP_id == snp & res$gene_id == gene, ]

chr = 5

load(paste0("DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/", "dgeSQTL_chr",chr,".RData"))


snpn <- which(dgeSQTL$SNPs$SNP_id == snp & dgeSQTL$SNPs$gene_id == gene)

dgeSQTL$SNPs[snpn, ]


dgeSQTL$fit[[snpn]]

expr <- dgeSQTL$counts[[gene]]

piH <- dgeSQTL$fit[[snpn]]$piH

rownames(piH) <- rownames(expr)





##########################################################################################
# Merge proportions piH from dgeSQTL
##########################################################################################

out.dir <- "DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

ngroups <- 3
lgroups <- c(0, 1, 2)

res <- mclapply(22:1, function(chr){ 
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
# Plot dispersion vs mean expression 
##########################################################################################

out.dir <- "DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


res <- mclapply(22:1, function(chr){ 
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









#######################################################
# !!!! Investigate the DM hist of p-values !!!!
#######################################################
res.dm.all <- read.table("DM_0_1_2_sQTL_analysis/Results_Data_DM_0_1_2_TagwiseDisp_gridCommonDispersion/CEU_results_all.txt", header = TRUE, as.is = TRUE)


head(res.dm.all)

dim(res.dm.all)

res.dm.all <- res.dm.all[complete.cases(res.dm.all), ]


table(res.dm.all$PValue == 1)

table(res.dm.all$PValue == 1 & res.dm.all$LR  < 0)

table(res.dm.all$LR < 0)

mean(res.dm.all$LR < 0)


min(res.dm.all$LR)
res.dm.all[which.min(res.dm.all$LR), ]



min(res.dm.all[res.dm.all$LR < 0, "df"])


pdf(paste0(out.dir, "hist_df.pdf"))
hist(res.dm.all[, "df"], breaks = 100)
hist(res.dm.all[res.dm.all$LR < 0, "df"], breaks = 100)
dev.off()



png(paste0(out.dir, "LR_LLnull.png"), 700, 700)
smoothScatter(log10(-res.dm.all$LLnull), res.dm.all$LR, nrpoints = Inf)
abline(h = 0, col = 2)
dev.off()


pdf(paste0(out.dir, "hist_LR.pdf"))
hist(res.dm.all$LR, breaks = 500)
dev.off()


pdf(paste0(out.dir, "hist_pvalues_DM.pdf"))
hist(res.dm.all[, "PValue"], col = "#FF7F00", breaks = 100, cex.lab=1.5, cex.axis = 1.5, main = "DM", cex.main = 3)
hist(res.dm.all[res.dm.all$LR > 0, "PValue"], col = "#FF7F00", breaks = 100, cex.lab=1.5, cex.axis = 1.5, main = "DM", cex.main = 3)
dev.off()



### results for genes with p-value = 1
res.pv1 <- res.dm.all[res.dm.all$PValue == 1, ]

table(table(res.pv1$gene_id))
### gene that has most of the p-values = 1
table(res.pv1$gene_id)[table(res.pv1$gene_id) == 7929] ## ENSG00000160183.8


### results for this one gene
res1g <- res.dm.all[res.dm.all$gene_id == "ENSG00000160183.8", ]
dim(res1g)

res1g <- res.pv1[res.pv1$gene_id == "ENSG00000160183.8" & res.pv1$LR < 0, ]
dim(res1g)
head(res1g)




### get the cases with 2 transripts for plotting the piH likelihoods
mex <- read.table("DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/meanExpr.txt", header = TRUE, sep = "\t")

g1tr <- mex[mex$nrTrans == 2, ]

g0lr <- res.dm.all[res.dm.all$LR < 0, ]

sum(unique(g0lr$gene_id) %in% unique(g1tr$gene_id))

g1tr0lr <- g0lr[g0lr$gene_id %in% g1tr$gene_id, ]

min(g1tr0lr$LR)



### plot the transcript expressions and ratios

out.dir <- "DM_0.1.2_sQTL_analysis/Plots_TagwiseDisp_gridCommonDispersion/DM_PValuesEqual1/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

### have to set up this values

### some results from top of the table
# plot.snps <- res1g[, c("gene_id", "SNP_id")]
# plot.names <- paste0("nr", 1:nrow(res1g))
# plot.main <- paste0("LR ", res1g[, c("LR")])

### genes with 2 transcripts
plot.snps <- g1tr0lr[, c("gene_id", "SNP_id")]
plot.names <- paste0("_2tr_nr", 1:nrow(g1tr0lr))
plot.main <- paste0("LR ", g1tr0lr[, c("LR")])


### SNP with most negative LR
minNeg <- res.dm.all[which.min(res.dm.all$LR),, drop= FALSE]
plot.snps <- minNeg[,  c("gene_id", "SNP_id")]
plot.names <- paste0("_minNeg_nr", 1:nrow(minNeg))
plot.main <- paste0("LR ", minNeg[, c("LR")])



for(i in 1:5){
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




##### load results from DM_v5

res5 <- read.table("DMv5_sQTL_analysis/dgeSQTL_results.txt", header = TRUE, sep = "\t")

table(res5$LR < 0)




### read in the complete results from DM_0.1.2


load("DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/dgeSQTL_chr21.RData")
load("DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/dgeSQTL_chr19.RData")

gene <- plot.snps[1,1]
snp <- plot.snps[1,2]


### the SNP that have the most negative LR - snp_5_179056159-ENSG00000169045.13
load("DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/dgeSQTL_chr5.RData")


gene <- plot.snps[1,1]
snp <- plot.snps[1,2]


### display all the fitting results 

dgeSQTL$fit[[snp]]

colSums(dgeSQTL$fit[[snp]]$piH)


dgeSQTL$fit.null[[snp]]

sort(as.numeric(dgeSQTL$fit.null[[snp]]$piH), decreasing = TRUE)[1:20]
sort(as.numeric(dgeSQTL$fit[[snp]]$piH[, 1]), decreasing = TRUE)[1:20]


dgeSQTL$meanExpr[gene]


dgeSQTL$counts[[gene]]

dim(dgeSQTL$counts[[gene]])






y <- dgeSQTL$counts[[gene]]
gamma0 <- dgeSQTL$dispersion[paste0(snp, "-", gene)]



fitDM <- DM::dmOneGeneGroup(y, gamma0, mode = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")[5], epsilon = 1e-05, maxIte = 1000, verbose = FALSE, plot = FALSE)



































