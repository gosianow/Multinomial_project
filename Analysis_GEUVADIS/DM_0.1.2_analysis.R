# BioC 3.0

# Created 16 Jan 2014
# Modyfied 16 Jan 2014

# Analyse data from DM_0.1.1_filtering

##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)


Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)



##########################################################################################
# Run the DMsQTL pipeline by chromosome // df <- DFnull * (length(fit.full[[snp]]$df) - 1)
##########################################################################################

out.dir <- "DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)



load(paste0("DM_0.1.1_sQTL_analysis/Data/dgeSQTL.RData"))

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

out.dir <- "DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/"
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




####### plot p-values distribution

out.dir <- "DM_0.1.2_sQTL_analysis/Plots_TagwiseDisp_gridCommonDispersion/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


pdf(paste0(out.dir, "hist_pvalues.pdf"))

hist(res[, "PValue"], col = 1, breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "")
  
dev.off()





####### check the significans of snps from sQTLseekeR paper


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












































