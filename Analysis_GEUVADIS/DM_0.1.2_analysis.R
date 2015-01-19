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
# Run the DMsQTL pipeline by chromosome 
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




######### 



snp <- "snp_19_41937095"

snp <- "snp_19_41936608"

gene <- "ENSG00000105341.11"


dgeSQTL$counts[[gene]]

dgeSQTL$fit[[snp]]



snps <- which(names(dgeSQTL$fit) == snp)

dgeSQTL$fit[snps]
dgeSQTL$fit.null[snps]

dgeSQTL$SNPs[snps, ]


dgeSQTL$tagwiseDispersion[snps]


snp <- snps[2]



sum(duplicated(dgeSQTL$SNPs$SNP_id))

























