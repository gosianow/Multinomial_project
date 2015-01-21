# BioC 3.0

# Created 14 Jan 2014
# Modyfied 15 Jan 2014

# Perform filtering similar to the one in sQTLseekeR analysis of GEUVADIS
# Use sQTLseeker functions 

# Create dgeSQTL object 

##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

library(GenomicRanges)

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(sQTLseekeR)
sQTLseekeR.files <- list.files("/home/gosia/R/R_Multinomial_project/Analysis_GEUVADIS/sQTLseekeR-master/R/", full.names = TRUE)
for(i in sQTLseekeR.files)
source(i)

out.dir <- "DM_0.1.1_sQTL_analysis/Data/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


##########################################################################################
# prepare genotypes and transcript expression with filters as in sQTLseeker
##########################################################################################

data.dir <- "sQTLseekeR20_analysis/Data/"

## Input files: transcript expression, gene location and genotype information
trans.exp.f = paste0(data.dir, "trExpRPKM.tsv")
gene.bed.f = paste0(data.dir, "genes_noChr.bed")
genotype.f = paste0(data.dir, "snps_CEU_full.tsv")
genotype.indexed.f <- paste0(data.dir, "snps_CEU_full.tsv.bgz")


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

write.table(tre.df, paste0(out.dir, "tre.df.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(tre.df, file = paste0(out.dir, "tre.df.RData"))

######################### Filter SNPs

## Check if:
## - less than 3 missing genotype values
## - more than 5 samples per genotype group
## - more than 5 different splicing pts per genotype group
## Return: TRUE if snp passed, FALSE if not.
check.genotype <- function(geno.df, tre.df=NULL){
  apply(geno.df, 1, function(geno.snp){
    if(sum(as.numeric(geno.snp)==-1)>2){
      return("Missing genotype")
    } 
    geno.snp.t = table(geno.snp[geno.snp>-1])
    if(sum(geno.snp.t >= 5) < 2){
      return("One group of >5 samples")
    }
    ### Skipp this part of filtering 
    if(!is.null(tre.df)){
      nb.diff.pts = sapply(names(geno.snp.t)[geno.snp.t>1], function(geno.i){
        nbDiffPt(tre.df[,which(geno.snp==geno.i)])
      })
      if(sum(nb.diff.pts >= 5) < 2){
        return("One group of >5 different splicing")
      }
    }
    return("PASS")
  })
}


#### load genes

gene.bed = read.table(gene.bed.f, as.is=TRUE, sep="\t")
colnames(gene.bed) = c("chr","start","end","geneId")
genes.keep <- unique(tre.df$geneId)
gene.bed <- gene.bed[gene.bed$geneId %in% genes.keep,  ]
gene.bed <- gene.bed[gene.bed$chr %in% 1:22, ]

all(gene.bed$geneId %in% genes.keep)

length(genes.keep)
length(gene.bed$geneId)

genic.window=5e3


## for each chromosome find SNPs-genes matches and prepare genotypes and SNPs tables

genotypesList <- mclapply(1:22, function(chr){
  # chr = 22 
  cat(chr, fill = TRUE)

  
  gene.bed.tmp <- gene.bed[gene.bed$chr == chr, ]  
  gnrng <- GenomicRanges::GRanges(gene.bed.tmp$chr, IRanges::IRanges(gene.bed.tmp$start, gene.bed.tmp$end))  
  gnrng  <- GenomicRanges::resize(gnrng, GenomicRanges::width(gnrng) + 2 * genic.window, fix="center")
  
  genotype <- read.table(paste0(data.dir, "snps_CEU_chr", chr ,".tsv.sort.tsv"), header = FALSE, sep = "\t", as.is = TRUE)
  colnames(genotype) <- read.table(paste0(data.dir, "snps_CEU_head.tsv"), header = FALSE, sep = "\t", as.is = TRUE)
  
  gtrng <- GenomicRanges::GRanges(genotype$chr, IRanges::IRanges(genotype$start, genotype$end)) 
  
  all(colnames(genotype)[-c(1:4)] == colnames(tre.df)[-c(1,2)])
  
  ## Match SNPs and genes
  variantMatch <- GenomicRanges::findOverlaps(gnrng, gtrng, select = "all")
  length(variantMatch)
  
  tre.df.split <- split(tre.df[, -c(1,2)], f = tre.df$geneId)
  tre.df.split <- tre.df.split[gene.bed.tmp$geneId]
  
  tre.df.rel.split <- split(tre.df.rel[, -c(1,2)], f = tre.df.rel$geneId)
  tre.df.rel.split <- tre.df.rel.split[gene.bed.tmp$geneId]
  
  all(names(tre.df.split) == gene.bed.tmp$geneId)
  
  
  genotype.keep <- plyr::ldply(lapply(unique(queryHits(variantMatch)), function(gene.i){
    # gene.i <- 4
    tre.df.gene <- tre.df.split[[gene.i]]    
    tre.df.rel.gene <- tre.df.rel.split[[gene.i]]
    
    genotype.gene.full <- data.frame(geneId = gene.bed.tmp$geneId[gene.i], genotype[subjectHits(variantMatch)[queryHits(variantMatch) == gene.i],])
    
    genotype.gene.names <- genotype.gene.full[, c(1:5)]
    genotype.gene <- genotype.gene.full[, -c(1:5)]
    
    ## NA for samples with non expressed genes
    genotype.gene[, is.na(tre.df.gene[1,])] <- NA
  
    snps.to.keep = check.genotype(genotype.gene[,!is.na(tre.df.gene[1,])], tre.df.rel.gene[,!is.na(tre.df.gene[1,])])
#     snps.to.keep = check.genotype(genotype.gene[,!is.na(tre.df.gene[1,])])
    sum(snps.to.keep == "PASS")

    genotype.gene.full <- data.frame(genotype.gene.names[snps.to.keep == "PASS", ], genotype.gene[snps.to.keep == "PASS", ], stringsAsFactors = FALSE)
    
return(genotype.gene.full)
    
  }), identity)

  
  genotype.keep[genotype.keep == -1] <- NA
  
gc()

write.table(genotype.keep, paste0(out.dir, "genotypes_chr",chr,".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


  return(genotype.keep)
  
}, mc.cores = 22)


# genotype.keep <- genotypesList[[19]]
# for(i in 1:ncol(genotype.keep))
#   cat(sum(is.na(genotype.keep[,i])), " ")



genotypes <- do.call(rbind, genotypesList)

write.table(genotypes, paste0(out.dir, "genotypes.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

genotypes$geneId <- as.character(genotypes$geneId)
genotypes$snpId <- as.character(genotypes$snpId)

save(genotypes, file = paste0(out.dir, "genotypes.RData"))



length(unique(genotypes$geneId))
dim(genotypes)
dim(unique(genotypes[, c("geneId", "snpId")]))





# ###### check how many of them overpal with SNPs from sQTLseekeR 2.0
# 
# res.seeker <- read.table("sQTLseekeR20_analysis/Results/CEU_results_all.txt", header = TRUE, as.is = TRUE)
# 
# filt.seeker <- res.seeker[, c("geneId", "snpId")]
# 
# filt.dm <- genotypes[, c("geneId", "snpId")]
# 
# setdiff(unique(filt.dm$geneId), unique(filt.seeker$geneId))




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

write.table(tre.df, paste0(out.dir, "tre.df.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(tre.df, file = paste0(out.dir, "tre.df.RData"))


length(unique(tre.df$geneId))

######################### SNPs and genotype

### get the SNP - gene match from seeker results 

res.seeker.all <- read.table("sQTLseekeR20_analysis/Results/CEU_results_all.txt", header = TRUE, as.is = TRUE)

snp2gene <- unique(res.seeker.all[, c("geneId", "snpId")])


snps_CEU_full <- read.table(genotype.f, header = TRUE, as.is = TRUE)

save(snps_CEU_full, file = paste0(data.dir, "snps_CEU_full.RData"))

snps_CEU_filtered <- snps_CEU_full[snps_CEU_full$snpId %in% unique(res.seeker.all$snpId), ]

save(snps_CEU_filtered, file = paste0(data.dir, "snps_CEU_filtered.RData"))


genotypes <- merge(snp2gene, snps_CEU_filtered, by = "snpId", all.x = TRUE, sort = FALSE)

genotypes[genotypes == -1]  <- NA


#### Check if there are NAs
for(i in 1:ncol(genotypes))
  cat(sum(is.na(genotypes[,i])), " ")

for(i in 1:ncol(genotypes))
  cat(sum(genotypes[,i] == -1), " ")


write.table(genotypes, paste0(out.dir, "genotypes_basedOnSQTLseekeRresults.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


genotypes$geneId <- as.character(genotypes$geneId)
genotypes$snpId <- as.character(genotypes$snpId)

save(genotypes, file = paste0(out.dir, "genotypes_basedOnSQTLseekeRresults.RData"))




##########################################################################################
# Create dgeSQTL
##########################################################################################


# genotypes <- read.table(paste0(out.dir, "genotypes.txt"), header = TRUE, as.is = TRUE)
# tre.df <- read.table(paste0(out.dir, "tre.df.txt"), header = TRUE, as.is = TRUE)

# load(paste0(out.dir, "genotypes.RData"))

load(paste0(out.dir, "genotypes_basedOnSQTLseekeRresults.RData"))
load(paste0(out.dir, "tre.df.RData"))


library(edgeR)

dgeSQTL <- DGEList(counts = tre.df[, -c(1, 2)], genes = tre.df[, c(1, 2)])
rownames(dgeSQTL$counts) <- tre.df[,"trId"]

dgeSQTL$counts <- round(dgeSQTL$counts * 100) ### RPKM -> counts

dgeSQTL$samples <- colnames(tre.df[, -c(1, 2)])
colnames(dgeSQTL$genes) <- c("ete_id", "gene_id")

dgeSQTL$genotypes <- as.matrix(genotypes[, -c(1:5)])
dgeSQTL$SNPs <- genotypes[, c("snpId", "geneId", "chr", "start")]
colnames(dgeSQTL$SNPs)[1:2] <- c("SNP_id", "gene_id")

all(colnames(dgeSQTL$counts) == colnames(dgeSQTL$genotypes))


dgeSQTL$counts <- split(data.frame(dgeSQTL$counts), factor(dgeSQTL$genes$gene_id, levels = unique(dgeSQTL$genes$gene_id)))
dgeSQTL$counts <- lapply(dgeSQTL$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


save(dgeSQTL, file=paste0(out.dir, "dgeSQTL.RData"))





### subset the dgeSQTL object 
keep.genes <- intersect(unique(dgeSQTL$genes$gene_id), unique(dgeSQTL$SNPs$gene_id))[1:20]

dgeSQTL$counts <- dgeSQTL$counts[keep.genes]
dgeSQTL$genes <- dgeSQTL$genes[dgeSQTL$genes$gene_id %in% names(dgeSQTL$counts), ]
dgeSQTL$genotypes <- dgeSQTL$genotypes[dgeSQTL$SNPs$gene_id %in% names(dgeSQTL$counts), ]
dgeSQTL$SNPs <- dgeSQTL$SNPs[dgeSQTL$SNPs$gene_id %in% names(dgeSQTL$counts), ]


dgeSQTL$counts$ENSG00000159733.9


























