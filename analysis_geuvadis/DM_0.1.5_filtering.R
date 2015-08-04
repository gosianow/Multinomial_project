# BioC 3.0

# Created 28 May 2015


##############################################################################################################

setwd("/home/Shared/data/seq/geuvadis/")

library(GenomicRanges)

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)

out.dir <- "DM_0_1_5_Data/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

library(BiocParallel)

##########################################################################################
# paths to data used for filtering
##########################################################################################

data.dir <- "Data/"

## Input files: transcript expression, gene location and genotype information
trans.exp.f = paste0(data.dir, "Expression/trExpCount_CEU.tsv")
gene.bed.f = paste0(data.dir, "Annotation/genes_noChr.bed")

## Getting the IDs of samples in CEU population
groups.f <- paste0(data.dir, "Metadata/sample-groups.tsv")
groups = read.table(groups.f, header=TRUE, as.is=TRUE)
groups = subset(groups,group=="CEU")


##########################################################################################
# 1) Filter transcripts and genes for they min expression 
##########################################################################################

BPPARAM <- MulticoreParam(workers = 5)

trans.exp <- read.table(trans.exp.f, header = TRUE, as.is = TRUE)

min_transcript_prop <- 0.1
min_samps_transcript_prop <- 5
min_samps_gene_expr <- 70
min_gene_expr <- 1 # in cpm


expr <-  round(as.matrix(trans.exp[, groups$sample]))
colnames(expr) <- groups$sampleShort
rownames(expr) <- trans.exp$trId
expr_info <- trans.exp[, c(1,2)]

### calculate cpm
dge <- DGEList(counts = expr, genes = expr_info)
expr_cpm <- cpm(dge)
rownames(expr_cpm) <- expr_info$trId

expr_cpm_spl <- split.data.frame(expr_cpm, expr_info$geneId, drop = TRUE) 

expr_spl <- split.data.frame(expr, expr_info$geneId, drop = TRUE) 


counts <- bplapply(names(expr_cpm_spl), function(gene){
  # gene <- "ENSG00000164308.12"
  
  expr_cpm_gene <- expr_cpm_spl[[gene]]
  expr_gene <- expr_spl[[gene]]
  
  ### no genes with one transcript
  if(dim(expr_cpm_gene)[1] == 1)
    return(NULL)
  
  ### genes with min expression
  samps2keep <- colSums(expr_cpm_gene) > min_gene_expr
  
  if(sum(samps2keep) < min_samps_gene_expr)
    return(NULL)
  
  ### transcripts with min proportion
  prop <- prop.table(expr_gene[, samps2keep], 2)
  
  trans2keep <- rowSums(prop > min_transcript_prop) >= min_samps_transcript_prop
  
  ### no genes with one transcript
  if(sum(trans2keep) <= 1)
    return(NULL)
  
  expr <- expr_gene[trans2keep, ] 
  expr[, !samps2keep] <- NA
  
  return(expr)
  
}, BPPARAM = BPPARAM)

names(counts) <- names(expr_cpm_spl)


counts2keep <- !sapply(counts, is.null)

counts <- counts[counts2keep]

save(counts, file = paste0(out.dir, "/counts.RData"))


tt <- sapply(counts, nrow)

pdf(paste0(out.dir, "/Hist_numberOfTranscripts.pdf"))
hist(tt, breaks = max(tt), col = "orangered", main = paste0(length(tt), " genes \n ", sum(tt) , " transcripts "), xlab = "Number of transcripts per gene")
dev.off()



##########################################################################################
# 2) Filter SNPs
##########################################################################################


#### load genes
library(rtracklayer)
library(plyr)

load(paste0(out.dir, "/counts.RData"))

gene_rngs = import(gene.bed.f)
names(gene_rngs) <- gene_rngs$name

genes_keep <- names(counts)

genic_window = 5e3

gene_rngs <- gene_rngs[gene_rngs$name %in% genes_keep]
gene_rngs <- gene_rngs[seqnames(gene_rngs) %in% 1:22]
gene_rngs <- GenomicRanges::resize(gene_rngs, GenomicRanges::width(gene_rngs) + 2 * genic_window, fix="center")


minSNP <- 5

BPPARAM <- MulticoreParam(workers = 5)

## for each chromosome find SNPs-genes matches and prepare genotypes and SNPs tables
genotypesList <- bplapply(1:22, function(chr){
  # chr = 19
  print(chr)
  
  ### keep the data for one chromosome only
  gene_rngs_chr <- gene_rngs[seqnames(gene_rngs) == chr, ]  
  counts_chr <- counts[names(gene_rngs_chr)]
  
  ### read genotypes
  genotypes <- read.table(paste0(data.dir, "Genotypes/snps_CEU_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)
  rownames(genotypes) <- genotypes$snpId
  
  genotypes_info <- genotypes[, 1:4]
  genotypes <- as.matrix(genotypes[, -c(1:4)])
  
  ### snps ranges 
  snp_rngs_chr <- GenomicRanges::GRanges(genotypes_info$chr, IRanges::IRanges(genotypes_info$start, genotypes_info$end)) 
  
  ## Match genes and SNPs
  variantMatch <- GenomicRanges::findOverlaps(gene_rngs_chr, snp_rngs_chr, select = "all")
  length(variantMatch)
  
  q <- queryHits(variantMatch)
  s <- subjectHits(variantMatch)
  
  genotypes <- split.data.frame(genotypes[s, ], names(gene_rngs_chr)[q], drop = TRUE)
  
  geneList <- intersect(names(genotypes), names(counts_chr))
  
  
  genotypes2keep <- lapply(geneList, function(gene){
    
    # gene <- "ENSG00000105341.11"
    # print(gene)
    
    counts_gene = counts_chr[[gene]]
    genotypes_gene = genotypes[[gene]]
    
    
    ## NA for samples with non expressed genes and missing genotype
    genotypes_gene[, is.na(counts_gene[1,])] <- NA
    genotypes_gene[genotypes_gene == -1] <- NA
    
    ##### Keep genotypes with at least minSNP number of variants per group; in other case replace them with NAs
    genotypes2keep_gene <- apply(genotypes_gene, 1 ,function(x){
      # x <- genotypes_gene[6,]
      
      tt <- table(x)
      
      if( length(tt)==1 )
        return(NULL)
      if( length(tt)==2 ){
        if(any(tt <= minSNP))
          return(NULL)
        return(x)
      }else{
        if(sum(tt <= minSNP) >= 2)
          return(NULL)
        x[x == names(tt[tt <= minSNP])] <- NA
        return(x)
      }    
    })
    
    if(!is.null(genotypes2keep_gene)){
      if(is.list(genotypes2keep_gene))
        genotypes2keep_gene <- do.call(rbind, genotypes2keep_gene)
      else
        genotypes2keep_gene <- t(genotypes2keep_gene)
    }
    
    return(genotypes2keep_gene)
    
  })
  
  names(genotypes2keep) <- geneList
  genotypes2keep <- genotypes2keep[!sapply(genotypes2keep, is.null)]
  
  save(genotypes2keep, file = paste0(out.dir, "genotypes_chr",chr ,".RData"))
  
  return(NULL)
  
}, BPPARAM = BPPARAM)




### merge genotypes from all chromosomes

genotypes <- lapply(1:22, function(chr){
  
  load(paste0(out.dir, "genotypes_chr",chr ,".RData"))
  
  return(genotypes2keep)
  
})

genotypes <- unlist(genotypes, recursive = FALSE)



tt <- sapply(genotypes, nrow)

pdf(paste0(out.dir, "/Hist_numberOfSnps.pdf"))
hist(tt, breaks = 100, col = "chartreuse2", main = paste0(length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "Number of SNPs per gene")
dev.off()



save(genotypes, file = paste0(out.dir, "genotypes.RData"))


##########################################################################################
# Create dgeSQTL object
##########################################################################################

### load genotypes
load(paste0(out.dir, "genotypes.RData"))

### load counts
load(paste0(out.dir, "counts.RData"))


### keep genes that passed both filters

geneList <- intersect(names(genotypes), names(counts))
length(geneList)






library(edgeR)

dgeSQTL <- DGEList()

dgeSQTL$counts <- counts[geneList]
dgeSQTL$genotypes <- genotypes[geneList]
dgeSQTL$samples <- data.frame(sample_names = colnames(counts[[1]]), stringsAsFactors = FALSE)

save(dgeSQTL, file=paste0(out.dir, "/dgeSQTL.RData"))


tt <- sapply(dgeSQTL$genotypes, nrow)

pdf(paste0(out.dir, "/dgeSQTL_Hist_numberOfSnps.pdf"))
hist(tt, breaks = 100, col = "chartreuse2", main = paste0(length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "Number of SNPs per gene")
dev.off()


tt <- sapply(dgeSQTL$counts, nrow)

pdf(paste0(out.dir, "/dgeSQTL_Hist_numberOfTranscripts.pdf"))
hist(tt, breaks = max(tt), col = "orangered", main = paste0(length(tt), " genes \n ", sum(tt) , " transcripts "), xlab = "Number of transcripts per gene")
dev.off()
























