# BioC 3.0

# Created 28 May 2015

##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

library(GenomicRanges)

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)

out.dir <- "DM_0_1_5_Data/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


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

trans.exp <- read.table(trans.exp.f, header = TRUE, as.is = TRUE)

min_samps <- 5
min_transcript_prop <- 0.1 # in cpm
min_gene_exp <- 1 # in cpm


expr <-  trans.exp[, groups$sample]
colnames(expr) <- groups$sampleShort
rownames(expr) <- trans.exp$trId
expr_info <- trans.exp[, c(1,2)]

dge <- DGEList(counts = expr, genes = expr_info)
expr_cpm <- cpm(dge)

rownames(expr_cpm) <- expr_info$trId

expr_cpm_spl <- split(data.frame(expr_cpm), expr_info$geneId, drop = TRUE) 
expr_cpm_spl <- lapply(expr_cpm_spl, as.matrix)


trans2keep <- lapply(expr_cpm_spl, function(expr_gene){
# expr_gene <- expr_cpm_spl[[1230]]
  
  ### no genes with one transcript
  if(dim(expr_gene)[1] == 1)
    return(NULL)
  
  ### genes with min expression in all samples
  if(any(colSums(expr_gene) < min_gene_exp))
    return(NULL)
  
  prop <- prop.table(expr_gene, 2)
  
  trans2keep <- rowSums(prop > min_transcript_prop) > min_samps
  
  ### no genes with one transcript
  if(sum(trans2keep) <= 1)
    return(NULL)
  
  return(names(which(trans2keep == TRUE)))
  
})

trans2keep <- unlist(trans2keep, use.names = FALSE)


expr <- expr[expr_info$trId %in% trans2keep, ]
expr_info <- expr_info[expr_info$trId %in% trans2keep, ]


### split expression by gene
counts <- split(expr, expr_info$geneId, drop = TRUE) 
counts <- lapply(counts, as.matrix)

counts_info <- expr_info[, c("geneId", "trId")]
colnames(counts_info) <- c("gene_id", "ete_id")


save(counts, counts_info, file = paste0(out.dir, "/counts.RData"))


tt <- table(counts_info$gene_id)

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

gene_rngs <- gene_rngs[gene_rngs$name %in% genes_keep]
gene_rngs <- gene_rngs[seqnames(gene_rngs) %in% 1:22]


genic_window=5e3
minSNP <- 5


## for each chromosome find SNPs-genes matches and prepare genotypes and SNPs tables
genotypesList <- mclapply(2:20, function(chr){
  # chr = 1
  
  cat(chr, fill = TRUE)
  
  ### keep the data for one chromosome only
  gene_rngs_chr <- gene_rngs[seqnames(gene_rngs) == chr, ]  
  counts_chr <- counts[names(gene_rngs_chr)]
  
  
  ### read genotypes
  genotypes <- read.table(paste0(data.dir, "Genotypes/snps_CEU_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)
  rownames(genotypes) <- genotypes$snpId
  
  genotypes_info <- genotypes[, 1:4]
  genotypes <- genotypes[, -c(1:4)]
  
  ### snps ranges 
  snp_rngs_chr <- GenomicRanges::GRanges(genotypes_info$chr, IRanges::IRanges(genotypes_info$start, genotypes_info$end)) 
  
  all(colnames(genotypes) == colnames(counts[[1]]))
  
  ## Match genes and SNPs
  variantMatch <- GenomicRanges::findOverlaps(gene_rngs_chr, snp_rngs_chr, select = "all")
  length(variantMatch)
  
  q <- queryHits(variantMatch)
  s <- subjectHits(variantMatch)
  
  genotypes <- split(genotypes[s, ], names(gene_rngs_chr)[q], drop = TRUE)
  genotypes <- lapply(genotypes, as.matrix)
  
  
  geneList <- intersect(names(genotypes), names(counts_chr))
  
  
  genotypes2keep <- mapply(function(counts_gene, genotypes_gene, gene){
  
    # gene = "ENSG00000203836.5"; counts_gene = counts_chr[[gene]]; genotypes_gene = genotypes[[gene]]
    # print(gene)
    
    ## NA for samples with non expressed genes and missing genotype
    genotypes_gene[, is.na(counts_gene[1,])] <- NA
    genotypes_gene[genotypes_gene == -1] <- NA
    
    ##### Keep genotypes with at least minSNP number of variants per group; in other case replace them with NAs
    
    genotypes2keep_gene <- apply(genotypes_gene, 1 ,function(x){
      # x <- genotypes_gene[6,]
      
      t <- table(x)
      
      if( length(t)==1 )
        return(NULL)
      if( length(t)==2 ){
        if(any(t <= minSNP))
          return(NULL)
        return(x)
      }else{
        if(sum(t <= minSNP) >= 2)
          return(NULL)
        x[x == names(t[t <= minSNP])] <- NA
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
    
  }, counts_chr[geneList], genotypes[geneList], geneList, SIMPLIFY = FALSE)
  
  
  names(genotypes2keep) <- geneList
  
  genotypes2keep <- genotypes2keep[!sapply(genotypes2keep, is.null)]
  
  tt <- unlist(lapply(genotypes2keep, nrow))
  
  pdf(paste0(out.dir, "/Hist_numberOfSnps_chr",chr,".pdf"))
  hist(tt, breaks = 100, col = "chartreuse2", main = paste0(length(tt), " genes \n ", sum(tt) , " SNPs "), xlab = "Number of SNPs per gene")
  dev.off()
  
  geneList <- names(genotypes2keep)

  genotypes2keep_info <- plyr::ldply(lapply(geneList, function(gene){
    # gene = "ENSG00000131795.7"
    # print(gene)
    
#     if(is.null(genotypes2keep[[gene]]))
#       return(NULL)
    
    data.frame(gene_id = gene, snp_id = rownames(genotypes2keep[[gene]]), chr = chr, stringsAsFactors = FALSE)
    
  }), identity)
  

  save(genotypes2keep, genotypes2keep_info, file = paste0(out.dir, "genotypes_chr",chr ,".RData"))
  
  gc()
  return(NULL)
  
}, mc.cores = 5)






##########################################################################################
# Create dgeSQTL for chromosome 5
##########################################################################################


genotypes <- read.table(paste0(out.dir, "/genotypes_chr5.txt"), header = TRUE, as.is = TRUE)
# load(paste0(out.dir, "genotypes.RData"))
head(genotypes)


# tre.df <- read.table(paste0(out.dir, "/tre.df.txt"), header = TRUE, as.is = TRUE)
load(paste0(out.dir, "/tre.df.RData"))
head(tre.df)

### keep only the genes from chr5
tre.df <- tre.df[tre.df$geneId %in% genotypes$geneId, ]


library(edgeR)

dgeSQTL <- DGEList(counts = tre.df[, -c(1, 2)], genes = tre.df[, c(1, 2)])
rownames(dgeSQTL$counts) <- tre.df[,"trId"]

dgeSQTL$counts <- round(dgeSQTL$counts) ### RPKM -> counts round(dgeSQTL$counts * 100)

dgeSQTL$samples <- colnames(tre.df[, -c(1, 2)])
colnames(dgeSQTL$genes) <- c("ete_id", "gene_id")

dgeSQTL$genotypes <- as.matrix(genotypes[, -c(1:5)])
dgeSQTL$SNPs <- genotypes[, c("snpId", "geneId", "chr", "start")]
colnames(dgeSQTL$SNPs)[1:2] <- c("SNP_id", "gene_id")

all(colnames(dgeSQTL$counts) == colnames(dgeSQTL$genotypes))


dgeSQTL$counts <- split(data.frame(dgeSQTL$counts), factor(dgeSQTL$genes$gene_id, levels = unique(dgeSQTL$genes$gene_id)))
dgeSQTL$counts <- lapply(dgeSQTL$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


# save(dgeSQTL, file=paste0(out.dir, "/dgeSQTL.RData"))
save(dgeSQTL, file=paste0(out.dir, "/dgeSQTL_chr5.RData"))





# ### subset the dgeSQTL object 
# keep.genes <- intersect(unique(dgeSQTL$genes$gene_id), unique(dgeSQTL$SNPs$gene_id))[1:20]
# 
# dgeSQTL$counts <- dgeSQTL$counts[keep.genes]
# dgeSQTL$genes <- dgeSQTL$genes[dgeSQTL$genes$gene_id %in% names(dgeSQTL$counts), ]
# dgeSQTL$genotypes <- dgeSQTL$genotypes[dgeSQTL$SNPs$gene_id %in% names(dgeSQTL$counts), ]
# dgeSQTL$SNPs <- dgeSQTL$SNPs[dgeSQTL$SNPs$gene_id %in% names(dgeSQTL$counts), ]
# 
# 
# dgeSQTL$counts$ENSG00000159733.9


























