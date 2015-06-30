# BioC 3.0

# Created 27 Jan 2014

# Updated 27 Jan 2014
# Perform filtering that should be optimal for DM

# Updated 13 May 2015
# Create dgeSQTL object 

##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

library(GenomicRanges)

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)

out.dir <- "DM_0_1_3_Data/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


##########################################################################################
# paths to data used for filtering
##########################################################################################

data.dir <- "Data/"

## Input files: transcript expression, gene location and genotype information
trans.exp.f = paste0(data.dir, "trExpCount_CEU.tsv")
gene.bed.f = paste0(data.dir, "Annotation/genes_noChr.bed")

## Getting the IDs of samples in CEU population
groups.f <- paste0(data.dir, "Metadata/sample-groups.tsv")
groups = read.table(groups.f, header=TRUE, as.is=TRUE)
groups = subset(groups,group=="CEU")


##########################################################################################
# 1) Filter transcripts and genes for they min expression 
##########################################################################################

trans.exp <- read.table(trans.exp.f, header = TRUE, as.is = TRUE)

min.samps <- 10
# min.transcript.exp <- 1 # in cpm
min.transcript.prop <- 0.05 # in cpm
min.gene.exp <- 1 # in cpm


expr <-  trans.exp[, groups$sample]
colnames(expr) <- groups$sampleShort
tr.gene <- trans.exp[, c(1,2)]

seq.depth <- colSums(expr)

dge <- DGEList(counts = expr, genes = tr.gene)
expr.cpm <- cpm(dge)

rownames(expr.cpm) <- tr.gene$trId

expr.cpm.spl <- split(data.frame(expr.cpm), tr.gene$geneId) 

# which(names(expr.cpm.spl) == "ENSG00000168918.9")


trans.to.keep <- lapply(seq(length(expr.cpm.spl)), function(g){
  # g = 24525
  
  ### no genes with one transcript
  if(dim(expr.cpm.spl[[g]])[1] == 1)
    return(NULL)
  
  ### genes with min expression in all samples
  if(any(colSums(expr.cpm.spl[[g]]) < min.gene.exp))
    return(NULL)
  
  ### transcripts with min expression in samples
  # trans.to.keep = apply(expr.cpm.spl[[g]], 1, function(t) sum(t > min.transcript.exp) >= 10 )
  
  tot <- colSums(expr.cpm.spl[[g]])
  
  trans.to.keep = apply(expr.cpm.spl[[g]], 1, function(t){
    # t = expr.cpm.spl[[g]][8, ]
    
    r <- t / tot
    
    sum(r > min.transcript.prop) >= min.samps 
    
  } )
  
  ### no genes with one transcript
  if(sum(trans.to.keep) <= 1)
    return(NULL)
  
  return(names(which(trans.to.keep == TRUE)))
  
})

trans.to.keep <- unlist(trans.to.keep)

length(trans.to.keep)


expr <- expr[tr.gene$trId %in% trans.to.keep, ]
tr.gene <- tr.gene[tr.gene$trId %in% trans.to.keep, ]

length(unique(tr.gene$geneId))

# tr.gene[tr.gene$geneId == "ENSG00000168918.9", ]
# tr.gene[tr.gene$geneId == "ENSG00000054654.10", ]


tre.df <- data.frame(tr.gene, expr, stringsAsFactors = FALSE, row.names = NULL)



write.table(tre.df, paste0(out.dir, "tre.df.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(tre.df, file = paste0(out.dir, "tre.df.RData"))

load(paste0(out.dir, "tre.df.RData"))


pdf(paste0(out.dir, "/Hist_numberOfTranscripts.pdf"))
tt <- table(tre.df$geneId)
hist(tt, breaks = max(tt), col = "orangered", main = paste0(length(tt), " genes \n ", sum(tt) , " transcripts "), xlab = "number of transcripts per gene")
dev.off()


##########################################################################################
# 2) Filter SNPs
##########################################################################################


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
minSNP <- 10

tre.df.split <- split(tre.df[, -c(1,2)], f = tre.df$geneId)

## for each chromosome find SNPs-genes matches and prepare genotypes and SNPs tables

genotypesList <- mclapply(22:1, function(chr){
  # chr = 1
	 
  cat(chr, fill = TRUE)

  ### gene ranges
  gene.bed.tmp <- gene.bed[gene.bed$chr == chr, ]  
  gnrng <- GenomicRanges::GRanges(gene.bed.tmp$chr, IRanges::IRanges(gene.bed.tmp$start, gene.bed.tmp$end))  
  gnrng  <- GenomicRanges::resize(gnrng, GenomicRanges::width(gnrng) + 2 * genic.window, fix="center")
  
  genotype <- read.table(paste0(data.dir, "Genotypes/snps_CEU_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)
  rownames(genotype) <- genotype$snpId

  ### snps positions 
  gtrng <- GenomicRanges::GRanges(genotype$chr, IRanges::IRanges(genotype$start, genotype$end)) 
  
  all(colnames(genotype)[-c(1:4)] == colnames(tre.df)[-c(1,2)])
  
  ## Match genes and SNPs
  variantMatch <- GenomicRanges::findOverlaps(gnrng, gtrng, select = "all")
  length(variantMatch)
  
  variantMatch.names <- data.frame(geneId = gene.bed.tmp$geneId[queryHits(variantMatch)], snpId = genotype$snpId[subjectHits(variantMatch)], stringsAsFactors = FALSE)
  
  tre.df.split.tmp <- tre.df.split[gene.bed.tmp$geneId]
  
  
  genotype.keep <- plyr::ldply(lapply(unique(variantMatch.names$geneId), function(gene.i){
    # gene.i <- "ENSG00000177663.8"
    
    tre.df.gene <- tre.df.split[[gene.i]]    
 
    genotype.gene.tmp <-  genotype[variantMatch.names[variantMatch.names$geneId == gene.i, "snpId"], ]
    
    genotype.gene.names <- genotype.gene.tmp[, c(1:4)]
    genotype.gene <- genotype.gene.tmp[, -c(1:4)]
    
    ## NA for samples with non expressed genes and missing genotype
    genotype.gene[, is.na(tre.df.gene[1,])] <- NA
    genotype.gene[genotype.gene == -1] <- NA
    genotype.gene <- as.matrix(genotype.gene)
    
    ##### Keep genotypes with at least minSNP number of variants per group; in other put NAs

    genotype.to.keep <- apply(genotype.gene, 1 ,function(x){
      # x <- genotype.gene[6,]
      
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
      
    genotype.to.keep <- do.call(rbind, genotype.to.keep)
    
    genotype.gene.full <- data.frame(geneId = gene.i, genotype.gene.names[rownames(genotype.to.keep), ], genotype.to.keep, row.names = NULL, stringsAsFactors = FALSE)
   
return(genotype.gene.full)
    
  }), identity)


gc()

# write.table(genotype.keep, paste0(out.dir, "genotypes_chr",chr,".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


  return(genotype.keep)
  
}, mc.cores = 22)


# genotype.keep <- genotypesList[[19]]
# for(i in 1:ncol(genotype.keep))
#   cat(sum(is.na(genotype.keep[,i])), " ")


genotypes <- do.call(rbind, genotypesList)
genotypes$geneId <- as.character(genotypes$geneId)
genotypes$snpId <- as.character(genotypes$snpId)


# write.table(genotypes, paste0(out.dir, "genotypes.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
save(genotypes, file = paste0(out.dir, "genotypes.RData"))


load(paste0(out.dir, "/genotypes.RData"))

tt <- table(genotypes$geneId)
pdf(paste0(out.dir, "/Hist_numberOfSnps.pdf"))
hist(tt, breaks = 100, col = "chartreuse2", main = paste0(length(tt), " genes \n ", sum(tt) , " snps "), xlab = "number of snps per gene")
dev.off()





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


























