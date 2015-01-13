################################################################################

# Created 02 Oct 2014
# Updated 14 Oct 2014
# BioC 14

# Keep the relavant SNPs from GEUVADIS data:
# a) Genes: protein coding, with at least two alternative transcript isoforms
# b) Samples:
# c) SNPs: +- 5Kb around genes 

################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")



################################################################################
# Genes from Genecode v12
################################################################################



library(GenomicRanges)
library(rtracklayer)

gtfFile <- "/home/Shared_penticton/data/annotation/Human/gencode/gencode.v12.annotation.gtf"

gtf0 <- import(gtfFile)

## keep protein coding genes
idx <- mcols(gtf0)$gene_type == "protein_coding" 
gtf0 <- gtf0[idx]

## keep genes with at least two isoforms
idx <- mcols(gtf0)$type == "transcript" 
gtft <- gtf0[idx]

idx <- mcols(gtf0)$type == "gene" 
gtfg <- gtf0[idx]

## 20,110
length(unique(mcols(gtft)$gene_id))


length(unique(mcols(gtft)$transcript_id))
length(mcols(gtft)$transcript_id)

oneIso <- (duplicated(mcols(gtft)$gene_id, fromLast = T) + duplicated(mcols(gtft)$gene_id, fromLast = F)) == 0

idx <- !oneIso
gtff <- gtft[idx]

## 16,581
length(unique(mcols(gtff)$gene_id))


goodGenes <- unique(mcols(gtff)$gene_id)


################################################################################
# Samples and Genes from expression quantification data
################################################################################

metadata <- read.table("analysis_results/E-GEUV-1.sdrf.txt", header = T, sep="\t", stringsAsFactors = FALSE)
metadata <- metadata[c("Assay.Name", "Characteristics.population.", "Factor.Value.laboratory.")]
metadata <- unique(metadata)

quantFile <- "analysis_results/GD660.TrQuantCount.txt"

qNames <- names(read.table(quantFile, header = T, nrows = 1, stringsAsFactors = FALSE))

## select only CEU columns 
mycols <- rep("NULL", length(qNames))
samps <- metadata[metadata$Characteristics.population. == "CEU", "Assay.Name"]
filter <- qNames %in% c("TargetID", "Gene_Symbol", "Chr", "Coord", samps )
mycols[filter] <- NA

quant <- read.table(quantFile, colClasses=mycols, header = T, stringsAsFactors = FALSE)

qGenes <- unique(quant$Gene_Symbol)

## keep quantification for genes with at least two transcripts
qKeep <- quant$Gene_Symbol %in% goodGenes
quant <- quant[qKeep, ]

## keep quantification for expressed transcripts
library(edgeR)
dge <- DGEList(counts=round(quant[,-c(1:4)]), genes=quant[,1:4])
qKeep2 <- rowSums(cpm(dge) > 1) >= ncol(dge$counts)/4
dge <- dge[qKeep2, ]

table(table(dge$genes$Gene_Symbol[!qKeep2]))

## keep genes with at least two isoforms 2
oneIso <- (duplicated(dge$genes$Gene_Symbol, fromLast = T) + duplicated(dge$genes$Gene_Symbol, fromLast = F)) == 0
dge <- dge[!oneIso, ]


dge ## FINAL


goodGenes2 <- unique(dge$genes$Gene_Symbol) ## FINAL good genes

length(goodGenes2) # in sQTLseekeR the have 10'012 genes left

# gtfGood <- gtf0[mcols(gtf0)$gene_id %in% goodGenes2]
# seqlevels(gtfGood) <- gsub("chr", "", seqlevels(gtfGood))


################################################################################
# SNPs from genotypes
################################################################################

library(VariantAnnotation)

## VCF for chr 19
vcfFile <- "genotypes/GEUVADIS.chr19.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"

## Create Tabix index
library(R.utils)

# gunzip(vcfFile, remove = FALSE)
# compressVcf <- bgzip(substr( vcfFile , start = 1 , stop = nchar(vcfFile) - 3), overwrite = TRUE)
# idx <- indexTabix(compressVcf, "vcf")

compressVcf <-  "genotypes/GEUVADIS.chr19.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.bgz"
idx <-  "genotypes/GEUVADIS.chr19.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.bgz.tbi" 


tab <- TabixFile(compressVcf, idx)


## Explore the file header with scanVcfHeader
hdr <- scanVcfHeader(tab)
hdr
info(hdr) 
geno(hdr) 

## prepare ranges +-5Kb for goodGenes2

idx <- mcols(gtfg)$gene_id %in% goodGenes2 
gtfg <- gtfg[idx]

gnrng <- gtfg ## good genes extended ranges!!!
start(gnrng) <- start(gnrng) - 5000
end(gnrng) <- end(gnrng) + 5000

seqlevels(gnrng)
seqlevels(gnrng) <- gsub("chr", "", seqlevels(gnrng))

gnrng <- gnrng[seqnames(gnrng) == "19"]


## Load VCF

sampsShort <- substr(samps, 1, 7)
all(sampsShort %in% samples(hdr))

param <- ScanVcfParam(which = gnrng, samples = sampsShort) # does not work
#param <- ScanVcfParam(which = gnrng)

## Extract the goodGenes2 ranges from the VCF file 
vcf <- readVcf(tab, "hg19", param)
vcf

head(fixed(vcf))
dim(fixed(vcf))
geno(hdr)
geno(vcf)
rowData(vcf)

head(geno(vcf)$GT)


# GT <- readGT(tab, nucleotides=FALSE, param)


## Keep only bi-allelic SNPs

# width of ref seq
rw <- width(ref(vcf))
# width of first alt seq
aw <- unlist(lapply(alt(vcf), function(x) {width(x[1])}))
# number of alternate genotypes
nalt <- elementLengths(alt(vcf))
# select only bi-allelic SNPs (monomorphic OK, so aw can be 0 or 1)
snp <- rw == 1 & aw <= 1 & nalt == 1
# subset vcf

vcfbi <- vcf[snp,]

## Convert genotype into number of mutated allels
geno <- geno(vcfbi)$GT
genoAB <- geno
genoAB[,] <- NA
genoAB[geno %in% c("0/0", "0|0")] <- 0 # AA = REF/REF
genoAB[geno %in% c("0/1", "0|1", "1/0", "1|0")] <- 1 # AB = REF/ALT
genoAB[geno %in% c("1/1", "1|1")] <- 2 # BB = ALT/ALT
# genoAB should be integer, not character
mode(genoAB) <- "integer"



## Keep the SNPs that have 2 different allels present in population  
minSNP <- 5


# snp2 <- apply(genoAB, MARGIN = 1, function(x){ # maybe use switch
#   ## assign NA to the allels with not enought counts that have TRUE
#   t <- table(x)
#   if( length(t)==1 )
#     return(FALSE)
#   if( length(t)==2 ){
#     if(any(t <= minSNP))
#       return(FALSE)
#     return(TRUE)
#   }else{
#     if(sum(t <= minSNP) == 2)
#       return(FALSE)
#     x[x == names(t[t <= minSNP])] <- NA
#     return(TRUE)
#   }    
# }) ## can not modify genoAB



genoABclean <- lapply(seq_len(nrow(genoAB)), function(i){ # maybe use switch
  ## assign NA to the allels with not enought counts that have TRUE
  x <- genoAB[i, , drop=FALSE]
  t <- table(x)
  if( length(t)==1 )
    return(NULL)
  if( length(t)==2 ){
    if(any(t <= minSNP))
      return(NULL)
    return(x)
  }else{
    if(sum(t <= minSNP) == 2)
      return(NULL)
    x[x == names(t[t <= minSNP])] <- NA
    return(x)
  }    
})

genoABclean <- do.call(rbind, genoABclean) ## FINAL

write.table(genoABclean, file="DM_sQTL_test/genoABclean.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

goodSNPs <- rownames(genoABclean) ## FINAL good SNPs

vcfClean <- vcf[goodSNPs, ]

writeVcf(vcfClean, "DM_sQTL_test/vcfClean.vcf", index = TRUE)


################################################################################
# prepare genotypes & SNPs objects for DGEList
################################################################################

## Check if variants match more then one gene
variantMatch <- findOverlaps(gnrng, rowData(vcfClean), select = "all")
length(variantMatch)

## check if values are monotonic
all(queryHits(variantMatch) == cummax(queryHits(variantMatch)))
all(subjectHits(variantMatch) == cummax(subjectHits(variantMatch)))

## Did any variants match more than one gene?
table(table(subjectHits(variantMatch)))

oneGeneSNPs <- duplicated(subjectHits(variantMatch), fromLast = T) + duplicated(subjectHits(variantMatch), fromLast = F) == 0
sum(oneGeneSNPs)

## How many variants the genes have?
table(table(queryHits(variantMatch)))




# save.image("DM_sQTL_test/workspace.RData")



### SNPs

SNPs <- data.frame(gene_id = mcols(gnrng)$gene_id[queryHits(variantMatch)], SNP_id = names(rowData(vcfClean))[subjectHits(variantMatch)], stringsAsFactors = FALSE)

geneOrder <- unique(SNPs$gene_id)

### genotypes

genotypes <- genoABclean[SNPs$SNP_id, ]


### genes 

genes <- data.frame(gene_id = dge$genes$Gene_Symbol, ete_id = dge$genes$TargetID, Chr = dge$genes$Chr, stringsAsFactors = FALSE)



### counts 

colnames(dge$counts) <- substr(colnames(dge$counts), 1, 7)
sampleOrder <- colnames(genotypes)
counts <- dge$counts[, sampleOrder]
rownames(counts) <- genes$ete_id
 


### samples 


samples <- data.frame(name = colnames(genotypes), stringsAsFactors = FALSE)
  

################################ dgeSQTL object
  
library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v5.R"))




dgeSQTL <- DGEList(counts = counts, genes = genes)

dgeSQTL$samples <- samples
dgeSQTL$genotypes <- genotypes
dgeSQTL$SNPs <- SNPs
  

dgeSQTL$counts <- split(data.frame(dgeSQTL$counts), factor(dgeSQTL$genes$gene_id, levels = unique(dgeSQTL$genes$gene_id)))
dgeSQTL$counts <- lapply(dgeSQTL$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


save(dgeSQTL, file="DM_sQTL_test/dgeSQTL.RData")
save(dgeSQTL, file="/home/Shared/tmp/gosia/Example_DM/dgeSQTL.RData")



#### split genotypes by gene and order acordingly 

dgeSQTL <- DGEList(counts = counts, genes = genes)

dgeSQTL$samples <- samples
dgeSQTL$genotypes <- genotypes
dgeSQTL$SNPs <- SNPs


### Keep data for the same genes
geneInters <- intersect(dgeSQTL$genes$gene_id, dgeSQTL$SNPs$gene_id)

dgeSQTL$counts <- dgeSQTL$counts[dgeSQTL$genes$gene_id %in% geneInters,]
dgeSQTL$genes <- dgeSQTL$genes[dgeSQTL$genes$gene_id %in% geneInters,]


dgeSQTL$genotypes <- dgeSQTL$genotypes[dgeSQTL$SNPs$gene_id %in% geneInters,]
dgeSQTL$SNPs <- dgeSQTL$SNPs[dgeSQTL$SNPs$gene_id %in% geneInters,]


### Split counts
dgeSQTL$counts <- split(data.frame(dgeSQTL$counts), factor(dgeSQTL$genes$gene_id, levels = unique(dgeSQTL$genes$gene_id)))
dgeSQTL$counts <- lapply(dgeSQTL$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


### Order genotypes and SNPs as counts
dgeSQTL$genotypes <- dgeSQTL$genotypes[order(match(dgeSQTL$SNPs$gene_id,names(dgeSQTL$counts))), ]
dgeSQTL$SNPs <- dgeSQTL$SNPs[order(match(dgeSQTL$SNPs$gene_id,names(dgeSQTL$counts))), ]

### Split genotypes
dgeSQTL$genotypes <- split(data.frame(dgeSQTL$genotypes), factor(dgeSQTL$SNPs$gene_id, levels = unique(dgeSQTL$SNPs$gene_id)))
dgeSQTL$genotypes <- lapply(dgeSQTL$genotypes, as.matrix)  ## !!! have to conver into matrix, othewise ERROR





save(dgeSQTL, file="DM_sQTL_test/dgeSQTL.RData")
save(dgeSQTL, file="/home/Shared/tmp/gosia/Example_DM/dgeSQTL.RData")


x <- dgeSQTL$SNPs
y <- names(dgeSQTL$counts)


z <- x[order(match(x$gene_id,y)), ]



x <- c(2, 2, 3, 4, 1, 4, 4, 3, 3)
y <- c(4, 2, 1, 3)


z <- y[sort(order(y)[x])]



x[order(match(x,y))]



















