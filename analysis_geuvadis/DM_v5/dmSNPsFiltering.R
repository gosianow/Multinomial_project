################################################################################

# Created 14 Oct 2014
# Updated 14 Oct 2014
# BioC 14

# Update of _back1 - wrapped into nice functions 

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
pcindex <- mcols(gtf0)$gene_type == "protein_coding" 
gtf0 <- gtf0[pcindex]

## keep genes with at least two isoforms
tindex <- mcols(gtf0)$type == "transcript" 
gtft <- gtf0[tindex]

gindex <- mcols(gtf0)$type == "gene" 
gtfg <- gtf0[gindex]

## 20,110
length(unique(mcols(gtft)$gene_id))


length(unique(mcols(gtft)$transcript_id))
length(mcols(gtft)$transcript_id)

oneIso <- (duplicated(mcols(gtft)$gene_id, fromLast = T) + duplicated(mcols(gtft)$gene_id, fromLast = F)) == 0

oindex <- !oneIso
gtff <- gtft[oindex]

## 16,581
length(unique(mcols(gtff)$gene_id))


### genes with at least two transcripts
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
sampsShort <- substr(samps, 1, 7)


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

length(goodGenes2) # here 9'619/9'249 genes, in sQTLseekeR the have 10'012 genes left


# gtfGood <- gtf0[mcols(gtf0)$gene_id %in% goodGenes2]
# seqlevels(gtfGood) <- gsub("chr", "", seqlevels(gtfGood))


################################################################################
# SNPs from genotypes
################################################################################

########### Convert to Big zip and create Tabix index for all vcf files
library(R.utils)
library(Rsamtools)
# 
# vcfFiles <- list.files(path = "genotypes", pattern = "genotypes.vcf.gz", full.names = TRUE, include.dirs = FALSE)
# 
# mclapply(11:length(vcfFiles), function(i){
#   # i=6
#   vcfFile <- vcfFiles[i]
#   print(vcfFile)
#    gunzip(vcfFile, remove = FALSE)
# vcfFileTmp <- substr( vcfFile , start = 1 , stop = nchar(vcfFile) - 3)
#   compressVcf <- bgzip(vcfFileTmp, overwrite = TRUE)
#   tbindex <- indexTabix(compressVcf, "vcf")
# #   system(cat("rm", vcfFileTmp, "\n"))
# 
#   return(vcfFile)
#   
# }, mc.cores=5)
  

########### Prepare filtered genotypes files 
library(VariantAnnotation)

outDir <- "DM_genotypes"
dir.create(outDir)

tbindex <- list.files(path = "genotypes", pattern = "genotypes.vcf.bgz.tbi", full.names = TRUE, include.dirs = FALSE)
compressVcf <- substr(tbindex, 1, nchar(tbindex)-4)

library(limma)
chr <- gsub("chr", "", strsplit2(compressVcf, split=".", fixed=TRUE)[,2])


## prepare ranges +-5Kb for goodGenes2  
gindex <- mcols(gtfg)$gene_id %in% goodGenes2 
gtfg <- gtfg[gindex]
gnrng <- gtfg ## good genes extended ranges!!!
range <- 5000
start(gnrng) <- start(gnrng) - range
end(gnrng) <- end(gnrng) + range
seqlevels(gnrng) <- gsub("chr", "", seqlevels(gnrng))


genoABcleanList <- list()
vcfCleanList <- list()


for(i in 1:length(compressVcf)){
  # i=5
  
  tab <- TabixFile(compressVcf[i], tbindex[i])
  gnrngTmp <- gnrng[seqnames(gnrng) == chr[i]]

  ## Explore the file header with scanVcfHeader
  hdr <- scanVcfHeader(tab)
  all(sampsShort %in% samples(hdr))
  
  
  ## Load VCF 
  param <- ScanVcfParam(which = gnrngTmp, samples = sampsShort) # does not work
  #param <- ScanVcfParam(which = gnrmgTmp)
  
  ## Extract the goodGenes2 ranges from the VCF file 
  vcf <- readVcf(tab, "hg19", param)
  vcf
  
  
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
  
  
  
  ## Keep the SNPs that have at least 2 different allels present at least 5 times in population  
  minSNP <- 5
    
  genoABcleanTmp <- lapply(seq_len(nrow(genoAB)), function(i){ # maybe use switch
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
  
  genoABclean <- do.call(rbind, genoABcleanTmp) ## FINAL
  genoABcleanList[[chr[i]]] <- genoABclean
  
  write.table(genoABclean, file=paste0(outDir, "/genoABclean_chr", chr[i], ".txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  goodSNPs <- rownames(genoABclean) ## FINAL good SNPs
  
	
	### keep vcf to access positions of SNPs
  vcfClean <- vcf[goodSNPs, ]
  vcfCleanList[[chr[i]]] <- vcfClean 
  
  writeVcf(vcfClean, paste0(outDir, "/vcfClean_chr", chr[i], ".vcf"), index = TRUE)
  

}




###############################
# prepare genotypes & SNPs (gene_id match SNP_id) for DGEList
###############################

SNPsList <- list()
genotypesList <- list()

for(i in 1:length(vcfCleanList)){
  # i = 11
  print(chr[i])
  gnrngTmp <- gnrng[seqnames(gnrng) == chr[i]]
  vcfClean <- vcfCleanList[[chr[i]]]
    
  ## Match SNPs and genes
  variantMatch <- findOverlaps(gnrngTmp, rowData(vcfClean), select = "all")
  length(variantMatch)
  
  ## check if values are monotonic
#   all(queryHits(variantMatch) == cummax(queryHits(variantMatch)))
#   all(subjectHits(variantMatch) == cummax(subjectHits(variantMatch)))
#   
  
  ### SNPs  
  SNPsList[[chr[i]]] <- data.frame(gene_id = mcols(gnrngTmp)$gene_id[queryHits(variantMatch)], SNP_id = names(rowData(vcfClean))[subjectHits(variantMatch)], stringsAsFactors = FALSE)
   
  ### genotypes  
  genotypesList[[chr[i]]] <- genoABcleanList[[chr[i]]][SNPsList[[chr[i]]]$SNP_id, ]  
  
}

SNPs <- do.call(rbind, SNPsList)
genotypes <- do.call(rbind, genotypesList)




###############################
# prepare counts & genes for DGEList
###############################



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


save(dgeSQTL, file="DM_genotypes/dgeSQTL.RData")






################################ dgeSQTL object with genotypes as list // DM_v5 functions are not adjusted for it yet


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
dgeSQTL$genotypes <- split(data.frame(dgeSQTL$genotypes), factor(dgeSQTL$SNPs$gene_id, levels = unique(dgeSQTL$SNPs$gene_id))) ## Taaaakes a lot of time!
dgeSQTL$genotypes <- lapply(dgeSQTL$genotypes, as.matrix)  ## !!! have to conver into matrix, othewise ERROR



save(dgeSQTL, file="DM_genotypes/dgeSQTL.RData")



















