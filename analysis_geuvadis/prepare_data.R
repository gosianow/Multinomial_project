# BioC 3.0
# Created 9 Jan 2014

# Updated 9 Oct 2015


setwd("/home/Shared/data/seq/geuvadis/")


##################################################################################
### Prepare data 
##################################################################################

out.dir <- "data/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


##########################  gene.bed.f

library(GenomicRanges)
library(rtracklayer)

gtfFile <- "geuvadis_annotation/gencode.v12.annotation.gtf"

gtf0 <- import(gtfFile)

## keep protein coding genes
pcindex <- mcols(gtf0)$gene_type == "protein_coding" 
gtf <- gtf0[pcindex]

gindex <- mcols(gtf)$type == "gene" 
gtfg <- gtf[gindex]

gene.bed.f <- data.frame(chr = seqnames(gtfg), start =  start(gtfg), end = end(gtfg), geneId = mcols(gtfg)$gene_id)

write.table(gene.bed.f, paste0(out.dir,"genes.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



gene.bed <- read.table(paste0(out.dir,"genes.bed"))

gene.bed[,1] <- gsub(pattern = "chr", replacement = "", x = gene.bed[,1])

write.table(gene.bed, paste0(out.dir,"genes_noChr.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


##########################  groups

groups <- read.table("GEUVADIS_analysis_results/E-GEUV-1.sdrf.txt", header = T, sep="\t", as.is=TRUE)

groups <- groups[c("Assay.Name", "Characteristics.population.")]
groups <- unique(groups)

table(groups$Characteristics.population.)

colnames(groups) <- c("sample", "group")
groups$sampleShort <- substr(groups$sample, 1, 7)

write.table(groups, paste0(out.dir, "sample-groups.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) 




##########################  trans.exp.f

### RPKM
quantFile <- "GEUVADIS_analysis_results/GD660.TrQuantRPKM.txt"

trans.exp.f0 <- read.table(quantFile, header = T, sep="\t", as.is = TRUE)

trans.exp.f <- trans.exp.f0[, c("TargetID", "Gene_Symbol", groups$sample)]
colnames(trans.exp.f) <- c("trId", "geneId", groups$sample)

write.table(trans.exp.f, paste0(out.dir,"trExpRPKM.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

for(i in unique(groups$group))
  write.table(trans.exp.f[,c("trId", "geneId", groups$sample[groups$group == i])], paste0(out.dir,"trExpRPKM_",i,".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



### counts
quantFile <- "GEUVADIS_analysis_results/GD660.TrQuantCount.txt"

trans.exp.f0 <- read.table(quantFile, header = T, sep="\t", as.is = TRUE)

trans.exp.f <- trans.exp.f0[, c("TargetID", "Gene_Symbol", groups$sample)]
colnames(trans.exp.f) <- c("trId", "geneId", groups$sample)

write.table(trans.exp.f, paste0(out.dir,"trExpCount.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

for(i in unique(groups$group))
  write.table(trans.exp.f[,c("trId", "geneId", groups$sample[groups$group == i])], paste0(out.dir,"trExpCount_",i,".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##### check the relation between counts and RPKM


trans.exp.rpkm <- read.table(paste0(out.dir,"trExpRPKM_CEU.tsv"), header = T, sep="\t", as.is = TRUE)

trans.exp.counts <- read.table(paste0(out.dir,"trExpCount_CEU.tsv"), header = T, sep="\t", as.is = TRUE)

all(trans.exp.rpkm$trId == trans.exp.counts$trId)


pdf(paste0(out.dir,"Relation_between_RPKM&Count.pdf"))
plot(trans.exp.rpkm[,3], trans.exp.counts[,3])
dev.off()






##########################  genotype.f

library(VariantAnnotation)

tbindex <- list.files(path = "GEUVADIS_genotypes", pattern = "genotypes.vcf.bgz.tbi", full.names = TRUE, include.dirs = FALSE)
compressVcf <- substr(tbindex, 1, nchar(tbindex)-4)

library(limma)
chr <- gsub("chr", "", strsplit2(compressVcf, split=".", fixed=TRUE)[,2])


### extended gene ranges
gnrng <- gtfg 
range <- 5000
start(gnrng) <- start(gnrng) - range
end(gnrng) <- end(gnrng) + range
seqlevels(gnrng) <- gsub("chr", "", seqlevels(gnrng))




for(j in c("CEU", "FIN", "GBR", "TSI", "YRI")){
  # j = "CEU"

for(i in 1:length(compressVcf)){
  # i=11
  cat(j, chr[i], fill = TRUE)
  
  tab <- TabixFile(compressVcf[i], tbindex[i])
  gnrngTmp <- gnrng[seqnames(gnrng) == chr[i]]
  
  
  ## Explore the file header with scanVcfHeader
  hdr <- scanVcfHeader(tab)
  all(sampleShort %in% samples(hdr))
  

    
    ## read VCF file 
    param <- ScanVcfParam(which = gnrngTmp, samples = sampleShort[groups$group == j])
    
    vcf <- readVcf(tab, "hg19", param)
    #   vcf
    
    
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
    
    rowdata <- rowData(vcfbi)
    
    ## Convert genotype into number of mutated allels
    geno <- geno(vcfbi)$GT
    geno01 <- geno
    geno01[,] <- -1
    geno01[geno %in% c("0/0", "0|0")] <- 0 # AA = REF/REF
    geno01[geno %in% c("0/1", "0|1", "1/0", "1|0")] <- 1 # AB = REF/ALT
    geno01[geno %in% c("1/1", "1|1")] <- 2 # BB = ALT/ALT
    # geno01 should be integer, not character
    mode(geno01) <- "integer"
    
    genotype <- unique(data.frame(chr = seqnames(rowdata), start = start(rowdata), end = end(rowdata), snpId = rownames(geno01), geno01))
    
  ### sorting
  genotype <- genotype[order(genotype[,2]), ]
  
  
    write.table(genotype, file=paste0(out.dir, "snps_",j,"_chr" ,chr[i], ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
  gc()
  
  }
  

}

















































