# BioC 3.0
# Created 9 Jan 2014
# Modyfied 9 Jan 2014


# Run sQTLseekeR 2.0 analysis 

setwd("/home/Shared/data/seq/GEUVADIS/")


##################################################################################
### Prepare data for sQTLSeekeR
##################################################################################

out.dir <- "sQTLseekeR20_analysis/Data/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)



############## sort by SNP position and merge chromosome files

snpsFiles <- list.files(path = out.dir, pattern = "snps_CEU", full.names = TRUE, include.dirs = FALSE)

snpsFiles <- snpsFiles[grepl("chr",snpsFiles) & !grepl("sort",snpsFiles)]

x <- gsub("^.*chr","", snpsFiles)
xx <- gsub(".tsv$", "", x)

snpsOrder <- as.numeric(xx)

snpsFiles <- snpsFiles[order(snpsOrder)]

snpsFiles


for(i in 1:length(snpsFiles)){
  # i = 1
cat(i, fill = TRUE)
  #### Does not sort in the right way -->> use tab as separator
#   cmd <- paste0("(head -n 1 ", snpsFiles[i], " && tail -n +2 ", snpsFiles[i], " | sort -k 2 ) > ", snpsFiles[i], ".sort.tsv")
#   cmd <- paste0( "tail -n +2 ", snpsFiles[i], " | sort -k 2  > ", snpsFiles[i], ".sort.tsv")  
#   system(cmd)
  
  d <- read.table(snpsFiles[i], header = TRUE, as.is = TRUE)
	
  o <- order(d[,2])
	
	print(table(o[-1] - o[-length(o)]))
	
  do <- d[o, ]
  
  # write.table(do, file=paste0(snpsFiles[i] ,".sort.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  gc()
  
}


cmd <- paste0("cat ", paste0(paste0(snpsFiles, ".sort.tsv"), collapse = " "), " > sQTLseekeR20_analysis/Data/snps_CEU.tsv")

system(cmd)

cmd <- paste0("head -n 1 ", snpsFiles[1], " > sQTLseekeR20_analysis/Data/snps_CEU_head.tsv")

system(cmd)

cmd <- paste0("cat sQTLseekeR20_analysis/Data/snps_CEU_head.tsv sQTLseekeR20_analysis/Data/snps_CEU.tsv > sQTLseekeR20_analysis/Data/snps_CEU_full.tsv")

system(cmd)


samples.head <- read.table("sQTLseekeR20_analysis/Data/snps_CEU_head.tsv")




##################################################################################
### Run sQTLseekeR 2.0 analysis 
##################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

library(sQTLseekeR)

out.dir <- "sQTLseekeR20_analysis/Results/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

data.dir <- "sQTLseekeR20_analysis/Data/"

## Input files: transcript expression, gene location and genotype information
trans.exp.f = paste0(data.dir, "trExpRPKM.tsv")
gene.bed.f = paste0(data.dir, "genes_noChr.bed")
genotype.f = paste0(data.dir, "snps_CEU_full.tsv")

## Getting the IDs of samples in CEU population
groups.f <- paste0(data.dir, "sample-groups.tsv")
groups = read.table(groups.f, header=TRUE, as.is=TRUE)
groupsCEU = subset(groups,group=="CEU")

### check if the samples order is the same
# all(groupsCEU$sampleShort == samples.head[-c(1:4)])


## 1) Index the genotype file (if not done externally before)
# genotype.indexed.f = index.genotype(genotype.f)
genotype.indexed.f <- paste0(data.dir, "snps_CEU_full.tsv.bgz")


## 2) Prepare transcript expression
# te.df.all = read.table(trans.exp.f, as.is=TRUE, header=TRUE, sep="\t")
# 
# te.df = te.df.all[,c("trId", "geneId", groupsCEU$sample)]
# colnames(te.df) <- c("trId", "geneId", groupsCEU$sampleShort)
# 
# tre.df = prepare.trans.exp(te.df)
# 
# write.table(tre.df, paste0(data.dir, "trExpRPKM_CEU_Seeker.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

tre.df <- tre.df.org <- read.table(paste0(data.dir, "trExpRPKM_CEU_Seeker.tsv"), header = TRUE, as.is=TRUE)


## 3) Test gene/SNP associations
gene.bed = read.table(gene.bed.f, as.is=TRUE, sep="\t")
colnames(gene.bed) = c("chr","start","end","geneId")



library(parallel)

genes.unq <- unique(tre.df.org$geneId)
dir.create(paste0(out.dir, "ByGene_CEU/"))


res.df.list <- mclapply(seq(length(genes.unq)), function(g){
  cat(g, fill = TRUE)
  res.df = sQTLseekeR::sqtl.seeker(tre.df.org[tre.df.org$geneId == genes.unq[g],  ], genotype.indexed.f, gene.bed)
  
  if(!is.null(res.df))
  write.table(res.df, paste0(out.dir, "ByGene_CEU/results_", g, genes.unq[g], ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  return(NULL)
  
}, mc.cores = 20)


##### read in the results 
res.df.files <- list.files(paste0(out.dir, "ByGene_CEU/"), full.names=TRUE)

res.df.list <- mclapply(seq(length(res.df.files)), function(g){

  res <- read.table(res.df.files[g], header = TRUE, as.is = TRUE)
  
}, mc.cores = 10)

res.df <- do.call(rbind, res.df.list)


write.table(res.df, paste0(out.dir, "CEU_results_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

dim(res.df)
dim(unique(res.df[,c("geneId", "snpId")]))

## 4) Get significant sQTLs
sqtls.df = sqtls(res.df, FDR=.05, out.pdf=paste0(out.dir, "CEU_results_FDR05.pdf"))


write.table(sqtls.df, paste0(out.dir, "CEU_results_FDR05.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)























































