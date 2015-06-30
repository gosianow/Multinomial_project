######################################################
# BioC 2.14
# Created 04 Nov 2014 

# Run DEXSeq_1.10.8

# Updated 11 June 2015
# Add Full_tumorVSnormal comparison
# Updated 22 June 2015
# Add kallisto counts

#######################################################



setwd("/home/Shared/data/seq/Kim_adenocarcinoma/")


library("DEXSeq")



######################################################################################################
# metadata
######################################################################################################

metadata <- read.table("3_metadata/Malgorzata Nowicka2014-11-04GSE37764.csv", stringsAsFactors=F, sep=",", header=T) 
metadata <- metadata[metadata$X == "RNA-seq",]

metadata$sampleName <- metadata$ids
metadata$condition <- metadata$Tissue.Type

metadata <- metadata[order(metadata$condition), ]

metadata


######################################################################################################
# htseq counts
######################################################################################################

counts_out <- "2_counts/htseq/"
dir.create(counts_out)


gtf= "/home/Shared_penticton/data/annotation/Human/Ensembl_GRCh37.71/gtf/Homo_sapiens.GRCh37.71.gtf"

DEXSeq.gff = "/home/Shared/data/annotation/Human/Ensembl_GRCh37.71/DEXSeq_1.10.8_gff/Homo_sapiens.GRCh37.71.DEXSeq.flattened.rNO.gff"

# python scripts
pkgDir <- system.file(package="DEXSeq")
pythDir <- file.path(pkgDir, "python_scripts")
list.files(pythDir)


### crerate gff file # disable aggregation with the option “-r no”
system(paste0("python ", pythDir, "/dexseq_prepare_annotation.py --help "))

python.cmd1 <- paste0("python ", pythDir, "/dexseq_prepare_annotation.py -r no ", gtf, " ", DEXSeq.gff)
cat(python.cmd1)
system(python.cmd1)



### add chr to gff
## awk '{print "chr"$0}' Homo_sapiens.GRCh37.71.DEXSeq.flattened.rNO.gff > Homo_sapiens.GRCh37.71.DEXSeq.flattened.rNO.chr.gff

DEXSeq.gff <- "/home/Shared/data/annotation/Human/Ensembl_GRCh37.71/DEXSeq_1.10.8_gff/Homo_sapiens.GRCh37.71.DEXSeq.flattened.rNO.chr.gff"


### exon counts # by name
system(paste0("python ", pythDir, "/dexseq_count.py --help "))

python.cmd2 <- with(metadata, paste0("python ", pythDir, "/dexseq_count.py -p yes -s no -f bam -r pos ", DEXSeq.gff, " bam_insilicodb/", ids, "_s.bam ", counts_out, ids , ".counts \n"))
cat(python.cmd2)

 for(i in 1:nrow(metadata)){
   # i = 1 
   system(python.cmd2[i])
   
 }


######################################################################################################
# DEXSeq models
######################################################################################################

library(BiocParallel)

BPPARAM = MulticoreParam( workers = 10)
allDEXSeq <- list()


DEXSeq.gff <- "/home/Shared/data/annotation/Human/Ensembl_GRCh37.71/DEXSeq_1.10.8_gff/Homo_sapiens.GRCh37.71.DEXSeq.flattened.rNO.chr.gff"


method_out <- "4_results/DEXSeq_1.10.8/"

count.methodList <- c("htseq", "kallisto")
count.method <- count.methodList[2]

counts_out <- paste0("2_counts/", count.method, "/")






# Null_normal1
model <- "Null_normal1"
dir.create(paste0(method_out, "/diff_out_", model))

metadata <-  metadataOrg[metadataOrg$condition == "normal", ]
countFiles <- paste0(counts_out, metadata$sampleName, ".counts")

sampleTable = data.frame(row.names = metadata$sampleName, condition = c(rep("C1", 3), rep("C2", 3)))

dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design = ~ sample + exon + condition:exon , flattenedfile = NULL ) 

exptData(dxd)

formulaFullModel = ~ sample + exon + condition:exon
formulaReducedModel =  ~ sample + exon 

allDEXSeq[[model]] <- DEXSeq(dxd, fullModel = formulaFullModel, reducedModel = formulaReducedModel, BPPARAM=BPPARAM, fitExpToVar="condition")




# Null_tumor1
model <- "Null_tumor1"
dir.create(paste0(method_out, "/diff_out_", model))

metadata <-  metadataOrg[metadataOrg$condition == "tumor", ]
countFiles <- paste0(counts_out, metadata$sampleName, ".counts")

sampleTable = data.frame(row.names = metadata$sampleName, condition = c(rep("C1", 3), rep("C2", 3)))

dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design = ~ sample + exon + condition:exon , flattenedfile = NULL ) 

exptData(dxd)

formulaFullModel = ~ sample + exon + condition:exon
formulaReducedModel =  ~ sample + exon 

allDEXSeq[[model]] <- DEXSeq(dxd, fullModel = formulaFullModel, reducedModel = formulaReducedModel, BPPARAM=BPPARAM, fitExpToVar="condition")




# Full_tumorVSnormal
model <- "Full_tumorVSnormal"
dir.create(paste0(method_out, "/diff_out_", model), showWarnings = FALSE)

countFiles <- paste0(counts_out, metadata$sampleName, ".counts")

sampleTable = data.frame(row.names = metadata$sampleName, condition = metadata$condition)

dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData = sampleTable, design = ~ sample + exon + condition:exon , flattenedfile = NULL) 

exptData(dxd)
head( counts(dxd), 5 )


formulaFullModel = ~ sample + exon + condition:exon
formulaReducedModel =  ~ sample + exon 

allDEXSeq[[model]] <- DEXSeq(dxd, fullModel = formulaFullModel, reducedModel = formulaReducedModel, BPPARAM=BPPARAM, fitExpToVar="condition")





### table of results

models <- names(allDEXSeq)
analysis <- paste0( method_out ,"/diff_out_", models)
analysis


for(i in 1:length(analysis)){
  # i=1
  dir.create(analysis[i], showWarnings=FALSE, recursive=TRUE)
  dxr <- allDEXSeq[[i]]
  
  dxrf <- as.data.frame(dxr)
  
  write.table(dxrf[, -ncol(dxrf)], paste0(analysis[i], "/DEXSeq_", count.method, "_exon_results.txt"), sep="\t", row.names=FALSE, quote=F) 
  
  ### gene level test
  perGeneRes <- perGeneQValue(dxr)
  
  write.table(data.frame(geneID=names(perGeneRes), pvalue=perGeneRes), paste0(analysis[i], "/DEXSeq_", count.method, "_gene_results.txt"), sep="\t", row.names=FALSE, quote=F) 
  
  
}












# When using kallisto counts

# Warning messages:
#   1: glm.fit: algorithm did not converge
# 2: glm.fit: algorithm did not converge
# 3: In estimateDispersionsFit(object, fitType = fitType, quiet = quiet) :
#   the parametric fit of dispersion estimates over the mean of counts
# failed, which occurs when the trend is not well captured by the
# function y = a/x + b. A local regression fit is automatically performed,
# and the analysis can continue. You can specify fitType='local' or 'mean'
# to avoid this message if re-running the same data.
# When using local regression fit, the user should examine plotDispEsts(dds)
# to make sure the fitted line is not sharply curving up or down based on
# the position of individual points.













































