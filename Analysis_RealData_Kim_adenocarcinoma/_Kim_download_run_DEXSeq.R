######################################################

# Created 04 Nov 2014 

# Analysis of Kim_adenocarcinoma data

#######################################################
# BioC 2.14


setwd("/home/Shared/data/seq/Kim_adenocarcinoma/")


metadata <- read.table("metadata/Malgorzata Nowicka2014-11-04GSE37764.csv", stringsAsFactors=F, sep=",", header=T) 

metadataOrg <- metadata <- metadata[metadata$X == "RNA-seq",]




##############################################################################
# download data
##############################################################################



#### download tophat results from insilicodb ####


link <- paste0("https://insilicodb.com/app/publicutilities/getfile?seriesName=GSE37764&platformName=GPL10999&file=source/GSE37764-measurements/", metadata$ids ,"/rnaseq/v0.9/tophat_out/accepted_hits.bam")
outPath <- paste0("tophat_insilicodb/", metadata$ids ,"/accepted_hits.bam")



for(i in 1:nrow(metadata)){
  
  dir.create(paste0("tophat_insilicodb/", metadata$ids[i] ,"/"), recursive = TRUE)
  download.file(link[i], outPath[i], method="wget")
  
}



#### organize BAM files with samtools ####

dir.create(paste0("bam_insilicodb/"), recursive = TRUE)

for (i in 1:nrow(metadata)){
  # i=1
  bam.in <- paste0("tophat_insilicodb/" , metadata$ids[i], "/accepted_hits.bam")
  bam.out <- paste0("bam_insilicodb/" ,metadata$ids[i])
  
  sort by position (for featureCounts)
  convert to SAM (for DEXSeq) not needed 
    cmdt = paste0("samtools sort ", bam.in, " ", bam.out, "_s", "\n",
                  "samtools view -o ", bam.out, "_s.sam ", bam.out, "_s.bam", "\n" )
  
#   sort by name (for DEXSeq)
#   cmdt <- paste0("samtools sort -n ", bam.in, " ", bam.out, "_sn", "\n")
#   
#   add index for IGV
#   cmdt <- paste0("samtools index ", bam.out, "_s.bam", "\n")
  
  
  cat(cmdt)
  system(cmdt)
  
}







##############################################################################
# reads counting - featureCounts with flattened gtf / wrapper_featureCounts_v3
##############################################################################



source("/home/gosia/R/R_Counting/featureCounts/wrapper_featureCounts_v3.R")


gtfFile <- "/home/Shared_penticton/data/annotation/Human/Ensembl_GRCh37.71/gtf/Homo_sapiens.GRCh37.71.gtf"
out.dir <- "featureCounts"


PE_bam_files <-  paste0("/home/Shared/data/seq/Kim_adenocarcinoma/bam_insilicodb/", metadata$ids, "_s.bam")
PE_bam_files


fc <- wrapper_featureCounts(gtfFile, PE_bam_files=PE_bam_files,  nthreads=20)


###################################################
# DEXSeq
###################################################
cat("Running DEXSeq \n")

library("DEXSeq")

out.dir <- "DEXSeq_1.10.8"
dir.create(out.dir)
dir.create("DEXSeq_1.10.8/Exon_counts", recursive=TRUE, showWarnings=FALSE)


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

python.cmd2 <- with(metadata, paste0("python ", pythDir, "/dexseq_count.py -p yes -s no -f bam -r pos ", DEXSeq.gff, " bam_insilicodb/", ids, "_s.bam ", out.dir, "/Exon_counts/", ids , ".counts \n"))
cat(python.cmd2)

 for(i in 1:nrow(metadata)){
   # i = 1 
   system(python.cmd2[i])
   
 }





### DEXSeq
library(BiocParallel)
BPPARAM = MulticoreParam(workers=10)
allDEXSeq <- list()



# Null_normal1
model <- "Null_normal1"
dir.create(paste0(out.dir, "/diff_out_", model))

metadata <-  metadataOrg[metadataOrg$Tissue.Type == "normal", ]
countFiles <- paste0(out.dir, "/Exon_counts/", metadata$ids, ".counts")

sampleTable = data.frame(row.names = metadata$ids, condition = c(rep("C1", 3), rep("C2", 3)))

dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design = ~ sample + exon + condition:exon , flattenedfile = DEXSeq.gff ) 

exptData(dxd)

formulaFullModel = ~ sample + exon + condition:exon
formulaReducedModel =  ~ sample + exon 

allDEXSeq[[model]] <- DEXSeq(dxd, fullModel = formulaFullModel, reducedModel = formulaReducedModel, BPPARAM=BPPARAM, fitExpToVar="condition")




# Null_tumor1
model <- "Null_tumor1"
dir.create(paste0(out.dir, "/diff_out_", model))

metadata <-  metadataOrg[metadataOrg$Tissue.Type == "tumor", ]
countFiles <- paste0(out.dir, "/Exon_counts/", metadata$ids, ".counts")

sampleTable = data.frame(row.names = metadata$ids, condition = c(rep("C1", 3), rep("C2", 3)))

dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData=sampleTable, design = ~ sample + exon + condition:exon , flattenedfile = DEXSeq.gff ) 

exptData(dxd)

formulaFullModel = ~ sample + exon + condition:exon
formulaReducedModel =  ~ sample + exon 

allDEXSeq[[model]] <- DEXSeq(dxd, fullModel = formulaFullModel, reducedModel = formulaReducedModel, BPPARAM=BPPARAM, fitExpToVar="condition")





### table of results

models <- names(allDEXSeq)
analysis <- paste0( out.dir ,"/diff_out_", models)

for(i in 1:length(analysis)){
  # i=1
  dir.create(analysis[i], showWarnings=FALSE, recursive=TRUE)
  dxr <- allDEXSeq[[i]]
  
  dxrf <- as.data.frame(dxr)
  
  write.table(dxrf[, -ncol(dxrf)], paste0(analysis[i], "/DEXSeq_exon_results.txt"), sep="\t", row.names=FALSE, quote=F) 
  
  ### gene level test
  perGeneRes <- perGeneQValue(dxr)
  
  write.table(data.frame(geneID=names(perGeneRes), pvalue=perGeneRes), paste0(analysis[i], "/DEXSeq_gene_results.txt"), sep="\t", row.names=FALSE, quote=F) 
  
  
}




###################################################
# compare featureCounts with HTSeq
###################################################


#### check with DEXSeq
DS <- read.table("DEXSeq_1.10.8/Exon_counts/GSM927308.counts")
colnames(DS) <- c("exons", "DEXSeq")

fc <- read.table("featureCounts/fc.txt", header=T)
FC <- data.frame(exons=fc[,"exon_id"], featureCounts = fc[, grep("GSM927308", colnames(fc))])

all( FC$exons %in% DS$exons )

DS [ which (! DS$exons %in% FC$exons ) , ]


counts <- merge(FC, DS, by="exons", all.x=TRUE)


pdf(paste0("DEXSeq_1.10.8/DEXSeq_vs_featureCounts.pdf"))
smoothScatter(counts[, "DEXSeq"], counts[, "featureCounts"], xlim=c(0, 1e3), ylim=c(0, 1e3), nrpoints = Inf)
dev.off()

























































