######################################################
# BioC 2.14
# Created 04 Nov 2014 

# Run featureCounts and compare to htseq counts

#######################################################



setwd("/home/Shared/data/seq/Kim_adenocarcinoma/")


metadata <- read.table("metadata/Malgorzata Nowicka2014-11-04GSE37764.csv", stringsAsFactors=F, sep=",", header=T) 

metadataOrg <- metadata <- metadata[metadata$X == "RNA-seq",]


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

























































