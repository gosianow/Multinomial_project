# 18 June 2014 (Last updated 23 June 2014)
library(limma)
files <- paste0(path, list.files(path, "results.txt"))

# load cuffdiff results
cuffdiff.iso <- read.table("./simulation_hsapiens/cuffdiff_results.txt", header=TRUE)
head(cuffdiff.iso)
cuffdiff.iso <- read.columns(files[grep("cuffdiff", files)], 
	required.col=c("p_value", "q_value", "transcript_id", "ensembl_gene_id", "TSS_group_id"))
cuffdiff.iso <- cuffdiff.iso[,c(5,4,1,2,3)]
colnames(cuffdiff.iso) <- c("GeneID", "TranscriptID", "TSSID", "P.Value", "FDR")
# o <- order(cuffdiff.iso$FDR)
# cuffdiff.gene <- cuffdiff.iso[o,]
# cuffdiff.gene <- cuffdiff.gene[!duplicated(cuffdiff.gene$GeneID),]

# load dexseq results
dexseq.exon <- read.columns(files[grep("dexseq_exon", files)], 
	required.col=c("groupID", "featureID", "pvalue", "padj", "log2fold_Condition2_Condition1"))
dexseq.exon$featureID <- gsub(":E", ":", dexseq.exon$groupID)
dexseq.exon$groupID <- strsplit2(dexseq.exon$groupID, split=":")[,1]
colnames(dexseq.exon) <- c("GeneID", "ExonID", "P.Value", "FDR", "logFC")
dexseq.gene <- read.table(files[grep("dexseq_gene", files)], header=TRUE)
colnames(dexseq.gene) <- c("GeneID", "FDR")

# load voomex results
voom.exon <- read.table(files[grep("voom_exon", files)], header=TRUE)[,-3]
colnames(voom.exon)[1:2] <- c("GeneID", "ExonID")
voom.gene <- read.table(files[grep("voom_gene", files)], header=TRUE)[,-1]

# load annotation for ds exons for nonnull simulations
if (length(grep("null", sim))==0) {
	annofile <- paste0(path, list.files(path, "true_exons"))
	annofile <- annofile[grep(".txt", annofile)]
	ds.anno <- read.columns(annofile, required.col=c("Gene", "Transcript", "exon", "status", "status_exon"))
	ds <-  list(Gene=NULL, Transcript=NULL, Exon=NULL)
	ds$Gene <- ds.anno$Gene[!is.na(ds.anno$Gene)]
	ds$Gene <- unique(as.character(ds$Gene))
	ds$Transcript <- ds.anno$Transcript[!is.na(ds.anno$Transcript)]
	ds$Transcript <- unique(as.character(ds$Transcript))
	ds$Exon <- ds.anno$exon[ds.anno$status_exon==1]
	ds$Exon <- as.character(ds$Exon[!is.na(ds$Exon)])
	if (sum(duplicated(ds$Exon))!=0) cat("WARNING: Non-unique exon IDs in annotation of DS features", "\n")
	if (length(ds$Gene)!=1000) cat("WARNING: Incorrect no. of DS genes", "\n")
	if (length(ds$Transcript)!=2000) cat("WARNING: Incorrect no. of DS transcripts", "\n")
	if (length(ds$Exon)<length(ds$Gene)) cat("WARNING: Fewer DS exons than DS genes", "\n")
}
