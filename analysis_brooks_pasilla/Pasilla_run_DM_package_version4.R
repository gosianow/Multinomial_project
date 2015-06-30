#######################################################
# 
# Created 04 Sep 2014 

# comprare the DM version 4 commonDispersion and adjustement performance with other methods and DM dirmult

# Update 04 Sep 2014:



#######################################################
# BioC 2.14

setwd("/home/gosia/Multinomial_project/RealData_brooks_pasilla/")

# create metadata file
sri.org <- read.table("/home/Shared/data/seq/brooks_pasilla/download/SraRunInfo.csv", stringsAsFactors=F, sep=",", header=T)
keep <- grep(paste("GSM4611", 76:82, sep="", collapse="|"), sri.org$SampleName)
sri <- sri.org[keep,]

# manual trimming of reads to have equal lenght 
which.trim <- sri[sri$SampleName=="GSM461179", c("Run")]

# cmd <- paste0("java -jar /usr/local/software/Trimmomatic-0.30/trimmomatic-0.30.jar  SE -threads 20 download/", which.trim , ".fastq.gz download/", which.trim , ".trimed.fastq.gz CROP:40 \n")
# cat(cmd)
# for(i in 1:length(cmd))
#   system(cmd[i])

sri[sri$Run %in% which.trim, "avgLength"] <- 40
sri[sri$Run %in% which.trim, "Run"] <- paste0(sri[sri$Run %in% which.trim, "Run"], ".trimed")

sri$LibraryName = gsub("S2_DRSC_","",sri$LibraryName) # trim label
metadata = unique(sri[,c("LibraryName","LibraryLayout", "SampleName", "avgLength" )])

for(i in seq_len(nrow(metadata))) {
  rw = (sri$LibraryName==metadata$LibraryName[i])
  if(metadata$LibraryLayout[i]=="PAIRED") {
    metadata$fastq1[i] = paste0("download/",sri$Run[rw],"_1.fastq.gz",collapse=",")
    metadata$fastq2[i] = paste0("download/",sri$Run[rw],"_2.fastq.gz",collapse=",")
    metadata$ReadLength[i] <- metadata$avgLength[i] / 2
    metadata$SRR[i] <- paste0(sri$Run[rw], collapse=",")
  } else {
    metadata$fastq1[i] = paste0("download/",sri$Run[rw],".fastq.gz",collapse=",")
    metadata$fastq2[i] = ""
    metadata$ReadLength[i] <- metadata$avgLength[i]
    metadata$SRR[i] <- paste0(sri$Run[rw], collapse=",")
  }
}

metadata$condition = "CTL"
metadata$condition[grep("RNAi",metadata$LibraryName)] = "KD"
metadata$shortname = paste( seq_len(nrow(metadata)), substr(metadata$condition,1,2),  substr(metadata$LibraryLayout,1,2), metadata$ReadLength, sep=".")

metadata$color[metadata$condition == "CTL"] <- "chartreuse3"
metadata$color[metadata$condition == "KD"] <- "darkorchid3"

metadata.org <- metadata



##########################################################################
# run new DM (common dispersion, and adjustement)  on FC
##########################################################################

library(edgeR)
library(parallel)
library(dirmult)
library(limma)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version4/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v4.R"))





out.dir <- "DM_v4/fc/"
dir.create(out.dir, showWarnings=F, recursive=T)


#### load FC

fc <- read.table("featureCounts/fc.txt", header = T)

colnames(fc) <- gsub("_s.bam","",gsub("sam.", "",colnames(fc)))

counts <- fc[, metadata$SampleName]
gene.id <- strsplit2(fc[,1], ":")[,1]
ete.id <- fc[,1]


dgeOrg <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))

name1 <- "fc"
dge <- dgeOrg
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)



###############################################################
name1 <- "fc_g0_s4_keep0s_subsetInf_DM4adj"
dge <- dgeOrg
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$genes$gene.id <- as.factor(as.character(dge$genes$gene.id))
nlevels(dge$genes$gene.id)

# group=NULL; adjust = TRUE; mode = "constrOptim2"; epsilon = 1e-05; maxIte = 1000; interval = c(0, 1e+5); tol = 1e-03; mcCores=30; verbose=FALSE

## run DM pipeline
dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=30, verbose=FALSE)

write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




#######################################################################################################
# MDS plot
#######################################################################################################


pdf(paste0(out.dir, "/",name1,"_MDS.pdf"))

mds <- plotMDS(dge, col=metadata$color, top=500, main="Method=logFC, prior.count=2", labels=metadata$shortname, cex=2, cex.lab=1.45, cex.axis=1.5, cex.main=1.5, xlim=c(-2, 2), ylim=c(-2, 2))
mds10 <- plotMDS(dge,col=metadata$color,main="Method=logFC, prior.count=10", prior.count=10, top=500, labels=metadata$shortname, cex=2, cex.lab=1.45, cex.axis=1.5, cex.main=1.5, xlim=c(-2, 2), ylim=c(-2, 2))
mdsb <- plotMDS(dge, col=metadata$color, main="Method=bcv", method="bcv", top=500, labels=metadata$shortname, cex=2, cex.lab=1.45, cex.axis=1.5, cex.main=1.5, xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5))

dev.off()



###################################################
# create annotation table
###################################################

out.dir <- "gene_info"
dir.create(out.dir, showWarnings=F, recursive=T)

library(rtracklayer)

gtfFile <- "/home/Shared_penticton/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf"

gtf0 <- import(gtfFile, asRangedData=FALSE)

idx <- mcols(gtf0)$type == "exon" 
gtf0 <- gtf0[idx]

ann <- data.frame(GeneID = mcols(gtf0)$gene_id, GeneName=mcols(gtf0)$gene_name)
ann.u <- unique(ann)

# library(limma)
# ann.u$gene_symbol <- alias2SymbolTable(alias=ann.u$gene_name, species = "Dm")

write.table(ann.u, paste0(out.dir, "/","Gene_symbols.xls"), quote=FALSE, sep="\t", row.names=FALSE)

### does not work
# library(org.Dm.eg.db)
# get( "FBgn0039634" , org.Dm.egGENENAME)


###################################################
# check the results
###################################################

### add gene 
table <- dgeDM$table
table <- merge(table, ann.u, by="GeneID", all.x=TRUE)


out.dir <- "DM_v4/fc/"
dir.create(out.dir, showWarnings=F, recursive=T)

write.table(table, paste0(out.dir, "/",name1,"_resultsGeneNames.xls"), quote=FALSE, sep="\t", row.names=FALSE)














