########################################################
# Casper transcript counting

# Created 24 Oct 2014
# BioC 3.0

########################################################

library(casper)

library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)

setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2")


#### find junctions

# gtf must be with chr
gtfFile <- "/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.chr.gtf"


genDB <- makeTranscriptDbFromGFF(file=gtfFile, format="gtf", exonRankAttributeName="exon_number")
genomeDBTmp <- genomeDB <- procGenome(genDB=genDB, genome="Dm70", mc.cores=10)


genDB2 <- import(gtfFile) ## somehow does not work
genomeDB2 <- procGenome(genDB2, genome="Dm70", mc.cores=10)


#### load bams (sorted & indexed)

name=""
mc.cores=20

library(Rsamtools)
what <- scanBamWhat()
what <- what[!(what %in% c("seq","qual"))]
flag <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE)
param <- ScanBamParam(flag=flag, what=what, tag="XS")


bamFiles <- paste0("/home/Shared_taupo/tmp/Simulation/simulation_drosophila_V2/1_reads/bam/sim_samp", 1:6, "_s.bam")
names(bamFiles) <- c(paste0("C1", "S", 1:3), paste0("C2", "S", 1:3))



esets <- lapply(bamFiles, function(bamFile){
  # bamFile <- bamFiles[1]
  
  bam0 <- scanBam(file=bamFile, param=param)[[1]]
  
  bam0 <- rmShortInserts(bam0, isizeMin=50)
  
  pbam0 <- procBam(bam0)
  #  head(getReads(pbam0))
  
  distrs <- getDistrs(genomeDB, bam=bam0, readLength=75)

  ### path counting
  pc <- pathCounts(pbam0, DB=genomeDB, mc.cores = mc.cores)
  # head(pc@counts[[1]])

  eset <- calcExp(distrs=distrs, genomeDB=genomeDB, pc=pc, readLength=75, rpkm=TRUE, mc.cores = mc.cores)
#   head(exprs(eset))
# head(fData(eset))
  
  cat(paste("Done for \n", bamFile, "\n"))
  
return(eset)
  
})

dir.create("casper_2.0.0")

load("casper_2.0.0/esets.RData")


for(i in 1:length(esets)){
  # i = 1 
  eset <- esets[[i]]
  
  rpkm <- 2^exprs(eset)
  colnames(rpkm) <- paste0("rpkm",i)
  counts <- fData(eset)
  
  out.table <- merge(counts, rpkm, by=0, all=TRUE)

  if(i == 1)
    rpkmAll <- out.table[,c("transcript", "gene_id", "rpkm1")]

  
  if( i > 1)
  rpkmAll <- merge(rpkmAll, rpkm, by.x="transcript", by.y=0)
  
  write.table(out.table[,-1], paste0("casper_2.0.0/casper_counts_rpkm_samp", i ,".xls"), row.names=F, sep="\t")
 
  
}

write.table(rpkmAll, paste0("casper_2.0.0/casper_counts_rpkm.xls"), row.names=F, sep="\t")

save(esets, file = "casper_2.0.0/esets.RData")







### compare with simulations



casper <- rpkmAll
# casper[,3:8] <- round(casper[,3:8])

simuInfo <- read.table("/home/gosia/Multinomial_project/Simulations_drosophila_V2/Simu_info/true_genes_simulation.txt", header = TRUE, stringsAsFactors = FALSE)

width <- simuInfo[,"width"]

simuInfo <- simuInfo[,c("transcript_id", paste0("rpk", 1:6))]

simuInfo[, paste0("rpk", 1:6)] <- simuInfo[, paste0("rpk", 1:6)] / (width/1000) / (colSums(simuInfo[, paste0("rpk", 1:6)])/10^6) ### Convert to RPKM 


casperSimu <- merge(simuInfo, casper, by.x = "transcript_id", by.y = "transcript")


pdf(paste0("casper_2.0.0/CaspervsSimulations.pdf"))

for(i in 1:6){
  # i=1
  
  smoothScatter(log(casperSimu[, paste0("rpk",i)]+1), log(casperSimu[, paste0("rpkm", i)]+1), nrpoints = Inf, main=paste("Sample", i), xlab="log simulations RPKM ", ylab="log casper RPKM", colramp = colorRampPalette(c("white", "cyan4")))
  abline(a=0, b=1, col="brown3")
  
}

dev.off()































