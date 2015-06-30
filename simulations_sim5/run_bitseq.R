# BioC 3.0

# Crated 16 June 2015



#########################################################################################
# Align reads to transcriptom with bowtie
#########################################################################################



# create bowtie reference index for the annotation
# $ bowtie-build -f --ntoa ensSelect1.fasta ensSelect1-index

cDNA.fasta <- "/home/Shared/data/annotation/Drosophila/Ensembl70/cDNA/Drosophila_melanogaster.BDGP5.70.cdna.all.fa"
index <- "/home/Shared/data/annotation/Drosophila/Ensembl70/bowtie_index_cDNA/Dme_BDGP5.70-index"
dir.create(dirname(index), recursive = TRUE, showWarnings = FALSE)


cmd <- paste("bowtie-build -f --ntoa ", cDNA.fasta, index , sep = " ")

system(cmd)



# align reads in data-c0b0.fastq against index
# $ bowtie -q -v 3 -3 0 -p 4 -a -m 100 --sam ensSelect1-index data-c0b0.fastq data-c0b0.sam

fastq.path <- "/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull/1_reads/reads/"

fastq1 <- paste0(fastq.path, "sample", 1:6, "/sample_", 1:6 ,"_1.fq")
fastq2 <- paste0(fastq.path, "sample", 1:6, "/sample_", 1:6 ,"_2.fq")


out.path <- "/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull/1_reads/Bowtie/"
dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

for(i in 1:6){
  # i = 1
  cmd <- paste0("bowtie -q -v 3 -3 0 -p 10 -a -m 100 --sam ", index, " -1 ", fastq1[i], " -2 ", fastq2[i], " ", out.path, "sample_", i ,".sam" )
  cat(cmd, fill = TRUE)
  
  system(cmd)
  
}



#########################################################################################
# Get the transcript expression with BitSeq getExpression
#########################################################################################

library(BitSeq)

out.path <- "/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull/2_counts/BitSeq_1.10.0/"
dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

setwd(out.path)

trSeqFile <- cDNA.fasta 

out <- mclapply(1:6, function(i){
  # i = 1
  
  alignFile <- paste0("/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull/1_reads/Bowtie/", "sample_", i ,".sam")

  res <- getExpression(alignFile, trSeqFile, outPrefix = paste0("sample_", i ), type="RPKM", log=FALSE, seed=47)
   
  save(res, file = paste0("BitSeq_results_sample_", i, ".RData"))
  
  return(NULL)
  
}, mc.cores = 6)


############## save results in tables



countsList <- lapply(1:6, function(i){
  
  load(paste0("BitSeq_results_sample_", i, ".RData"))
  
  counts <- data.frame(res$trInfo[, c("gene", "transcript")], counts = res$counts, stringsAsFactors = FALSE)
  colnames(counts) <- c("gene_id", "transcript_id", paste0("sample_", i))
  
  counts_dexseqLike <- data.frame(group_id = paste0(counts[, "gene_id"], ":", counts[, "transcript_id"]), counts[, paste0("sample_", i)], stringsAsFactors = FALSE)
  
  write.table(counts_dexseqLike, paste0("bitseq", i, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(counts_dexseqLike)
  
})
  


counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), countsList)


write.table(counts, paste0("bitseq_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)












































