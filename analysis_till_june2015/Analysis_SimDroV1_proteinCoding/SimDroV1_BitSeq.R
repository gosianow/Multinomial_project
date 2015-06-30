
# Crated 27 Jan 2015
# Midified

# BioC 3.0


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

fastq.path <- "/home/Shared/tmp/Simulation2/simulation_drosophila_V1/1_reads/created_files/"

fastq1 <- paste0(fastq.path, "samp", 1:6, "/sim_1.fastq.gz")
fastq2 <- paste0(fastq.path, "samp", 1:6, "/sim_2.fastq.gz")


out.path <- "/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding/Fastq/"
dir.create(out.path, recursive = TRUE, showWarnings = FALSE)


fastq1.new <- paste0(out.path, "samp", 1:6, "/sim_1.fastq")
fastq2.new <- paste0(out.path, "samp", 1:6, "/sim_2.fastq")


for(i in 2:6){
  # i = 1
  dir.create(dirname(fastq1.new[i]), recursive = TRUE, showWarnings = FALSE)
  
  cmd <- paste0("gunzip -c ", fastq1[i], " > ", fastq1.new[i])
  cat(cmd, fill = TRUE)
  system(cmd)
  
  cmd <- paste0("gunzip -c ", fastq2[i], " > ", fastq2.new[i])
  cat(cmd, fill = TRUE)
  system(cmd)
  
}


out.path <- "/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding/Bowtie/"
dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

for(i in 2:6){
  # i = 1
  cmd <- paste0("bowtie -q -v 3 -3 0 -p 10 -a -m 100 --sam ", index, " -1 ", fastq1.new[i], " -2 ", fastq2.new[i], " ", out.path, "samp", i ,".sam" )
  cat(cmd, fill = TRUE)
  
  system(cmd)
  
}



#########################################################################################
# Get the transcript expression with BitSeq getExpression
#########################################################################################

library(BitSeq)

out.path <- "/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding/BitSeq_1.10.0/"
dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

setwd(out.path)

trSeqFile <- cDNA.fasta 

out <- mclapply(1:6, function(i){
  # i = 1
  
  alignFile <- paste0("/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding/Bowtie/", "samp", i ,".sam")

  res <- getExpression(alignFile, trSeqFile, outPrefix = paste0("samp", i ), type="RPKM", log=FALSE, seed=47)
   
  save(res, file = paste0("BitSeq_results_samp", i, ".RData"))
  
  gc()
  
  return(NULL)
  
}, mc.cores = 6)


############## save results in tables

counts <- list()

for(i in 1:6){
  # i = 1
  
  load(file = paste0("BitSeq_results_samp", i, ".RData"))
  
  counts[[i]] <- data.frame(res$trInfo[, c("gene", "transcript")], counts = res$counts, stringsAsFactors = FALSE)
  
  colnames(counts[[i]]) <- c("gene", "transcript", paste0("counts", i))

  write.table(data.frame(paste0(counts[[i]][, "gene"], ":",counts[[i]][, "transcript"]), counts[[i]][, paste0("counts", i)]), paste0("BitSeq_counts",i,".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

  
}


counts.m <- counts[[1]]

for(i in 2:6){  
  counts.m <- merge(counts.m, counts[[i]][,c("transcript", paste0("counts", i))], by = "transcript", all = TRUE, sort = FALSE)
}

write.table(counts.m, paste0("BitSeq_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)





#########################################################################################
### compare BitSeq expression with simulations
#########################################################################################


simuInfo <- read.table("/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding/Simu_info/true_genes_simulation.txt", header = TRUE, as.is = TRUE)

simuInfo <- simuInfo[,c("transcript_id", paste0("rpk", 1:6))]

BitSeqSimu <- merge(simuInfo, counts.m, by.x = "transcript_id", by.y = "transcript")


pdf(paste0(out.path, "BitSeqvsSimulations.pdf"))

for(i in 1:6){
  # i=1
  
  smoothScatter(log(BitSeqSimu[, paste0("rpk",i)]+1), log(BitSeqSimu[, paste0("counts", i)]+1), nrpoints = Inf, main=paste("Sample", i), xlab="log simulations", ylab="log BitSeq expected counts")
  abline(a=0, b=1, col="magenta")
  
}

dev.off()









































