########################################################
# Run RSEM

# Created 24 Oct 2014
# BioC 3.0

########################################################

setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2/RSEM/")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


########################################################
### I. Preparing Reference Sequences
########################################################

### rsem-prepare-reference [options] reference_fasta_file(s) reference_name

reference_fasta_file <- "/home/Shared/data/annotation/Drosophila/Ensembl70/cDNA/Drosophila_melanogaster.BDGP5.70.cdna.all.fa"
reference_name <- "/home/Shared/data/annotation/Drosophila/Ensembl70/RSEM/Dm70"

cmd <-  paste("rsem-prepare-reference ", "--bowtie2", reference_fasta_file, reference_name)


system(cmd)


########################################################
### II. Calculating Expression Values
########################################################

### rsem-calculate-expression [options]  upstream_read_file(s) downstream_read_file(s) reference_name sample_name 


upstream_read_files <- paste0("/home/Shared/tmp/Simulation/simulation_drosophila_V2/created_files/samp", metadata$SampleName1, "/sim_1.fastq")
downstream_read_files <- paste0("/home/Shared/tmp/Simulation/simulation_drosophila_V2/created_files/samp", metadata$SampleName1, "/sim_2.fastq")
reference_name <- "/home/Shared/data/annotation/Drosophila/Ensembl70/RSEM/Dm70"
sample_name <- paste0("samp",metadata$SampleName1)


cmd <- paste("rsem-calculate-expression", "--paired-end", "--bowtie2", "--num-threads 10", upstream_read_files, downstream_read_files, reference_name, sample_name)



for(i in 1:length(cmd))
  system(cmd[i])




########################################################
### Merge results
########################################################



for(i in 1:6){
  
  res <- read.table(paste0("samp",i,".isoforms.results"), header = TRUE)
  print(head(res))
  
  if(i == 1){
    RSEM_FPKM <- res[, c("transcript_id", "length")]
    RSEM_exp_count <- res[, c("transcript_id", "length")]
  }
  
  
  RSEM_FPKM <- merge(RSEM_FPKM, res[, c("transcript_id", "FPKM")], by = "transcript_id")
  
  RSEM_exp_count <- merge(RSEM_exp_count, res[,c("transcript_id", "expected_count")], by = "transcript_id")
  
  
}


colnames(RSEM_FPKM) <- c("transcript_id", "length", paste0("FPKM", 1:6))

colnames(RSEM_exp_count) <- c("transcript_id", "length", paste0("expected_count", 1:6))


### add gene_id


library(GenomicRanges)
library(rtracklayer)

gtfFile <- "/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf"

gtf0 <- import(gtfFile)

t2g <- unique(mcols(gtf0)[,c("gene_id", "transcript_id")])


RSEM_FPKM <- merge(t2g, RSEM_FPKM, by = "transcript_id")
RSEM_exp_count <- merge(t2g, RSEM_exp_count, by = "transcript_id")



### save results

write.table(RSEM_FPKM, "RSEM_FPKM.xls", quote = F, sep = "\t", row.names = FALSE)

write.table(RSEM_exp_count, "RSEM_exp_count.xls", quote = F, sep = "\t", row.names = FALSE)


### compare with simulations



rsem <- data.frame(RSEM_exp_count)
rsem[,4:9] <- round(rsem[,4:9])

simuInfo <- read.table("/home/gosia/Multinomial_project/Simulations_drosophila_V2/Simu_info/true_genes_simulation.txt", header = TRUE, stringsAsFactors = FALSE)

simuInfo <- simuInfo[,c("transcript_id", paste0("rpk", 1:6))]

rsemSimu <- merge(simuInfo, rsem, by = "transcript_id")


pdf(paste0("RSEMvsSimulations.pdf"))

for(i in 1:6){
  # i=1
  
  smoothScatter(log(rsemSimu[, paste0("rpk",i)]+1), log(rsemSimu[, paste0("expected_count", i)]+1), nrpoints = Inf, main=paste("Sample", i), xlab="log simulations", ylab="log RSEM expected counts")
  abline(a=0, b=1, col="magenta")
  
#     smoothScatter(rsemSimu[, paste0("rpk",i)], rsemSimu[, paste0("expected_count", i)], xlim=c(0, 1e3), ylim=c(0, 1e3), nrpoints = Inf, main=paste("Sample", i), xlab="Simulations", ylab="RSEM expected counts")
#     abline(a=0, b=1, col="red")
#   
}

dev.off()






rsem <- data.frame(RSEM_FPKM)
rsem[,4:9] <- round(rsem[,4:9])

simuInfo <- read.table("/home/gosia/Multinomial_project/Simulations_drosophila_V2/Simu_info/true_genes_simulation.txt", header = TRUE, stringsAsFactors = FALSE)

simuInfo <- simuInfo[,c("transcript_id", paste0("rpk", 1:6))]

rsemSimu <- merge(simuInfo, rsem, by = "transcript_id")


pdf(paste0("RSEMvsSimulations_FPKM.pdf"))

for(i in 1:6){
  # i=1
  
  smoothScatter(log(rsemSimu[, paste0("rpk",i)]+1), log(rsemSimu[, paste0("FPKM", i)]+1), nrpoints = Inf, main=paste("Sample", i), xlab="log simulations", ylab="log RSEM FPKM")
  abline(a=0, b=1, col="red")
  
  #     smoothScatter(rsemSimu[, paste0("rpk",i)], rsemSimu[, paste0("FPKM", i)], xlim=c(0, 1e3), ylim=c(0, 1e3), nrpoints = Inf, main=paste("Sample", i), xlab="Simulations", ylab="RSEM FPKM")
  #     abline(a=0, b=1, col="red")
  #   
}

dev.off()




























