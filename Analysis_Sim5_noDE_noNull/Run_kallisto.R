# BioC 3.1

# Crated 17 June 2015



#########################################################################################
# Building an index
#########################################################################################

# kallisto index -i transcripts.idx transcripts.fasta.gz

cDNA.fasta <- "/home/Shared/data/annotation/Drosophila/Ensembl70/cDNA/Drosophila_melanogaster.BDGP5.70.cdna.all.fa"
index <- "/home/Shared/data/annotation/Drosophila/Ensembl70/kallisto/Drosophila_melanogaster.BDGP5.70.cdna.all.idx"
dir.create(dirname(index), recursive = TRUE, showWarnings = FALSE)


cmd <- paste("kallisto index -i",  index, cDNA.fasta, sep = " ")

system(cmd)



#########################################################################################
# Quantification
#########################################################################################

# kallisto quant -i transcripts.idx -o output -b 100 -t 10 reads_1.fastq.gz reads_2.fastq.gz
# kallisto quant -i index -o output pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq


out.path <- paste0("/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull/2_counts/kallisto/sample_", 1:6, "/")


fastq.path <- "/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull/1_reads/reads/"

fastq <- paste0(fastq.path, "sample", 1:6, "/sample_", 1:6 ,"_1.fq ", fastq.path, "sample", 1:6, "/sample_", 1:6 ,"_2.fq")



for(i in 2:6){
	# i = 1
	dir.create(out.path[i], recursive = TRUE, showWarnings = FALSE)
	cmd <- paste("kallisto quant -i", index, "-o", out.path[i], "-b 1 -t 10", fastq[i])	
	
	system(cmd)
	
}





#########################################################################################
# Adjust the abundance.txt file to DM needs 
#########################################################################################

library(rtracklayer)

gtf.path = "/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf"

gtf <- import(gtf.path)

geneTrans <- unique(mcols(gtf)[, c("gene_id", "transcript_id")])
rownames(geneTrans) <- geneTrans$transcript_id



############## save results in tables

out.path.tmp <- "/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull/2_counts/kallisto/"

countsList <- lapply(1:6, function(i){
  # i = 1

	abundance <- read.table(paste0(out.path[i], "abundance.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  
  counts <- data.frame(gene_id = geneTrans[abundance$target_id, "gene_id"], transcript_id = abundance$target_id, counts = round(abundance$est_counts), stringsAsFactors = FALSE)
	
  colnames(counts) <- c("gene_id", "transcript_id", paste0("sample_", i))
  
  counts_dexseqLike <- data.frame(group_id = paste0(counts[, "gene_id"], ":", counts[, "transcript_id"]), counts[, paste0("sample_", i)], stringsAsFactors = FALSE)
  
	colnames(counts_dexseqLike) <- c("group_id", paste0("sample_", i))
	
  write.table(counts_dexseqLike, paste0(out.path.tmp, "kallisto", i, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(counts_dexseqLike)
  
})
  


counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), countsList)


write.table(counts, paste0(out.path.tmp, "kallisto_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)












































