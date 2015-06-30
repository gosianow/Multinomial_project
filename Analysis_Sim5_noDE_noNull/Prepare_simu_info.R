##############################################################################

# BioC 3.1
# Created 16 June 2015:

##############################################################################

### drosophila
# setwd("/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull")
# DEXSeq_gff_path = "/home/Shared/data/annotation/Drosophila/Ensembl70/DEXSeq_1.10.8_gff/Drosophila_melanogaster.BDGP5.70.DEXSeq.flattened.rNO.gff"

### hsapiens noDE
setwd("/home/gosia/Multinomial_project/Simulations_hsapiens_Sim5_noDE_noNull")
DEXSeq_gff_path = "/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/hsapiens/reference_files/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff"


### hsapiens withDE
setwd("/home/gosia/Multinomial_project/Simulations_hsapiens_Sim5_withDE_noNull")
DEXSeq_gff_path = "/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/hsapiens/reference_files/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.nomerge.gff"







### simulation_details

simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, as.is = TRUE)

simulation_details <- simulation_details[, c("gene_id", "transcript_id", "expected_count", "gene_ds_status", "transcript_ds_status", "gene_de_status")]

simulation_details_ds <- simulation_details[simulation_details$transcript_ds_status == 1, ]

### DEXSeq gtf file 
library(rtracklayer)

DEXSeq_gff <- import(DEXSeq_gff_path)



group <- as.character(mcols(DEXSeq_gff)[, "group"])

library(limma)

group_tmp <- strsplit2(group, "gene_id ")

groupg <- group_tmp[, 2]
groupg <- gsub("\"", "", groupg)


group_tmp <- strsplit2(group_tmp[,1], "exonic_part_number")

groupb <- group_tmp[, 2]
groupb <- gsub("\"| |; ", "", groupb)


groupt <- group_tmp[, 1]
groupt <- gsub("transcripts \"|\"; ", "", groupt)


groupt_spl <- strsplit(groupt, "+", fixed = TRUE)

groupgb <- paste0(groupg, ":", groupb)


group_df <- data.frame(gene_id = groupg, transcripts_id = groupt, bin_id = groupgb, stringsAsFactors = FALSE)

group_df <- group_df[group_df$gene_id %in% simulation_details_ds$gene_id, ]

group_spl <- split(group_df, group_df$gene_id)


simulation_details_ds_spl <- split(simulation_details_ds, simulation_details_ds$gene_id)


simu_infoList <- lapply(names(simulation_details_ds_spl), function(g){
  # g = 1
  bins <- group_spl[[g]]
  sim <- simulation_details_ds_spl[[g]]
  

  exon1 <- bins[grepl(sim$transcript_id[1], bins$transcripts_id) & !grepl(sim$transcript_id[2], bins$transcripts_id), "bin_id"]
  if(length(exon1) == 0)
    exon1 <- NA
  
  sim1 <- data.frame(sim[rep(1, length(exon1)), ,drop = FALSE], exon_id = exon1)
    

  exon2 <- bins[grepl(sim$transcript_id[2], bins$transcripts_id) & !grepl(sim$transcript_id[1], bins$transcripts_id), "bin_id"]
  if(length(exon2) == 0)
    exon2 <- NA
  
  sim2 <- data.frame(sim[rep(2, length(exon2)),,drop = FALSE], exon_id = exon2)
  
  return(rbind(sim1, sim2))
  
})



table(sapply(simu_infoList, nrow))



simu_info <- do.call(rbind, simu_infoList)


write.table(simu_info, "3_truth/simu_info.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
































