# BioC 2.14

########################################################
# wrapper_featureCounts_v3
########################################################


setwd("/home/gosia/Multinomial_project/SimulationsV2")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


source("/home/gosia/R/R_Multinomial_project/Analysis_SimV2/wrapper_featureCounts_v3.R")


gtfFile <- "/home/Shared_penticton/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf"
out.dir <- "featureCounts_v3"


PE_bam_files <-  paste0("/home/Shared/tmp/Simulation/simulation_drosophila_V2/sam/",  paste0("sim_samp", metadata$SampleName1, "_s.bam"))
PE_bam_files


fc <- wrapper_featureCounts(gtfFile, PE_bam_files=PE_bam_files,  nthreads=25)











