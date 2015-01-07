#######################################################
# 
# Created 23 Oct 2014 / Last updated 23 Oct 2014

# 
#######################################################
# BioC 3.0


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


#######################################################
# run Ballgown
####################################################### 

######## tablemaker

merged.gtf <- "cuffdiff/merged/merged.gtf"

cmd <- with(metadata, paste0("/home/gosia/Programs/tablemaker-2.1.1.Linux_x86_64/tablemaker -p 20 -q -W -G ", merged.gtf, " -o ", "tablemaker_2.1.1/samp", SampleName1 , " /home/Shared/tmp/Simulation/simulation_drosophila_V2/bam/sim_samp",SampleName1 ,"_s.bam \n"))


for(i in 1:length(cmd))
  system(cmd[i])




######## ballgown

library(ballgown)

# data_directory <- "tablemaker_2.1.1/"
# 
# bg <- ballgown(dataDir=data_directory, samplePattern="samp", meas=c("FPKM"))

# Error in ballgown(dataDir = data_directory, samplePattern = "samp", meas = c("FPKM")) : Transcript ids were either not the same or not in the same order across samples. Double check t_data.ctab for each sample.

# transcript_fpkm <- texpr(bg, 'FPKM')


sampleNames(bg)

## phenotype information
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0), each=10))



# plotTranscripts(gene='XLOC_000454', gown=bg, samples='sample12', meas='FPKM', colorby='transcript', main='transcripts from gene XLOC_000454: sample 12, FPKM')
# 
# plotMeans('XLOC_000454', bg, groupvar='group', meas='FPKM', colorby='transcript')


### Differential expression analysis

stat_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group')
head(stat_results)














