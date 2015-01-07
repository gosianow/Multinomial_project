
# BioC 14

# Created 19 Oct 2014
# Modyfied 19 Oct 2014


##############################################################################################################

# run DM SQTL pipeline on GEUVADIS data

##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

# load("DM_genotypes/dgeSQTL.RData")


library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v5.R"))

out.dir <- "DM_sQTLanalysis/"
dir.create(out.dir, showWarnings=F, recursive=T)



time_ECD <- system.time(dgeSQTL <- dmSQTLEstimateCommonDisp(dgeSQTL, adjust = TRUE, subset=1e+04, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=10, verbose=FALSE))


time_ECD["elapsed"]/60 ### takes ~ 88 min


time_F <- system.time(dgeSQTL <- dmSQTLFit(dgeSQTL, model="full", dispersion=NULL, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 10))


time_F["elapsed"]/60 ### takes ~ 139 min


save(dgeSQTL, file=paste0(out.dir, "/dgeSQTL.RData"))


load(paste0(out.dir, "/dgeSQTL.RData"))

time_T <- system.time(dgeSQTL <- dmSQTLTest(dgeSQTL, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 10))


time_T["elapsed"]/60/60 ### takes 8.49 hours


name1 <- "dgeSQTL"
write.table(dgeSQTL$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)

write.table(dgeSQTL$table, paste0(out.dir, "/",name1,"_results.txt"), quote=F, sep="\t", row.names=F, col.names=T)
















