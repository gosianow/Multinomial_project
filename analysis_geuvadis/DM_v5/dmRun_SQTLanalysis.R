
# BioC 2.14

# Created 19 Oct 2014
# Modyfied 19 Oct 2014


##############################################################################################################

# run DM SQTL pipeline on GEUVADIS data

##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

##### dgeSQTL produced with dmSNPsFiltering.R
# load("DMv5_sQTL_analysis/DM_genotypes/dgeSQTL.RData")


library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v5.R"))

out.dir <- "DMv5_sQTL_analysis/"
dir.create(out.dir, showWarnings=F, recursive=T)



time_ECD <- system.time(dgeSQTL <- dmSQTLEstimateCommonDisp(dgeSQTL, adjust = TRUE, subset=1e+04, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=10, verbose=FALSE))

dgeSQTL$commonDispersion ## 3.409244


time_ECD["elapsed"]/60 ### takes ~ 88 min


time_F <- system.time(dgeSQTL <- dmSQTLFit(dgeSQTL, model="full", dispersion=NULL, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 10))


time_F["elapsed"]/60 ### takes ~ 139 min


save(dgeSQTL, file=paste0(out.dir, "/dgeSQTL.RData"))

load(paste0(out.dir, "/dgeSQTL.RData"))


time_T <- system.time(dgeSQTL <- dmSQTLTest(dgeSQTL, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 10))


time_T["elapsed"]/60/60 ### takes 8.49 hours ### WHY SO SLOW?


name1 <- "dgeSQTL"
write.table(dgeSQTL$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)

write.table(dgeSQTL$table, paste0(out.dir, "/",name1,"_results.txt"), quote=F, sep="\t", row.names=F, col.names=T)
















