

setwd("/home/gosia/Multinomial_project/DM_package_devel/")
out.dir <- "DM_output/"
dir.create(out.dir, showWarnings=F, recursive=T)


### load DM
library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")


### Source all R files in DM package
Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])
