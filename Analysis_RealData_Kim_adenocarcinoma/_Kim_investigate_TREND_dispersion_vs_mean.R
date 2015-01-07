#######################################################
# 
# Created 19 Nov 2014 

# Update 19 Nov 2014:

# plots of TREND dispersion vs. mean 


#######################################################
# BioC 2.14

setwd("/home/Shared/data/seq/Kim_adenocarcinoma/")


library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")

source(paste0(Rdir, "dmFunctions_v5.R"))



out.dir <- "PLOTS_DM_v5_TREND_dispVSmean/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)



#######################################################
### calculate the mean gene expression
#######################################################

# model <- "diff_out_Null_normal1"
model <- "diff_out_Null_tumor1"

dir.create(paste0(out.dir,"/", model))

load(paste0("DM_v5/fc/", model ,"/fc_g0_s4_keep0s_subsetInf_DM5adj_dgeDM.RData"))


meanExpr <- sapply(dgeDM$counts, function(g){ mean(colSums(g)) } )

meanExpr <- data.frame(gene_id = names(meanExpr), meanExpr = meanExpr)

head(meanExpr)

table <- meanExpr

#######################################################
# plot dispersion vs mean
#######################################################



### load common dispersions
cDisp <- read.table(paste0("DM_v5/fc/",model,"/fc_g0_s4_keep0s_subsetInf_DM5adj_commonDispersion.txt"))


files <- list.files(path = paste0("DM_v5/fc/", model), pattern = "_results.xls" )
files <- files[grepl(pattern = "TG", files)]

TGmethods <- gsub(pattern = "_results.xls", replacement = "" , files)


for( i in 1:length(TGmethods)){
  # i = 1
  
  tDisp <- read.table(paste0("DM_v5/fc/", model ,"/",TGmethods[i],"_tagwiseDispersion.txt"))
  tName <- paste0(TGmethods[i],"_tagwiseDispersion")
  colnames(tDisp) <- c("gene_id", tName)
  
  
  table <- unique(merge(table, tDisp, by = "gene_id", all.x=TRUE))
  
  pdf(paste0(out.dir, "/", model, "/TREMD_mean_vs_gamma-",TGmethods[i],".pdf"))
  
  smoothScatter(log10(table$meanExpr), log10(table[,tName]), xlab="log10 mean gene expression", ylab="log10 gamma +", nrpoints = Inf, colramp=colorRampPalette(c("white", "grey40")), pch = 19, cex=0.6)
  abline(h = log10(cDisp), col = "red")

  dev.off()
  
}



















































