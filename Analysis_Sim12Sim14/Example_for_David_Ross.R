setwd("/home/gosia/Multinomial_project/Simulations")

library("edgeR")
library("DEXSeq")



load(paste0("PLOTS/", "dge_design.RData"))

## to run the analysis // All my functions are in this file
source("/home/gosia/R/R_Multinomial_project/glmFitDM.R")

# glmfit <- glmFitDM(dge, design, optimization="optim") # or optimization="nlm"
# glmlrt <- glmLRTDM(glmfit, coef=2)

## to load workspace: glmlrt
name <- "optim_"
#name <- "nlm_"

load(paste0("PLOTS/optim/", name, "glmlrt.RData"))

head(glmlrt$table)


# table with simulation info status: 1 spliced, 0 not spliced
table <- read.table(paste0("PLOTS/optim/",name,"table_MISO.xls"), header=T)


pdf("PLOTS/LR_hist.pdf")
hist(table$LR)
dev.off()



#######################################################
# check the statistics
# $loglikh
# $Y
# $bh - estimates of Betas
# $gh - gammas from Dirichlet exp(X %*% t(bh))
# $th - probabilities gh / rowSums(gh)
#######################################################

# I have chosen genes from table

# genes with very big LR
# FP 
gene <- "FBgn0002283"
gene <- "FBgn0040284"


# FN
gene <- "FBgn0033232"
gene <- "FBgn0002626"


# LR < 0

# TN
gene <- "FBgn0025725"
gene <- "FBgn0032987"
# FN
gene <- "FBgn0031228"
gene <- "FBgn0259937"


glmlrt$FitDM[[gene]]

glmlrt$fit.null[[gene]]









