# Created 27 Nov 2014
# BioC 2.14



#### metadata


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2")

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata

#############################
# load packages
#############################

library(edgeR)
library(parallel)
library(dirmult)
library(limma)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v5.R"))




load("DM_v5/fc/fc_g0_s4_keep0s_subsetInf_DM5adj_dgeDM.RData")


#############################
# Exon/bin expression
#############################












pdf("PLOTS/ExploreDM/UniqTruePositives2.pdf", width = 14, height = 7)

for(g in uniq.tp[1]){
  
  dge.g.keep <- dgec[dgec$gene.id==g,,drop=FALSE]
  dge.g <- dgecnf[dgecnf$gene.id==g,,drop=FALSE]
  ds.ex.g <- ds.ex[ds.ex$Gene==g, , drop=FALSE]$exon
  
  col.axis <- ifelse(dge.g$ete.id %in% ds.ex.g, col3, "grey50")
  
  
  labels <- strsplit2(dge.g$ete.id, ":")[,2]
  genes.keep <- dge.g$ete.id %in% dge.g.keep$ete.id
  labels[genes.keep] <- paste0("*" ,labels[genes.keep], "*")
  
  n <- length(dge.g$ete.id)
  
  par(mar=c(8, 5, 4, 2) +  0.1)
  plot(rep(1, n), type="n", ylim=c(0, 10), xaxt="n", ylab="",xlab="", main=g,  cex.main=2, cex.axis=2)
  title(ylab = "log2(counts + 1)", mgp = c(3, 1, 0), cex.lab=1.5)
  title(xlab = "Flattened exon ID", mgp = c(6, 1, 0), cex.lab=1.5)
  
  for(j in 1:length(dge.g$ete.id))
    axis(side=1, at=j, labels=labels[j], las=2, col.axis = col.axis[j], cex.axis=2)
  
  for(i in 1:3)
    lines(log2(dge.g[,2+i]+1), col=col1, type="l", lwd=3)
  for(i in 4:6)
    lines(log2(dge.g[,2+i]+1), col=col2, type="l", lwd=3)
   
  legend("topright", c("Condition1", "Condition2"), lty=1, lwd=3, col=c(col1, col2), cex=2)
  
}
dev.off()

