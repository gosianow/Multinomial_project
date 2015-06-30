#######################################################
# 
# Created 29 Oct 2014 

# Aim: to find out what is the origin of FP called by DM

# Update 03 Nov 2014:

# plots of TREND dispersion vs. mean 
# plots of entropy 

#######################################################
# BioC 2.14

setwd("/home/Shared/data/seq/brooks_pasilla")


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

load("DM_v5/fc/fc_g0_s4_keep0s_subsetInf_DM5adj_dgeDM.RData")


meanExpr <- sapply(dgeDM$counts, function(g){ mean(colSums(g)) } )

meanExpr <- data.frame(gene_id = names(meanExpr), meanExpr = meanExpr)

head(meanExpr)

table <- meanExpr

#######################################################
# plot dispersion vs mean
#######################################################


### load common dispersions
cDisp <- read.table("DM_v5/fc/fc_g0_s4_keep0s_subsetInf_DM5adj_commonDispersion.txt")


files <- list.files(path = paste0("DM_v5/fc/"), pattern = "_results.xls" )
files <- files[grepl(pattern = "TG", files)]

TGmethods <- gsub(pattern = "_results.xls", replacement = "" , files)


for( i in 1:length(TGmethods)){
  # i = 1
  
  tDisp <- read.table(paste0("DM_v5/fc/",TGmethods[i],"_tagwiseDispersion.txt"))
  tName <- paste0(TGmethods[i],"_tagwiseDispersion")
  colnames(tDisp) <- c("gene_id", tName)
  
  
  table <- unique(merge(table, tDisp, by = "gene_id", all.x=TRUE))
  
  pdf(paste0(out.dir, "/TREMD_mean_vs_gamma-",TGmethods[i],".pdf"))
  
  smoothScatter(log10(table$meanExpr), log10(table[,tName]), xlab="log10 mean gene expression", ylab="log10 gamma +", nrpoints = Inf, colramp=colorRampPalette(c("white", "grey40")), pch = 19, cex=0.6)
  abline(h = log10(cDisp), col = "red")

  dev.off()
  
}




pdf(paste0(out.dir, "/tagwiseDispersions.pdf"), 10, 10)

my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1, col = "red", lwd = 2 )
}

pairs(log10(table[,paste0(TGmethods,"_tagwiseDispersion")]), pch=19, upper.panel = my_line , lower.panel = NULL)

dev.off()






















#######################################################
# plots of entropy //  not sure if make sense bcs this is an estimate of entropy and max entropy is different for different nr of bins 
#######################################################


### density of Dirichlet-multinomial
dDmultinom <- function (x, gamma, log = FALSE){  
  
  size <- sum(x)
  gammaP <- sum(gamma)
  
  r <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(gammaP) - lgamma(size + gammaP) + sum(lgamma(x + gamma) - lgamma(gamma)) 
  
  if (log)
    r
  else exp(r)
}



entropy <- function(Y, gamma){  ##### something must be wrong with this calculations because Entropy H > 0 
  
  lp <- rep(0, ncol(Y))
  
  for(i in 1:ncol(Y)){
    
    lp[i] <- dDmultinom(x=Y[,i], gamma=gamma, log = TRUE)
    
  }
  
  H <- - sum(exp(lp) * lp)
  
  return(H)
  
}




entropyNull <- sapply(1:length(dgeDM$counts), function(g){  
  # g <- 1
  
  if(is.null(dgeDM$fit.null[[g]]))
    return(NA)
  
  Y <- dgeDM$counts[[g]]
  gamma <- dgeDM$fit.null[[g]]$piH * dgeDM$fit.null[[g]]$gamma0

  return(entropy(Y, gamma))
  
} )


entropyNull <- data.frame(gene_id = names(dgeDM$counts), entropyNull = entropyNull)

head(entropyNull)

table <- unique(merge(table, entropyNull, by = "gene_id", all.x=TRUE))

library(beanplot)


pdf(paste0(out.dir, "/TREMD_entropy.pdf"), width = 14, height = 7)

beanplot(entropyNull~num, table[table$num < 15, ], what=c(0,1,1,0), xlab="df", ylab="log Entropy")
boxplot(entropyNull~num, table[table$num < 15, ], xlab="df", ylab="Entropy")
boxplot(entropyNull~num, table[table$num < 15, ], log = "y", xlab="df", ylab="log Entropy")

dev.off()



TGmethods <- allMethods[grepl("TG", allMethods) & grepl("adj", allMethods)]


for( i in 1:length(TGmethods)){
  # i = 1
  
 tName <- paste0(TGmethods[i],"_tagwiseDispersion")
 
 pdf(paste0(out.dir, "/TREMD_entropy_vs_dispersion-",TGmethods[i],".pdf"))
 
 smoothScatter(table$entropyNull, log10(table[,tName]), xlab="entropy", ylab="log10 gamma +", nrpoints = Inf, colramp=colorRampPalette(c("white", "grey40")), pch = 19, cex=0.6)

 dev.off()
 

}



#######################################################
# plots of exon counts per gene for some interesting genes / TODO
#######################################################





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










