#######################################################
# 
# Created 29 Oct 2014 

# Aim: to find out what is the origin of FP called by DM

# Update 28 Jan 2015:

# plots of TREND dispersion vs. mean 
# plots of entropy 

# plot of expression for FP genes

#######################################################
# BioC 2.14

setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2")


library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")

source(paste0(Rdir, "dmFunctions_v5.R"))



out.dir <- "PLOTS_DM_v5_FP/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)



# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata



#######################################################
# load results
#######################################################


tableOrg <- read.table(paste0("PLOTS_DM_v5/Table_all_results.xls"), header = T, stringsAsFactors = F)
tableNames <- colnames(tableOrg)

library(RColorBrewer)


allMethods <- gsub("adjPValue_", "",tableNames[grep(pattern = "adjPValue_", tableNames)])
allMethods

ramp <- colorRampPalette(brewer.pal(9,"Set1"))(length(allMethods))

colors.org <- ramp
n.colors.org <- allMethods
names(colors.org) <- n.colors.org

colors <- colors.org
n.colors <- n.colors.org

name <- ""



#######################################################
### calculate the mean gene expression
#######################################################

load("DM_v5/fc/fc_g0_s4_keep0s_subsetInf_DM5adj_dgeDM.RData")


meanExpr <- sapply(dgeDM$counts, function(g){ mean(colSums(g)) } )

meanExpr <- data.frame(gene_id = names(meanExpr), meanExpr = meanExpr)

head(meanExpr)

tableOrg <- unique(merge(tableOrg, meanExpr, by = "gene_id", all.x=TRUE))


#######################################################
# plot dispersion vs mean & FP
#######################################################

table <- tableOrg

### load common dispersions
cDisp <- read.table("DM_v5/fc/fc_g0_s4_keep0s_subsetInf_DM5adj_commonDispersion.txt")


TGmethods <- allMethods[grepl("TG", allMethods) & grepl("adj", allMethods)]


for( i in 1:length(TGmethods)){
  # i = 1
  
  tDisp <- read.table(paste0("DM_v5/fc/",TGmethods[i],"_tagwiseDispersion.txt"))
  tName <- paste0(TGmethods[i],"_tagwiseDispersion")
  colnames(tDisp) <- c("gene_id", tName)
  
  
  table <- unique(merge(table, tDisp, by = "gene_id", all.x=TRUE))
  
  
  filtFPc <- table$status == 0 & (table$adjPValue_fc_g0_s4_keep0s_subsetInf_DM5adj < 0.05) == 1
  filtFPc[is.na(filtFPc)] <- FALSE
  h(table[filtFPc,])
  
  filtFPt <- table$status == 0 & (table[, paste0("adjPValue_", TGmethods[i])] < 0.05) == 1
  filtFPt[is.na(filtFPt)] <- FALSE
  h(table[filtFPt,])
  
  
  
  pdf(paste0(out.dir, "/TREMD_mean_vs_gamma-",TGmethods[i],".pdf"))
  
  
  smoothScatter(log10(table$meanExpr), log10(table[,tName]), xlab="log10 mean gene expression", ylab="log10 gamma +", nrpoints = Inf, colramp=colorRampPalette(c("white", "white")), pch = 19, cex=0.6)
  points(log10(table$meanExpr[filtFPc]), log10(table[filtFPc, tName]), pch = 19, cex=1, col= colors["fc_g0_s4_keep0s_subsetInf_DM5adj"])
  points(log10(table$meanExpr[filtFPt]), log10(table[filtFPt, tName]), pch = 17, cex=0.8, col= colors[TGmethods[i]])
  abline(h = log10(cDisp), col = "red")
  
  
  
  smoothScatter(log10(table$meanExpr), log10(table[,tName]), xlab="log10 mean gene expression", ylab="log10 gamma +", nrpoints = Inf, colramp=colorRampPalette(c("white", "grey40")), pch = 19, cex=0.6)
  abline(h = log10(cDisp), col = "red")
  

  dev.off()
  
}


################## plot the correlation between different tagwise methods

pdf(paste0(out.dir, "/tagwiseDispersions.pdf"), 10, 10)

my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1, col = "red", lwd = 2 )
}

pairs(log10(table[,paste0(TGmethods,"_tagwiseDispersion")]), pch=19, upper.panel = my_line , lower.panel = NULL)

dev.off()



#######################################################
# plots of entropy 
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



############################################################################
# plots of exon counts per gene for FP with high significance 
############################################################################


table <- tableOrg

disp.meth <- "adjPValue_fc_g0_s4_keep0s_subsetInf_DM5adj"

### table with most significant  FP on top  
filtFP <- table$status == 0 & table[,disp.meth] < 0.05 & !is.na(table[,disp.meth])
tableFP <- table[filtFP, c("gene_id",disp.meth), drop = FALSE]
tableFP <- tableFP[order(tableFP[, disp.meth], decreasing = FALSE), ]

### load DM results 

load("DM_v5/fc/fc_g0_s4_keep0s_subsetInf_DM5adj_dgeDM.RData")
dgeDM




### plot exon expression in both conditions 
pdf(paste0(out.dir, "ExonExpression_topFP.pdf"), width = 14, height = 7)

for(g in 1:10){
  # g = 1
  gene <- tableFP[g, "gene_id"]
  expr <- dgeDM$counts[[gene]]
  colnames(expr) <- metadata$SampleName
  rownames(expr) <- subset(dgeDM$genes, gene_id==gene)$ete_id
  
  labels <- strsplit2(rownames(expr), ":")[,2]
  n <- nrow(expr)
  
  
  par(mar=c(8, 5, 4, 2) +  0.1)
  plot(rep(1, n), type="n", ylim=c(0, max(expr)), xaxt="n", ylab="",xlab="", main=gene,  cex.main=2, cex.axis=1.7)
  title(ylab = "Counts", mgp = c(3, 1, 0), cex.lab=2)
  title(xlab = "Flattened exon IDs", mgp = c(6, 1, 0), cex.lab=2)
  axis(side=1, at=1:n, labels=labels, las=2, cex.axis=2)
  
  for(i in 1:3)
    lines(expr[,i], col="firebrick1", type="l", lwd=5)
  for(i in 4:6)
    lines(expr[,i], col="dodgerblue", type="l", lwd=5)

  
  legend("topright", c("Condition1", "Condition2"), lty=1, lwd=5, col=c("firebrick1", "dodgerblue"), cex=2, bty = "n")
  
  
}

dev.off()




### plot exon proportions in both conditions 
pdf(paste0(out.dir, "ExonProportions_topFP.pdf"), width = 14, height = 7)

for(g in 1:10){
  # g = 1
  gene <- tableFP[g, "gene_id"]
  expr <- dgeDM$counts[[gene]]
  colnames(expr) <- metadata$SampleName
  rownames(expr) <- subset(dgeDM$genes, gene_id==gene)$ete_id    
  tot <- colSums(expr)
  prop.smp <- data.frame(t(apply(expr, 1, function(t){ t / tot })))  
  labels <- strsplit2(rownames(expr), ":")[,2]
  n <- nrow(expr)  
  
  
  par(mar=c(8, 5, 4, 2) +  0.1)
  plot(rep(1, n), type="n", ylim=c(0, max(prop.smp)), xaxt="n", ylab="",xlab="", main=gene,  cex.main=2, cex.axis=1.7)
  title(ylab = "Proportions", mgp = c(3, 1, 0), cex.lab=2)
  title(xlab = "Flattened exon IDs", mgp = c(6, 1, 0), cex.lab=2)
  axis(side=1, at=1:n, labels=labels, las=2, cex.axis=2)
  
  for(i in 1:3)
    lines(prop.smp[,i], col="firebrick", type="l", lwd=5)
  for(i in 4:6)
    lines(prop.smp[,i], col="dodgerblue4", type="l", lwd=5)
  
  prop.est <- dgeDM$fit[[gene]]$piH
  
  points(1:n, prop.est[,1], pch = 19, cex = 3, col = "firebrick1")
  points(1:n, prop.est[,2], pch = 18, cex = 3, col = "dodgerblue")
  
  
  legend("topright", c("Condition1", "Condition2"), lty=1, lwd=5, col=c("firebrick1", "dodgerblue"), cex=2, bty = "n")
  
  
}

dev.off()


library(ggplot2)
library(reshape2)
library(gridExtra)


### plot exon proportions in both conditions 
pdf(paste0(out.dir, "ExonProportions_topFP_ggplot.pdf"), width = 7, height = 3)

for(g in 1:20){
  # g = 1
  gene <- tableFP[g, "gene_id"]
  expr <- dgeDM$counts[[gene]]
  colnames(expr) <- metadata$SampleName
  rownames(expr) <- subset(dgeDM$genes, gene_id==gene)$ete_id    
  tot <- colSums(expr)
  labels <- strsplit2(rownames(expr), ":")[,2]
  prop.smp <- data.frame( ete_id =  labels, t(apply(expr, 1, function(t){ t / tot })))  
  n <- nrow(expr)  
  prop.est <- data.frame(ete_id = labels, dgeDM$fit[[gene]]$piH)
  
  prop.smp.m <- melt(prop.smp, id.vars = "ete_id", variable.name = "Samples", value.name = "Proportions")  
  Condition <- prop.smp.m$Samples 
  levels(Condition) <- substr(levels(Condition), 1,2)
  
  prop.est.m <- melt(prop.est, id.vars = "ete_id", variable.name = "Samples", value.name = "Proportions")

  ggp <- ggplot() +
    theme_bw() +
    geom_line(data = prop.smp.m, aes(x = ete_id, y = Proportions, group = factor(Samples), colour = Condition )) +
#     scale_color_manual(values=gg_color_variants(length(var.counts))) 
#     theme(axis.text.x  = element_text(angle=80, vjust=0.5, size=12, colour = gg_color_hue(nlevels(m.prop$Transcript)), face="bold"), panel.background = element_blank(),  axis.line = element_line(colour = "grey")) 

  geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, group = factor(Samples), colour = Samples ), size = 3) +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.5 )) +
  ggtitle(paste0(gene))

  print(ggp)

}

dev.off()




































