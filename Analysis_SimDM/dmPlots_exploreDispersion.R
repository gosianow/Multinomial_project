##############################################################################

# BioC 3.0
# Created 27 Nov 2014:
# For IMLS presentation

# Simulate data from DM
# Run DM_0.1.1

# See how the tagwise dispersion estimates depends on sample size and total counts 


##############################################################################


setwd("/home/gosia/Multinomial_project/Simulations_DM/")

out.dir <- "ExploreDispersion/"
dir.create(out.dir, showWarnings=F, recursive=T)


library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)

#### function to simulate data from DM
source("/home/gosia/R/R_Multinomial_project/Analysis_SimDM/simulate_from_DM.R")


##################################################################################
# Simulate data from two group null distribution I) with common dispersion
# Scenario: different total counts : nM
# consider different sample size: sample.size = 3, 10
##################################################################################


### Scenario parameters

scenario.nM <- c(20, 1000)
sname <- "scenario.nM"
s <- paste0(sname, scenario.nM)
out.dir.s <- paste0(out.dir, "/", sname, "/")
dir.create(out.dir.s, showWarnings=F, recursive=T)

nr.genes <- 1000
mcCores <- 10
LRlist <- list()
tableList <- list()


#### Simulate data
simList <- list()

for(sc in 1:length(s)){
  # sc = 1
  cat("*****",s[sc], fill = TRUE)
  
  #### Parameters
  nBins <- 3
  pi.org <- rep(1, nBins)/nBins
  g0.org <- 300
  nM <- scenario.nM[sc]
  nD <- nM
  sample.size <- 3
  
  simPar <- list(sample.size = sample.size*2, s = s[sc] , pi.org = pi.org , g0.org = g0.org, nr.genes = nr.genes, nM = nM, tot = "uni", nD = nD, out.dir = out.dir.s)
  
  sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, mc.cores=mcCores, save = FALSE)
  
  sim$dge$samples$group <- as.factor(c(rep("C1", simPar$sample.size/2), rep("C2", simPar$sample.size/2)))
  
  simList[[s[sc]]] <- sim
  
}

save(simList, file=paste0(out.dir.s, "/simList.RData"))



##################################################################################
#### Run DM common dispersion & tagwise dispersoin - different modeDisp and trend
##################################################################################


### commonDispersion

commonDispersion <- numeric(length(s))
names(commonDispersion) <- s

for(sc in seq(length(s))){
  # sc = 1
  
  sim <- simList[[s[sc]]]
  dge <- sim$dge
  
  dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-00, mcCores=mcCores, verbose=FALSE)
  
  commonDispersion[sc] <- dgeDM$commonDispersion
  
}

write.table(commonDispersion, paste0(out.dir.s, "/cmnDisp.txt"), quote=F, sep="\t", row.names=F, col.names=T)


##### tagwiseDispersion 
mcCores <- 10
tagwiseDispersionList <- list()


# name0 <- c("optimize", "optim", "constrOptim", "grid-none", "grid_commonDispersion", "grid-tagwiseDispersion")
# name <- paste0("tgDisp_", name0)
# 
# modeDisp <- c("optimize", "optim", "constrOptim", "grid", "grid", "grid")
# trend <- c("none", "none", "none", "none", "commonDispersion", "trendedDispersion")

name0 <- c("constrOptim")
name <- paste0("tgDisp_", name0)

modeDisp <- c("constrOptim")
trend <- c("none")


priorDf <- sample.size * 2 - 2

for(i in 1:length(name)){
  cat(name[i], fill = T)
  tagwiseDispersion <- matrix(0, nr.genes, length(s))
  colnames(tagwiseDispersion) <- s
  
  for(sc in seq(length(s))){
    # sc = 1
    cat(s[sc], fill = T)
    
    sim <- simList[[s[sc]]]
    dge <- sim$dge
    dge$commonDispersion <- commonDispersion[sc]
    
    dgeDM <- dmEstimateTagwiseDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2G", epsilon = 1e-05, maxIte = 1000, modeDisp=modeDisp[i], interval = c(0, 1e+5), tol = 1e-00,  initDisp = 10, initWeirMoM = TRUE, gridLength = 15, gridRange = c(-7, 7), trend = trend[i], priorDf = priorDf, span = 0.3, mcCores = mcCores, verbose=FALSE)
    
    tagwiseDispersion[,sc] <- dgeDM$tagwiseDispersion
    
  }
  
  tagwiseDispersionList[[name[i]]] <- tagwiseDispersion
  write.table(tagwiseDispersion, paste0(out.dir.s, "/",name[i],".txt"), quote=F, sep="\t", row.names=F, col.names=T)
  
}


save(tagwiseDispersionList, file=paste0(out.dir.s, "/tgDisp.RData"))


##################################################################################
#### Plots 
# Scenario: different total counts : nM
##################################################################################

load(paste0(out.dir.s, "/tgDisp.RData"))


library(RColorBrewer)
ramp <- colorRampPalette(brewer.pal(12,"Paired"))(length(tagwiseDispersionList))

col <- list()
for(op in 1:length(tagwiseDispersionList))
  col[[op]] <- c(rep(ramp[op], 3), "black")


library(beanplot)


pdf(paste0(out.dir.s, "/tagwiseDispersion_modeDisp.pdf"), width = 14, height = 7)

for(sc in 1:length(s)){
  # sc = 1 
  
  tgD <- lapply(tagwiseDispersionList, function(i) i[,sc])
  tgD <- do.call(cbind, tgD)
  
  beanplot(data.frame(log10(tgD)), what=c(0,1,1,0), main=s[sc], log="", cex=3, cex.lab=1.5, cex.axis=1, cex.main=1.5, col = col, ylab = "log10 tagwise dispersion", names = gsub("tgDisp_", "", colnames(tgD)))
  abline(h = log10(g0.org), lwd=4, col="grey", lty=3)
  
  
  #     boxplot(data.frame(log10(tgD)), cex=3, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, col = ramp[1:length(s)])
  
  
}

dev.off()




pdf(paste0(out.dir.s, "/tagwiseDispersion_constrOptim.pdf"), width = 7, height = 7)

for(sc in 1:length(s)){
  # sc = 1 
  
  tgD <- tagwiseDispersionList[["tgDisp_constrOptim"]][, s[sc]]
  nM <- rep(scenario.nM[sc], length(tgD))
  
  nM <- jitter(nM, factor = 1, amount = NULL)

  plot(nM, log10(tgD), col="grey", lwd=2, cex=1.5, cex.lab=1.45, cex.axis=1.5, xlab="", ylab="log10 gamma +", pch=19, xaxt="n", ylim=c(1, 4))
  abline(h = log10(g0.org), col="deeppink", lwd=5, lty=1)
  
  
  plot(nM, tgD, col="grey", lwd=2, cex=1.5, cex.lab=1.45, cex.axis=1.5, xlab="", ylab="gamma +", pch=19, xaxt="n", ylim=c(0, 1e4))
  abline(h = g0.org, col="deeppink", lwd=5, lty=1)
  
  
}

dev.off()



pdf(paste0(out.dir.s, "/tagwiseDispersion_constrOptim_commonDisp.pdf"), width = 7, height = 7)

for(sc in 1:length(s)){
  # sc = 1 
  
  tgD <- tagwiseDispersionList[["tgDisp_constrOptim"]][, s[sc]]
  nM <- rep(scenario.nM[sc], length(tgD))
  
  nM <- jitter(nM, factor = 1, amount = NULL)
  
  plot(nM, log10(tgD), col="grey", lwd=2, cex=1.5, cex.lab=1.45, cex.axis=1.5, xlab="", ylab="log10 gamma +", pch=19, xaxt="n", ylim=c(1, 4))
  abline(h = log10(g0.org), col="deeppink", lwd=5, lty=1)
  abline(h = log10(commonDispersion[s[sc]]), col="orange", lwd=10, lty=2)
  
  
  plot(nM, tgD, col="grey", lwd=2, cex=1.5, cex.lab=1.45, cex.axis=1.5, xlab="", ylab="gamma +", pch=19, xaxt="n", ylim=c(0, 1e4))
  abline(h = g0.org, col="deeppink", lwd=5, lty=1)
  abline(h = commonDispersion[s[sc]], col="orange", lwd=10, lty=2)
  
  
}

dev.off()




















