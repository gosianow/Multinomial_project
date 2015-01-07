##############################################################################

# BioC 14
# Created 04 Sep 2014:

# For ECCB 

# DM package version 4
# Check if the LR null distribution is ChiSquared
# Simulate null two group data from DM

# Updated 06 Sep 2014:


##############################################################################

setwd("/home/gosia/Multinomial_project/DM_package")

library(edgeR)
library(parallel)
library(dirmult)

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version4/"
source(paste0(Rdir, "dmFunctions_v4.R"))

out.dir <- "PLOTS4/Test_LR_distribution_ECCB/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)


############################################################################
# function to simulate data from DM
############################################################################



simulate_from_DM <- function(sample.size = 3, s = "1" , pi.org = c(0.4, 0.6), g0.org = 30, nr.genes = 100, nM = 200, tot="nbinom", nD=3, out.dir, mc.cores=10, save=TRUE){
  
  dir.create(out.dir, recursive = T, showWarnings = FALSE)
  name1 <- paste0("SIM",s,"_ss",sample.size,"_g0", g0.org, "_")
  # simParams <- as.list(match.call())
  # simParams <- c(mget(names(formals()), envir = as.environment(-1)), sapply(as.list(substitute({ ... })[-1]), deparse))
  simParams <- mget(names(formals()), envir = as.environment(-1))
  sink(file=paste0(out.dir, "/", name1, "simParams",".txt"), split=TRUE)
  print(simParams)
  sink(file=NULL)
  
  g.dir.org <- pi.org * g0.org # gamma
  
  sim <- mclapply(1:nr.genes, function(i){
    # i=1
    # simulate dirichlets
    g.dir <- rdirichlet( sample.size, g.dir.org )
    
    if(tot=="nbinom")
      t <- rnbinom(sample.size, mu=nM, size=nD)
    if(tot=="norm")
      t <- round(rnorm(sample.size, mean=nM, sd=nD))  
    if(tot=="uni")
      t <- rep(nM, sample.size)
    
    # simulate multinomial
    d <- sapply(1:sample.size, function(u) rmultinom(1, prob=g.dir[u,], size=t[u]))    
    genes <- paste0("g", rep(i, length(g.dir.org)), ":e", 1:length(g.dir.org))   
    rownames(d) <- genes
    
    return(d)
    
  }, mc.cores=mc.cores)
  
  sim <- do.call(rbind, sim)
  ids <- strsplit2(rownames(sim), ":", fixed = TRUE)
  dge <- DGEList( counts=sim, group = rep(1, sample.size), genes=data.frame(gene.id=ids[,1], ete.id = rownames(sim)) )
  
  ## add 1
  # dge$counts <- dge$counts + 1 
  # dge$counts[ dge$counts == 0 ] <- 1
  
  
  if(save){
    dir.create(out.dir, recursive = T, showWarnings = F)
    save(dge, name1, g.dir.org, file=paste0(out.dir, "/", name1, "dge",".RData"))   
    save(simParams, file = paste0(out.dir, "/", name1, "simParams",".RData"))
  }
  
  
  return(list(dge=dge, name1=name1, g.dir.org=g.dir.org))
  
  
}


############################################################################
# Simulate data from two group null distribution
############################################################################

###################### Simulations - parameters
## from true gene
# load("/home/gosia/Multinomial_project/Simulations_drosophila_V2/DM_v4/fc/fc_g0_s4_keep0s_subsetInf_DM4adj_dgeDM.RData")
# y <- split(as.data.frame(dgeDM$counts), dgeDM$genes$gene.id)
# gene <- "FBgn0000014"
# y[[gene]]
# pi.org <- dgeDM$fit.null$fit[[gene]]$piH
# g0.org <- round(dgeDM$fit.null$fit[[gene]]$gamma0)
# nM <- mean(colSums(y[[gene]]))
# nD <- nM


pi.org <- c(0.5, 0.3, 0.2)
g0.org <- 100
nM <- 200
nD <- nM


simPar <- list(sample.size = 6, s = paste0(1) , pi.org = pi.org , g0.org = g0.org, nr.genes = 1000, nM = nM, tot = "nbinom", nD = nD, out.dir = out.dir)


###################### Simulate data


sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, save = FALSE)


sim$dge$samples$group <- as.factor(c(rep("C1", simPar$sample.size/2), rep("C2", simPar$sample.size/2)))

dge <- sim$dge
name1 <- sim$name1
g.dir.org <- sim$g.dir.org

save(sim, file=paste0(out.dir, "/", name1, "simNull",".Rdata"))  


############################################################################
# run DM version4 pipeline
############################################################################


## run DM pipeline
# dgeDM <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=20, verbose=FALSE)


# ### take original dispersion
# name1 <- paste0(sim$name1, "dispOrg_")
# dgeDM <- dge
# dgeDM$commonDispersion <- g0.org


### take original dispersion
name1 <- paste0(sim$name1, "2dispOrg_")
dgeDM <- dge
dgeDM$commonDispersion <- 2*g0.org



write.table(dgeDM$commonDispersion, paste0(out.dir, "/",name1,"_commonDispersion.txt"), quote=F, sep="\t", row.names=F, col.names=F)

dgeDM <- dmFit(dgeDM, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

dgeDM <- dmTest(dgeDM, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20)

write.table(dgeDM$table, paste0(out.dir, "/",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dgeDM, file=paste0(out.dir, "/",name1,"_dgeDM.RData"))




############################################################################
# plot LR 
############################################################################

table <- dgeDM$table

pdf(paste0(out.dir, "/", name1, "hist_pvalues.pdf"))

  hist(table[, "PValue"], col="darkorange", breaks=50, cex.main=2, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main="Histogram of p-values")
  
dev.off()


pdf(paste0(out.dir, "/", name1, "hist_LR.pdf"))
LR <- table[, "LR"]
hist(LR, col="dodgerblue", breaks=50, cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Likelihood ratio statistic", main=paste0("Histogram of LR for null model"), freq=FALSE)

x <- seq(0, max(LR), by=0.1)
yChiSq <- dchisq(x, df = length(g.dir.org)-1)
lines(x, yChiSq, col="deeppink", lwd=8)

dev.off()


pdf(paste0(out.dir, "/", name1, "qqplot_LR.pdf"))

#qqplot(LR, rchisq(n = 10000, df = length(g.dir.org)-1), pch=20, cex=2, main="QQ plot")
qqplot(qchisq(p = x/max(x), df = length(g.dir.org)-1), LR, pch=20, cex=3, main="QQ plot", cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Chi-square")
abline(a = 0, b = 1, col="red", lwd=4)

dev.off()























