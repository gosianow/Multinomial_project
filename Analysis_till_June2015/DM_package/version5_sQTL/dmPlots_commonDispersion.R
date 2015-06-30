
##############################################################################

# BioC 14
# Created 10 Oct 2014:

# simulation form DM

# common dispersion plots

##############################################################################


setwd("/home/gosia/Multinomial_project/DM_package")

library(edgeR)
library(parallel)
library(dirmult)

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"
source(paste0(Rdir, "dmFunctions_v5.R"))

out.dir <- "PLOTS5/commonDispersion/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)


############################################################################
# function to simulate data from DM
############################################################################


### updated for new dmFunctions_v5 - dge$counts is a list
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
  dge <- DGEList( counts=sim, group = rep(1, sample.size), genes=data.frame(gene_id=ids[,1], ete_id = rownames(sim)) )
  
  dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels=unique(dge$genes$gene_id)))
  dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR
  
  
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
# Simulate data from null distribution
############################################################################

###################### Simulations - parameters

# load("/home/gosia/Multinomial_project/Simulations_drosophila_V2/DM_v4/fc/fc_g0_s4_keep0s_subsetInf_DM4adj_dgeDM.RData")
# y <- split(as.data.frame(dgeDM$counts), dgeDM$genes$gene.id)
# gene <- "FBgn0262731"
# y[[gene]]
# pi.org <- dgeDM$fit.null$fit[[gene]]$piH
# g0.org <- round(dgeDM$fit.null$fit[[gene]]$gamma0) ## too high 
# nM <- mean(colSums(y[[gene]]))


pi.org <- c(0.5, 0.3, 0.2)
g0.org <- 100
nM <- 200
nD <- nM

### low number of samples 
simPar <- list(sample.size = 3, s = paste0(1) , pi.org = pi.org , g0.org = g0.org, nr.genes = 1000, nM = nM, tot = "nbinom", nD = nD, out.dir = out.dir)

### high number of samples 
# simPar <- list(sample.size = 15, s = paste0(2) , pi.org = pi.org , g0.org = g0.org, nr.genes = 1000, nM = nM, tot = "nbinom", nD = nD, out.dir = out.dir)


###################### Simulate data

R <- 50
# R <- 1 

simList <- list()

for(r in 1:R){
  # r=25
  simList[[r]] <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, save = FALSE)
}

name1 <- simList[[1]]$name1
g.dir.org <- simList[[1]]$g.dir.org
gamma0org <- sum(g.dir.org) 


save(simList, file=paste0(out.dir, "/", name1, "simList",".Rdata"))  


############################################################################
# common dispersion 
############################################################################


commonDipEsts <- data.frame(equalTheta = rep(0, R))

for(r in 1:R){
  # r=1
  cat("r:", r, fill = T)
  dge <- simList[[r]]$dge  
  #### constrOptim2
  commonDipEsts$constrOptim2[r] <- dmEstimateCommonDisp(dge, group=NULL, adjust = FALSE, subset=Inf, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=30, verbose=FALSE)$commonDispersion
}  


for(r in 1:R){
  # r=1
  cat("r:", r, fill = T)
  dge <- simList[[r]]$dge 
  #### constrOptim2ADJ
  commonDipEsts$constrOptim2ADJ[r] <- dmEstimateCommonDisp(dge, group=NULL, adjust = TRUE, subset=Inf, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=30, verbose=FALSE)$commonDispersion
}  


write.table(commonDipEsts, paste0(out.dir, "/", name1, "commonDipEsts",".xls"), quote = F, row.names = F, col.names = T, sep = "\t")  


library(beanplot)

commonDipEsts <- commonDipEsts[,c("constrOptim2", "constrOptim2ADJ")]
gamma0org <- sum(g.dir.org)
SE <- commonDipEsts - gamma0org
colnames(SE) <- c( "CPL", "ACPL")



pdf(paste0(out.dir, "/", name1, "MSEcommonDisp",".pdf"))

beanplot(SE, what=c(0,1,1,0), main=paste0("Standard Error","\n gamma: ", paste(g.dir.org, collapse=", ")), ylim=c(min(SE)-gamma0org^0.5, max(SE)+gamma0org^0.5), cex=3, cex.lab=2, cex.axis=1.5, cex.main=1.5, col = list(c("chartreuse3", "chartreuse3", "chartreuse3", "black"), c("orange", "orange", "orange", "black")) )
abline(h=0, lty=2, lwd=4)

dev.off()


# install.packages("vioplot", lib="/home/gosia/R/libraries/3.1.0/")


library(vioplot)


pdf(paste0(out.dir, "/", name1, "MSEcommonDisp_violin",".pdf"))

vioplot(SE$CPL, SE$ACPL, names=c("CPL", "ACPL"), col= c("chartreuse3", "orange") )
title(paste0("Standard Error","\n gamma: ", paste(g.dir.org, collapse=", ")))

dev.off()






















