##############################################################################

# BioC 14
# Created 28 Aug 2014

# Updated 29 Aug 2014:

# Check for errors in optimization

##############################################################################

setwd("/home/gosia/Multinomial_project/DM_package")
Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/"

library(edgeR)
library(parallel)
library(dirmult)

out.dir <- "PLOTS/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)

source("/home/gosia/R/R_Multinomial_project/DM_adjusted/dirmult_code.R")
source(paste0(Rdir, "dmFunctions.R"))


######################################
# function to simulate data from DM
######################################



simulate_from_DM <- function(sample.size = 3, s = "1", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 30, nr.genes = 100, nM = 200, tot="uni", nD=3, out.dir, mc.cores=20, save=TRUE){
  
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


######################################
# MSE of common dispersion estimates
######################################

source("/home/gosia/R/R_Multinomial_project/DM_adjusted/dirmult_code.R")
source(paste0(Rdir, "dmFunctions.R"))

out.dir <- "PLOTS/Test_for_errors/"

###################### Simulations - parameters

simPar <- list(sample.size = 10, s = "1b", pi.org = c(0.9, 0.1), g0.org = 30, nr.genes = 100, nM = 200, tot = "nbinom", nD = 3, out.dir = out.dir)


###################### Simulate data

R <- 50

simList <- list()

for(r in 1:R){
  # r=25
  simList[[r]] <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, save = FALSE)
}

name1 <- simList[[1]]$name1
g.dir.org <- simList[[1]]$g.dir.org

save(simList, file=paste0(out.dir, "/", name1, "simList",".Rdata"))  

load(paste0(out.dir, "/","SIM1b_ss10_g030_simList.Rdata"))

###################### Estimation of common dispersion

# source(paste0(Rdir, "dmFunctions.R"))

commonDipEsts <- data.frame(equalTheta = rep(0, R))


for(r in 1:R){
  # r=25
  cat("r:", r, fill = T)  
  dge <- simList[[r]]$dge  
  #### equalTheta
  theta <- round(1 / (1 + simPar$g0.org), 2)
  yt <- lapply(split(as.data.frame(dge$counts), dge$genes$gene.id, drop=TRUE), t) # t - transpose
  et <- equalTheta(yt, theta=theta, trace = FALSE)
  commonDipEsts$equalTheta[r] <- sum(et$gamma[[1]])
}  

for(r in 1:R){
  # r=25
  cat("r:", r, fill = T)
  dge <- simList[[r]]$dge
  y <- split(as.data.frame(dge$counts), dge$genes$gene.id)   
  #### constrOptim
  commonDipEsts$constrOptim[r] <- dmEstimateCommonDisp(y, adjust = FALSE, mode = "constrOptim", mcCores=10, interval = c(0, 1e+3), tol = 1e-05)
}  

for(r in 1:R){
  # r=25
  cat("r:", r, fill = T)
  dge <- simList[[r]]$dge
  y <- split(as.data.frame(dge$counts), dge$genes$gene.id)   
  #### constrOptim2
  commonDipEsts$constrOptim2[r] <- dmEstimateCommonDisp(y, adjust = FALSE, mode = "constrOptim2", mcCores=10, interval = c(0, 1e+3), tol = 1e-05)
}  

for(r in 1:R){
  # r=25
  cat("r:", r, fill = T)
  dge <- simList[[r]]$dge
  y <- split(as.data.frame(dge$counts), dge$genes$gene.id)  
  #### FisherObs
  commonDipEsts$FisherObs[r] <- dmEstimateCommonDisp(y, adjust = FALSE, mode = "obs", mcCores=10, interval = c(0, 1e+3), tol = 1e-05, verbose=FALSE)
}  

for(r in 1:R){
  # r=25
  cat("r:", r, fill = T)
  dge <- simList[[r]]$dge
  y <- split(as.data.frame(dge$counts), dge$genes$gene.id)   
  #### FisherExp
  commonDipEsts$FisherExp[r] <- dmEstimateCommonDisp(y, adjust = FALSE, mode = "exp", mcCores=10, interval = c(0, 1e+3), tol = 1e-05, verbose=TRUE)  
}

for(r in 1:R){
  # r=1
  cat("r:", r, fill = T)
  dge <- simList[[r]]$dge
  y <- split(as.data.frame(dge$counts), dge$genes$gene.id)   
  #### optim2
  commonDipEsts$optim2[r] <- dmEstimateCommonDisp(y, adjust = FALSE, mode = "optim2", mcCores=10, interval = c(0, 1e+3), tol = 1e-05, verbose=TRUE)  
}

for(r in 1:R){
  # r=1
  cat("r:", r, fill = T)
  dge <- simList[[r]]$dge
  y <- split(as.data.frame(dge$counts), dge$genes$gene.id)   
  #### optim2NM
  commonDipEsts$optim2NM[r] <- dmEstimateCommonDisp(y, adjust = FALSE, mode = "optim2NM", mcCores=10, interval = c(0, 1e+3), tol = 1e-05, verbose=FALSE)  
}


write.table(commonDipEsts, paste0(out.dir, "/", name1, "commonDipEsts2",".xls"), quote = F, row.names = F, col.names = T, sep = "\t")  



###################### Plots of SE
library(beanplot)

gamma0org <- sum(g.dir.org)
SE <- commonDipEsts - gamma0org

pdf(paste0(out.dir, "/", name1, "MSEcommonDisp",".pdf"))

beanplot(commonDipEsts, col="orangered", what=c(0,1,1,0), main=paste0("gamma0H","\n gamma: ", paste(g.dir.org, collapse=", ")), ylim=c(min(commonDipEsts)-gamma0org^0.25, max(commonDipEsts)+gamma0org^0.25), log="y")
abline(h=sum(g.dir.org), lty=2)

beanplot(SE, col="deepskyblue2", what=c(0,1,1,0), main=paste0("SE: gamma0H - gamma0","\n gamma: ", paste(g.dir.org, collapse=", ")), ylim=c(min(SE)-gamma0org^0.25, max(SE)+gamma0org^0.25))
abline(h=0, lty=2)

dev.off()
































