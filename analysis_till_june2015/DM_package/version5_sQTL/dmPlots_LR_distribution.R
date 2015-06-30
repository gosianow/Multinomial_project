##############################################################################

# BioC 14
# Created 10 Nov 2014:

# Simulate null two group data from DM
# Run DM package version 5
# Check if the LR null distribution is ChiSquared for different scenarios: dispersion, count, sample size, number of bins, pi distribution

# + qqplots 

# Updated 10 Nov 2014:


##############################################################################

setwd("/home/gosia/Multinomial_project/DM_package/PLOTS5/")

library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v5.R"))

out.dir <- "Test_LR_distribution/"
dir.create(out.dir, showWarnings=F, recursive=T)


############################################################################
# function to simulate data from DM
############################################################################



simulate_from_DM <- function(sample.size = 5, s = "1" , pi.org = c(1/3, 1/3, 1/3), g0.org = 100, nr.genes = 100, nM = 150, tot="nbinom", nD=3, out.dir, mc.cores=10, save=TRUE){

  
  dir.create(out.dir, recursive = T, showWarnings = FALSE)
  name1 <- paste0("SIM", s , "_")
  # simParams <- as.list(match.call())
  # simParams <- c(mget(names(formals()), envir = as.environment(-1)), sapply(as.list(substitute({ ... })[-1]), deparse))
  simParams <- mget(names(formals()), envir = as.environment(-1))
  sink(file=paste0(out.dir, "/", name1, "simParams",".txt"), split=TRUE)
  print(simParams)
  sink(file=NULL)
  
  if(length(g0.org == 1))
    g0.org <- rep(g0.org, nr.genes)

  
  sim <- mclapply(1:nr.genes, function(i){
    # i=1
    
    g.dir.org <- pi.org * g0.org[i] # gamma
    
    # simulate dirichlets
    g.dir <- rdirichlet( sample.size, g.dir.org)
    
    # simulate total counts
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
  
  
  g.dir.org <- pi.org * g0.org[1] ### Normally do not need this param to return but I keep it by now

  if(save){
    dir.create(out.dir, recursive = T, showWarnings = F)
    save(dge, name1, g.dir.org, file=paste0(out.dir, "/", name1, "dge",".RData"))   
    save(simParams, file = paste0(out.dir, "/", name1, "simParams",".RData"))
  }
  
  
  return(list(dge=dge, name1=name1, g.dir.org=g.dir.org))
  
  
}



############################################################################

# Simulate data from two group null distribution + run DM version4 pipeline
# + perform KS test 

############################################################################



#########################################
# Scenario: different original dispersion : g0.org
#########################################


### Scenario parameters

scenario.g0.org <- c(10, 50, 100, 500, 1000)
sname <- "scenario.g0.org"
s <- paste0(sname, scenario.g0.org)
out.dir <- "Test_LR_distribution/scenario.g0.org/"
dir.create(out.dir, showWarnings=F, recursive=T)

R <- 1
nr.genes <- 1e+5
mcCores <- 10
LRlist <- list()
tableList <- list()

# Simulate data + calculate LR statistics

for(sc in 1:length(s)){
  # sc = 1
  cat("*****",s[sc], fill = TRUE)
  
  #### Parameters
  nBins = 3
  pi.org <- rep(1, nBins)/nBins
  g0.org <- scenario.g0.org[sc]
  nM <- 150
  nD <- nM
  sample.size <- 5
  
  simPar <- list(sample.size = sample.size*2, s = s[sc] , pi.org = pi.org , g0.org = g0.org, nr.genes = nr.genes, nM = nM, tot = "uni", nD = nD, out.dir = out.dir)
  
  
  #### Simulate data + run DM
  
  LR <- matrix(0, nrow = simPar$nr.genes, ncol = R)
  
  for(r in 1:R){
    # r=1
    cat("***** R", r, fill=TRUE)
    
    sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, mc.cores=mcCores, save = FALSE)
    
    sim$dge$samples$group <- as.factor(c(rep("C1", simPar$sample.size/2), rep("C2", simPar$sample.size/2)))
    

    name1 <- sim$name1
    g.dir.org <- sim$g.dir.org
    gamma0org <- sum(g.dir.org) 
    
    
    dgeDM <- sim$dge
    dgeDM$commonDispersion <- g0.org
    
    
    dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
    
    dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
    
    tableList[[s[sc]]] <- dgeDM$table
    
    LR[, r] <- dgeDM$table$LR
    
  }
  
  LRlist[[s[sc]]] <- LR
  
  
}


save(tableList, file = paste0(out.dir, "/", "tableList.RData"))




pdf(paste0(out.dir, "/qqplot_LR.pdf"))

for(sc in 1:length(s)){
  
  LR <- LRlist[[s[sc]]]
  
     qqplot(rchisq(n = 1e+5, df = length(pi.org)-1), LR, pch=20, cex=3, main= paste("QQ plot \n ", s[sc]), cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Chi-square")
     abline(a = 0, b = 1, col="red", lwd=4)
           
    }

dev.off()




for(sc in 1:length(s)){
  
  LR <- LRlist[[s[sc]]]
  
  png(paste0(out.dir, "/qqplot_LR", s[sc],".png"))
  
  qqplot(rchisq(n = 1e+5, df = length(pi.org)-1), LR, pch=20, cex=3, main= paste("QQ plot"), cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Chi-square")
  abline(a = 0, b = 1, col="red", lwd=4)
  
  dev.off()
  
}


tb <- tableList[[1]]

plot(tb$PValue, tb$FDR)
dev.off()

sum(tb$FDR < 0.05)

#########################################
# Scenario: different sample size : sample.size
#########################################


### Scenario parameters

scenario.sample.size <- c(3, 5, 10, 20)
sname <- "scenario.sample.size"
s <- paste0(sname, scenario.sample.size)
out.dir <- "Test_LR_distribution/scenario.sample.size/"
dir.create(out.dir, showWarnings=F, recursive=T)

R <- 1
nr.genes <- 1e+5
mcCores <- 10
LRlist <- list()
tableList <- list()

# Simulate data + calculate LR statistics

for(sc in 1:length(s)){
  # sc = 1
  cat("*****",s[sc], fill = TRUE)
  
  #### Parameters
  nBins <- 3
  pi.org <- rep(1, nBins)/nBins
  g0.org <- 100
  nM <- 150
  nD <- nM
  sample.size <- scenario.sample.size[sc]
  
  simPar <- list(sample.size = sample.size*2, s = s[sc] , pi.org = pi.org , g0.org = g0.org, nr.genes = nr.genes, nM = nM, tot = "uni", nD = nD, out.dir = out.dir)
  
  
  #### Simulate data + run DM
  
  LR <- matrix(0, nrow = simPar$nr.genes, ncol = R)
  
  for(r in 1:R){
    # r=1
    cat("***** R", r, fill=TRUE)
    
    sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, mc.cores = mcCores, save = FALSE)
    
    sim$dge$samples$group <- as.factor(c(rep("C1", simPar$sample.size/2), rep("C2", simPar$sample.size/2)))
    
    
    name1 <- sim$name1
    g.dir.org <- sim$g.dir.org
    gamma0org <- sum(g.dir.org) 
    
    
    dgeDM <- sim$dge
    dgeDM$commonDispersion <- g0.org
    
    
    dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
    
    dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
        
    tableList[[s[sc]]] <- dgeDM$table
    
    LR[, r] <- dgeDM$table$LR
    
  }
  
  LRlist[[s[sc]]] <- LR
  
  
}


save(tableList, file = paste0(out.dir, "/", "tableList.RData"))



pdf(paste0(out.dir, "/qqplot_LR.pdf"))

for(sc in 1:length(s)){
  
  LR <- LRlist[[s[sc]]]
  
  qqplot(rchisq(n = 1e+5, df = length(pi.org)-1), LR, pch=20, cex=3, main= paste("QQ plot \n ", s[sc]), cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Chi-square")
  abline(a = 0, b = 1, col="red", lwd=4)
  
}

dev.off()



#########################################
# Scenario: different number of bins: nBins
#########################################


### Scenario parameters

scenario.nBins <- c(3, 5, 10, 20)
sname <- "scenario.nBins"
s <- paste0(sname, scenario.nBins)
out.dir <- "Test_LR_distribution/scenario.nBins/"
dir.create(out.dir, showWarnings=F, recursive=T)

R <- 1
nr.genes <- 1e+5
mc.cores <- 10
LRlist <- list()
tableList <- list()

# Simulate data + calculate LR statistics

for(sc in 1:length(s)){
  # sc = 1
  cat("*****",s[sc], fill = TRUE)
  
  #### Parameters
  nBins <- scenario.nBins[sc]
  pi.org <- rep(1, nBins)/nBins
  g0.org <- 100
  nM <- nBins*50
  nD <- nM
  sample.size <- 5
  
  simPar <- list(sample.size = sample.size*2, s = s[sc] , pi.org = pi.org , g0.org = g0.org, nr.genes = nr.genes, nM = nM, tot = "uni", nD = nD, out.dir = out.dir)
  
  
  #### Simulate data + run DM
  
  LR <- matrix(0, nrow = simPar$nr.genes, ncol = R)
  
  for(r in 1:R){
    # r=1
    cat("***** R", r, fill=TRUE)
    
    sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, mc.cores=mc.cores, save = FALSE)
    
    sim$dge$samples$group <- as.factor(c(rep("C1", simPar$sample.size/2), rep("C2", simPar$sample.size/2)))
    
    
    name1 <- sim$name1
    g.dir.org <- sim$g.dir.org
    gamma0org <- sum(g.dir.org) 
    
    
    dgeDM <- sim$dge
    dgeDM$commonDispersion <- g0.org
    
    
    dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mc.cores)
    
    dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mc.cores)
    
    tableList[[s[sc]]] <- dgeDM$table
    
    LR[, r] <- dgeDM$table$LR
    
  }
  
  LRlist[[s[sc]]] <- LR
  
  
}



save(tableList, file = paste0(out.dir, "/", "tableList.RData"))


pdf(paste0(out.dir, "/qqplot_LR.pdf"))

for(sc in 1:length(s)){
  
  LR <- LRlist[[s[sc]]]
  
  nBins <- scenario.nBins[sc]
  pi.org <- rep(1, nBins)/nBins
  
  qqplot(rchisq(n = 1e+5, df = length(pi.org)-1), LR, pch=20, cex=3, main= paste("QQ plot \n ", s[sc]), cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Chi-square")
  abline(a = 0, b = 1, col="red", lwd=4)
  
}

dev.off()



for(sc in 1:length(s)){
  
  nBins <- scenario.nBins[sc]
  pi.org <- rep(1, nBins)/nBins
  
  LR <- LRlist[[s[sc]]]
  tb <- tableList[[s[sc]]]
  cat(sum(tb$FDR < 0.05), fill = TRUE)
  
  png(paste0(out.dir, "/qqplot_LR", s[sc],".png"))
  
  qqplot(rchisq(n = 1e+5, df = length(pi.org)-1), LR, pch=20, cex=3, main= paste("QQ plot"), cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Chi-square")
  abline(a = 0, b = 1, col="red", lwd=4)
  
  dev.off()
  
}





#########################################
# Scenario: different estimated dispersion : g0.est
#########################################


### Scenario parameters

scenario.g0.est <- c(0.8, 0.9, 1, 1.1, 1.2)
sname <- "scenario.g0.est"
s <- paste0(sname, scenario.g0.est)
out.dir <- "Test_LR_distribution/scenario.g0.est/"
dir.create(out.dir, showWarnings=F, recursive=T)

R <- 1
nr.genes <- 1e+5
mcCores <- 10
LRlist <- list()
tableList <- list()

# Simulate data + calculate LR statistics

for(sc in 1:length(s)){
  # sc = 1
  cat("*****",s[sc], fill = TRUE)
  
  #### Parameters
  nBins = 3
  pi.org <- rep(1, nBins)/nBins
  g0.org <- 100
  nM <- 150
  nD <- nM
  sample.size <- 5
  
  simPar <- list(sample.size = sample.size*2, s = s[sc] , pi.org = pi.org , g0.org = g0.org, nr.genes = nr.genes, nM = nM, tot = "uni", nD = nD, out.dir = out.dir)
  
  
  #### Simulate data + run DM
  
  LR <- matrix(0, nrow = simPar$nr.genes, ncol = R)
  
  for(r in 1:R){
    # r=1
    cat("***** R", r, fill=TRUE)
    
    sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, mc.cores=mcCores, save = FALSE)
    
    sim$dge$samples$group <- as.factor(c(rep("C1", simPar$sample.size/2), rep("C2", simPar$sample.size/2)))
    
    
    name1 <- sim$name1
    g.dir.org <- sim$g.dir.org
    gamma0org <- sum(g.dir.org) 
    
    
    dgeDM <- sim$dge
    dgeDM$commonDispersion <- g0.org * scenario.g0.est[sc]
    
    
    dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
    
    dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
        
    tableList[[s[sc]]] <- dgeDM$table
    
    LR[, r] <- dgeDM$table$LR
    
  }
  
  LRlist[[s[sc]]] <- LR
  
  
}



save(tableList, file = paste0(out.dir, "/", "tableList.RData"))


pdf(paste0(out.dir, "/qqplot_LR.pdf"))

for(sc in 1:length(s)){
  
  LR <- LRlist[[s[sc]]]
  
  qqplot(rchisq(n = 1e+5, df = length(pi.org)-1), LR, pch=20, cex=3, main= paste("QQ plot \n ", s[sc]), cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Chi-square")
  abline(a = 0, b = 1, col="red", lwd=4)
  
}

dev.off()





for(sc in 1:length(s)){
  
  LR <- LRlist[[s[sc]]]
  tb <- tableList[[s[sc]]]
  cat(sum(tb$FDR < 0.05), fill = TRUE)
  
  png(paste0(out.dir, "/qqplot_LR", s[sc],".png"))
  
  qqplot(rchisq(n = 1e+5, df = length(pi.org)-1), LR, pch=20, cex=3, main= paste("QQ plot"), cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Chi-square")
  abline(a = 0, b = 1, col="red", lwd=4)
  
  dev.off()
  
}



#########################################
# Scenario: different total counts : nM // did not run
#########################################


### Scenario parameters

scenario.nM <- c(25, 75, 150, 300, 1000)
sname <- "scenario.nM"
s <- paste0(sname, scenario.nM)
out.dir <- "Test_LR_distribution/scenario.nM/"
dir.create(out.dir, showWarnings=F, recursive=T)

R <- 1
nr.genes <- 1e+5
mcCores <- 10
LRlist <- list()
tableList <- list()

# Simulate data + calculate LR statistics

for(sc in 1:length(s)){
  # sc = 1
  cat("*****",s[sc], fill = TRUE)
  
  #### Parameters
  nBins <- 3
  pi.org <- rep(1, nBins)/nBins
  g0.org <- 100
  nM <- scenario.nM[sc]
  nD <- nM
  sample.size <- 5
  
  simPar <- list(sample.size = sample.size*2, s = s[sc] , pi.org = pi.org , g0.org = g0.org, nr.genes = nr.genes, nM = nM, tot = "uni", nD = nD, out.dir = out.dir)
  
  
  #### Simulate data + run DM
  
  LR <- matrix(0, nrow = simPar$nr.genes, ncol = R)
  
  for(r in 1:R){
    # r=1
    cat("***** R", r, fill=TRUE)
    
    sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, mc.cores = mcCores, save = FALSE)
    
    sim$dge$samples$group <- as.factor(c(rep("C1", simPar$sample.size/2), rep("C2", simPar$sample.size/2)))
    
    
    name1 <- sim$name1
    g.dir.org <- sim$g.dir.org
    gamma0org <- sum(g.dir.org) 
    
    
    dgeDM <- sim$dge
    dgeDM$commonDispersion <- g0.org
    
    
    dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
    
    dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
    
    tableList[[s[sc]]] <- dgeDM$table
    
    LR[, r] <- dgeDM$table$LR
    
  }
  
  LRlist[[s[sc]]] <- LR
  
  
}



save(tableList, file = paste0(out.dir, "/", "tableList.RData"))


pdf(paste0(out.dir, "/qqplot_LR.pdf"))

for(sc in 1:length(s)){
  
  LR <- LRlist[[s[sc]]]
  
  qqplot(rchisq(n = 1e+5, df = length(pi.org)-1), LR, pch=20, cex=3, main= paste("QQ plot \n ", s[sc]), cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Chi-square")
  abline(a = 0, b = 1, col="red", lwd=4)
  
}

dev.off()


#########################################
# Scenario: dispersion ~ N(100, sigma) ; in DM using common dispersion // did not run
#########################################


### Scenario parameters

scenario.sigma <- c(1, 3, 6, 10, 20)
sname <- "scenario.sigma"
s <- paste0(sname, scenario.sigma)
out.dir <- "Test_LR_distribution/scenario.sigma/"
dir.create(out.dir, showWarnings=F, recursive=T)

R <- 1
nr.genes <- 1e+5
mcCores <- 10
LRlist <- list()
tableList <- list()

# Simulate data + calculate LR statistics

for(sc in 1:length(s)){
  # sc = 1
  cat("*****",s[sc], fill = TRUE)
  
  #### Parameters
  nBins = 3
  pi.org <- rep(1, nBins)/nBins
  g0.mean <- 100
  g0.org <- abs(rnorm(nr.genes,  g0.mean, scenario.sigma[sc]))
  nM <- 150
  nD <- nM
  sample.size <- 5
  
  simPar <- list(sample.size = sample.size*2, s = s[sc] , pi.org = pi.org , g0.org = g0.org, nr.genes = nr.genes, nM = nM, tot = "uni", nD = nD, out.dir = out.dir)
  
  
  #### Simulate data + run DM
  
  LR <- matrix(0, nrow = simPar$nr.genes, ncol = R)
  
  for(r in 1:R){
    # r=1
    cat("***** R", r, fill=TRUE)
    
    sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, mc.cores=mcCores, save = FALSE)
    
    sim$dge$samples$group <- as.factor(c(rep("C1", simPar$sample.size/2), rep("C2", simPar$sample.size/2)))
    
    
    name1 <- sim$name1
    g.dir.org <- sim$g.dir.org
    gamma0org <- sum(g.dir.org) 
    
    
    dgeDM <- sim$dge
    dgeDM$commonDispersion <-  g0.mean
    
    
    dgeDM <- dmFit(dgeDM, group=NULL, dispersion="commonDispersion", mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
    
    dgeDM <- dmTest(dgeDM, mode="constrOptim2G", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)
    
    tableList[[s[sc]]] <- dgeDM$table
    
    LR[, r] <- dgeDM$table$LR
    
  }
  
  LRlist[[s[sc]]] <- LR
  
  
}


save(tableList, file = paste0(out.dir, "/", "tableList.RData"))



pdf(paste0(out.dir, "/qqplot_LR.pdf"))

for(sc in 1:length(s)){
  
  LR <- LRlist[[s[sc]]]
  
  qqplot(rchisq(n = 1e+5, df = length(pi.org)-1), LR, pch=20, cex=3, main= paste("QQ plot \n ", s[sc]), cex.main=1.5, cex.lab=1.45, cex.axis = 1.5, xlab="Chi-square")
  abline(a = 0, b = 1, col="red", lwd=4)
  
}

dev.off()




















