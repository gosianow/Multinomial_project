##############################################################################

# BioC 14
# Created 19 Aug 2014



# Updated 22 Aug 2014:


##############################################################################

setwd("/home/gosia/Multinomial_project/DM_package")
Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/"

library(edgeR)
library(parallel)
library(dirmult)

out.dir <- "PLOTS/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)

######################################
# simulate data from DM
######################################

##### function

simulate_from_DM <- function(sample.size = 3, s = "1", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 30, nr.genes = 100, nM = 200, tot="uni", nD=3, out.dir, mc.cores=20, save=TRUE){
  
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
  
  dge
  
  ## add 1
  # dge$counts <- dge$counts + 1 
  # dge$counts[ dge$counts == 0 ] <- 1
  
  name1 <- paste0("SIM",s,"_ss",sample.size,"_g0", g0.org,"_T", tot, "_nM",nM ,"_")
  
  if(save)
  save(dge, name1, g.dir.org, file=paste0(out.dir, "/", name1, "dge",".RData"))
  
  return(list(dge=dge, name1=name1, g.dir.org=g.dir.org))
  
  
}


##### run

# sim <- simulate_from_DM(sample.size = 3, s = "1", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 30, nr.genes = 100, nM = 200, tot="uni", nD=3, out.dir=out.dir)


#sim <- simulate_from_DM(sample.size = 100, s = "1", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 30, nr.genes = 100, nM = 300, tot="nbinom", nD=3, out.dir=out.dir)



sim <- simulate_from_DM(sample.size = 10, s = "1", pi.org = c(0.7, 0.3), g0.org = 30, nr.genes = 10, nM = 200, tot="nbinom", nD=3, out.dir=out.dir)



dge <- sim$dge
name1 <- sim$name1
g.dir.org <- sim$g.dir.org


######################################
# 
######################################


y <- split(as.data.frame(dge$counts), dge$genes$gene.id)
gamma0 <- rep(sum(g.dir.org), length(y))

y <- y[[2]]
gamma0 <- gamma0[1]
pi <- rowSums(y)/sum(y)

y <- t(y)

source("/home/gosia/R/R_Multinomial_project/DM_adjusted/dirmult_code.R")
source(paste0(Rdir, "dmFunctions.R"))



dmLogLik(pi, gamma0, y)
dmLogLikG(pi, gamma0, y)


loglik(x = y, t = pi*gamma0)


dmFitGene(t(y), gamma0, mode = "constrOptim")

dmFitGene(t(y), gamma0, mode = "obs", epsilon=1e-05)


######################################
# plots of likelihood functions
######################################

### plot of likelihood function

piX <- seq(0, 1, by = 0.01)
piX <- piX[c(-1, -length(piX))]
piX <- cbind(piX, 1-piX)


f <- rep(0, nrow(piX))
for(i in 1:length(f))
  f[i] <- dmLogLik(piX[i,], gamma0, y)


fG <- rep(0, nrow(piX))
for(i in 1:length(fG))
  fG[i] <- dmLogLikG(piX[i,], gamma0, y)


fGkm1 <- rep(0, nrow(piX))
for(i in 1:length(fGkm1))
  fGkm1[i] <- dmLogLikGkm1(piX[i,1], gamma0, y)



pdf(paste0(out.dir,"logLik.pdf"))

plot(piX[,1], f, type="l", col="green", lwd=4)
lines(piX[,2], f, type="l", col="red", lwd=4)

lines(piX[,1], fG, type="l", col="darkgreen", lty=2, lwd=4)
lines(piX[,2], fG, type="l", col="darkred", lty=2, lwd=4)

lines(piX[,1], fGkm1, type="l", col="pink", lty=3, lwd=4)

dev.off()


######################################
# plots of score functions
######################################

piX <- seq(0, 1, by = 0.01)
piX <- piX[c(-1, -length(piX))]
piX <- cbind(piX, 1-piX)


### plot of score function

ff <- matrix(0, nrow(piX), 2)
for(i in 1:nrow(ff))
  ff[i,] <- dmScoreFun(piX[i,], gamma0, y)


ffG <- matrix(0, nrow(piX), 2)
for(i in 1:nrow(ffG))
  ffG[i,] <- dmScoreFunG(piX[i,], gamma0, y)


ffGkm1 <- matrix(0, nrow(piX), 1)
for(i in 1:nrow(ffGkm1))
  ffGkm1[i,] <- dmScoreFunGkm1(piX[i,1], gamma0, y)


ffkm1 <- matrix(0, nrow(piX), 1)
for(i in 1:nrow(ffkm1))
  ffkm1[i,] <- dmScoreFunkm1(piX[i,1], gamma0, y)


pdf(paste0(out.dir,"score.pdf"))

plot(piX[,1], ff[,1], type="l", col="green", ylim=c(-2000, 2000), lwd=4)
lines(piX[,2], ff[,2], type="l", col="red", lwd=4)

lines(piX[,1], ffG[,1], type="l", col="darkgreen", lty=2, lwd=4)
lines(piX[,2], ffG[,2], type="l", col="darkred", lty=2, lwd=4)

lines(piX[,1], ffGkm1, type="l", col="pink", lty=3, lwd=4)
lines(piX[,1], ffkm1, type="l", col="magenta", lty=3, lwd=5)

abline(h = 0, lty=2)

dev.off()



### plot of score function from dirmult


ff <- matrix(0, nrow(piX), 2)

for(i in 1:nrow(ff))
  ff[i,] <- u(y, piX[i,]*gamma0)


pdf(paste0(out.dir,"score_dirmult.pdf"))

plot(piX[,1], ff[,1], type="l", col="green")
lines(piX[,2], ff[,2], type="l", col="red")
abline(h = 0, lty=2)

dev.off()




######################################
# common dispersion estimates
######################################

source("/home/gosia/R/R_Multinomial_project/DM_adjusted/dirmult_code.R")
source(paste0(Rdir, "dmFunctions.R"))

y <- split(as.data.frame(dge$counts), dge$genes$gene.id)

dmEstimateCommonDisp(y, adjust = FALSE, mode = "constrOptim", mcCores=20, interval = c(0, 1e+5), tol = 1e-05)


### plot profile log lik
gamma0Int <- seq(0, 200, by=1)[-1]
profLik <- sapply(gamma0Int, dmAdjustedProfileLik, y = y, adjust = FALSE, mode = "constrOptim", mcCores = 20, common = TRUE)
 

pdf(paste0(out.dir,"profLogLik.pdf"))
plot(gamma0Int, profLik, type="l")
dev.off()




yt <- lapply(split(as.data.frame(dge$counts), dge$genes$gene.id, drop=TRUE), t) # t - transpose

et <- equalTheta(yt, theta=round(1/(1+g0.org), 2))

sum(et$gamma[[1]])



######################################
# MSE of common dispersion estimates
######################################
library(beanplot)

source("/home/gosia/R/R_Multinomial_project/DM_adjusted/dirmult_code.R")
source(paste0(Rdir, "dmFunctions.R"))



simPar <- list(sample.size = 10, s = "1", pi.org = c(0.7, 0.3), g0.org = 30, nr.genes = 100, nM = 200, tot = "nbinom", nD = 3, out.dir = out.dir)


simPar <- list(sample.size = 3, s = "2", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 30, nr.genes = 100, nM = 200, tot = "nbinom", nD = 3, out.dir = out.dir)


simPar <- list(sample.size = 3, s = "2", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 3, nr.genes = 100, nM = 200, tot = "nbinom", nD = 3, out.dir = out.dir)


simPar <- list(sample.size = 3, s = "3", pi.org = c(0.3, 0.2, 0.5), g0.org = 10, nr.genes = 100, nM = 200, tot = "nbinom", nD = 3, out.dir = out.dir)


R <- 50

commonDipEsts <- matrix(0, nrow = R, ncol = 3)
colnames(commonDipEsts) <- c("equalTheta", "constrOptim", "FisherObs")

for(r in 1:R){
  # r=25
  print(r)
  sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nD, out.dir = simPar$out.dir, save = FALSE)
  dge <- sim$dge

  theta <- round(1 / (1 + simPar$g0.org), 2)
  yt <- lapply(split(as.data.frame(dge$counts), dge$genes$gene.id, drop=TRUE), t) # t - transpose
  et <- equalTheta(yt, theta=theta, trace = FALSE)
  commonDipEsts[r, 1] <- sum(et$gamma[[1]])
  
  
  y <- split(as.data.frame(dge$counts), dge$genes$gene.id) 
  commonDipEsts[r, 2] <- dmEstimateCommonDisp(y, adjust = FALSE, mode = "constrOptim", mcCores=20, interval = c(0, 1e+5), tol = 1e-05)
   
  y <- split(as.data.frame(dge$counts), dge$genes$gene.id) 
  commonDipEsts[r, 3] <- dmEstimateCommonDisp(y, adjust = FALSE, mode = "obs", mcCores=20, interval = c(0, 1e+5), tol = 1e-05)
  
  
}

commonDipEsts <- commonDipEsts[1:24,]

name1 <- sim$name1
g.dir.org <- sim$g.dir.org

commonDipEsts <- as.data.frame(commonDipEsts)
SE <- commonDipEsts - sum(g.dir.org)
  
  



pdf(paste0(out.dir,"/commonDisp_", name1, ".pdf"))

beanplot(commonDipEsts, col="orangered", what=c(0,1,1,0), main=paste0("gamma0H","\n gamma: ", paste(g.dir.org, collapse=", ")), ylim=c(min(commonDipEsts)-1, max(commonDipEsts)+1))
abline(h=sum(g.dir.org), lty=2)

beanplot(SE, col="deepskyblue2", what=c(0,1,1,0), main=paste0("SE: gamma0H - gamma0","\n gamma: ", paste(g.dir.org, collapse=", ")), ylim=c(min(SE)-1, max(SE)+1))
abline(h=0, lty=2)

dev.off()




















