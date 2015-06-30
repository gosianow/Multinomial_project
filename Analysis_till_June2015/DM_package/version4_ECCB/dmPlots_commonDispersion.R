
##############################################################################

# BioC 14
# Created 06 Sep 2014:

# simulation form DM

# common dispersion plots for ECCB W1

##############################################################################


setwd("/home/gosia/Multinomial_project/DM_package")

library(edgeR)
library(parallel)
library(dirmult)

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version4/"
source(paste0(Rdir, "dmFunctions_v4.R"))

out.dir <- "PLOTS4/commonDispersion_ECCB/"
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
# simPar <- list(sample.size = 3, s = paste0(1) , pi.org = pi.org , g0.org = g0.org, nr.genes = 1000, nM = nM, tot = "nbinom", nD = nD, out.dir = out.dir)


simPar <- list(sample.size = 15, s = paste0(2) , pi.org = pi.org , g0.org = g0.org, nr.genes = 1000, nM = nM, tot = "nbinom", nD = nD, out.dir = out.dir)


###################### Simulate data

# R <- 50
R <- 1 

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
  # r=25
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



pdf(paste0(out.dir, "/", name1, "MSEcommonDisp_PRES",".pdf"))

beanplot(SE, what=c(0,1,1,0), main=paste0("Standard Error","\n gamma: ", paste(g.dir.org, collapse=", ")), ylim=c(min(SE)-gamma0org^0.5, max(SE)+gamma0org^0.5), cex=3, cex.lab=2, cex.axis=1.5, cex.main=1.5, col = list(c("chartreuse3", "chartreuse3", "chartreuse3", "black"), c("orange", "orange", "orange", "black")) )
abline(h=0, lty=2, lwd=4)

dev.off()




############################################################################
# tagwise dispersion from dirmult 
############################################################################


##### function

source("/home/gosia/R/R_Multinomial_project/DM_weighted_likelihood/dirmult_code.R")

fit_DM <- function(dge, name1, out.dir, mc.cores=20){
  
  split_counts <- split(seq_len(nrow(dge$counts)), dge$genes$gene.id, drop=TRUE)
  gene.id.unique <- names(split_counts)
  group <- dge$samples$group
  group <- as.factor(group)
  
  
  FitDM <- mclapply(seq_len(length(gene.id.unique)), function(g){
    # g=9992
    print(g)
    counts.tmp <- dge$counts[split_counts[[g]], , drop=FALSE]
    
    if(nrow(counts.tmp)==1)
      return(list(loglikh=NULL, Y=NULL, gh=NULL, g0h = NULL ,pih=NULL, th=NULL, df=NULL , gMom=NULL , g0Mom=NULL, piMom=NULL, tMom=NULL, meanY=NULL))
    
    Y <- t(counts.tmp)
    colnames(Y) <- dge$genes$ete.id[split_counts[[g]]]
    q <- ncol(Y)
    
    if(any(rowSums(Y) == 0))
      return(list(loglikh=NULL, Y=NULL, gh=NULL, g0h = NULL ,pih=NULL, th=NULL, df=NULL , gMom=NULL , g0Mom=NULL, piMom=NULL, tMom=NULL, meanY=NULL))
    
    null <- NULL
    try (null <- dirmult( Y , trace=FALSE), silent=TRUE)
    
    if(is.null(null))
      return(list(loglikh=NULL, Y=NULL, gh=NULL, g0h = NULL ,pih=NULL, th=NULL, df=NULL , gMom=NULL , g0Mom=NULL, piMom=NULL, tMom=NULL, meanY=NULL))
    
    meanY <- mean(rowSums(Y))
    
    loglikh <- null$loglik
    gh <- null$gamma
    names(gh) <- dge$genes$ete.id[split_counts[[g]]]
    g0h <- sum(gh)
    pih <- null$pi
    th <- null$theta
    names(pih) <- dge$genes$ete.id[split_counts[[g]]]
    
    gMom <- null$g.mom
    g0Mom <- null$g0.mom
    tMom <- null$t.mom
    piMom <- null$pi.mom
    
    gc()
    
    return(list(loglikh=loglikh, Y=Y, gh=gh, g0h = g0h ,pih=pih, th=th, df=(q-1) , gMom=gMom , g0Mom=g0Mom, piMom=piMom, tMom=tMom, meanY=meanY))
    
  }, mc.cores=mc.cores)
  
  names(FitDM) <- gene.id.unique
  
  save(FitDM, file=paste0(out.dir, "/", name1, "FitDM",".RData"))
  
  
  params <- unlist(FitDM, recursive=FALSE)
  params.n <- names(params)
  
  g0h <- grep("g0h", params.n) 
  g0h <- unlist(params[g0h])
  
  th <- grep("th", params.n) 
  th <- unlist(params[th])
  
  loglikh <- grep("loglikh", params.n) 
  loglikh <- unlist(params[loglikh])
  
  meanY <- grep("meanY", params.n) 
  meanY <- unlist(params[meanY])
  
  table <- data.frame(gene.id=gene.id.unique, loglikh=loglikh, g0h = g0h, th=th, meanY=meanY, row.names = NULL)
  
  save(FitDM, table, file=paste0(out.dir, "/", name1, "FitDM",".RData"))
  
  return(list(FitDM=FitDM, table=table))
  
}


##### run

r <- 1

dge <- simList[[r]]$dge 

fit <- fit_DM(dge, name1=name1, out.dir=out.dir, mc.cores=40)

FitDM <- fit$FitDM
table <- fit$table



##### plot g0 estimates

library(beanplot)

name <- "dirmult"



pdf(paste0(out.dir,"/", name1, "tagwiseDispersion",name,".pdf"))

plot(table$meanY, log(table$g0h), col="slateblue2", lwd=6, cex=1.5, cex.lab=1.45, cex.axis=1.5, xlab="mean gene expression", ylab="log gamma +", pch=19)
abline(h = log(gamma0org), col=1, lwd=10, lty=2)

# plot(table$meanY, log(table$g0h), col="slateblue2", lwd=6, cex=1.5, cex.lab=1.45, cex.axis=1.5, xlab="mean gene expression", ylab="log gamma +", pch=19)
# abline(h = log(gamma0org), col=1, lwd=10, lty=2)
# abline(h = log(commonDipEsts[r, "constrOptim2"]), col="chartreuse3", lwd=8)
# abline(h = log(commonDipEsts[r, "constrOptim2ADJ"]), col="orange", lwd=8)

# beanplot(table$g0h, col="orangered", log="y", what=c(1,1,1,0), main=paste0("gamma: ", paste(g.dir.org, collapse=",")), ylim=c(10, 1e+7))
# abline(h=sum(g.dir.org), lty=2)
# 
# smoothScatter(table$meanY, log(table$g0h), col="orangered", main=paste0("gamma: ", paste(g.dir.org, collapse=",")), xlab="mean gene expression", ylab="log gamma +", nrpoints = Inf, colramp=colorRampPalette(c("white", blues9)), ylim=c(log(10),log(1e+7)))
# abline(h=log(sum(g.dir.org)), lty=2)

dev.off()









