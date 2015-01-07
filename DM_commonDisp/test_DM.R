
##############################################################################

# BioC 14
# Created 12 Aug 2014
# Updated 12 Aug 2014

# smulate data from known DIR-MULTI 
#   - estimation with dirmult()
#   - estimation with equalTheta() - common dispersion
# MSE

# update 12 Aug: 

##############################################################################

setwd("/home/gosia/Multinomial_project/DM_adjusted")
Rdir <- "/home/gosia/R/R_Multinomial_project/DM_adjusted/"


library(dirmult)
library(edgeR)
library(parallel)

out.dir <- "PLOTS/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)

######################################
# simulate data
######################################

##### function
sample.size = 3; s = "1"; pi.org = c(0.3, 0.3, 0.3, 0.1); g0.org = 30; nr.genes = 100; nM = 200; tot="uni"; nD=3; mc.cores=30

simulate_from_DM <- function(sample.size = 3, s = "1", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 30, nr.genes = 100, nM = 200, tot="uni", nD=3, out.dir, mc.cores=30){
  
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
  genes <- paste0("g", rep(i, length(g.dir.org)))   
  rownames(d) <- genes
  
  return(d)
  
}, mc.cores=mc.cores)
  
sim <- do.call(rbind, sim)
genes <- rownames(sim)
  dge <- DGEList( counts=sim, group = rep(1, sample.size), genes=data.frame(gene.id=genes, ete.id = paste0(genes, ":e", rep(1:length(g.dir.org), nr.genes)) ) )
  
  dge
  
  ## add 1
  # dge$counts <- dge$counts + 1 
  # dge$counts[ dge$counts == 0 ] <- 1
  
  name1 <- paste0("SIM",s,"_ss",sample.size,"_g0", g0.org,"_T", tot, "_nM",nM ,"_")
  
  save(dge, name1, g.dir.org, file=paste0(out.dir, "/", name1, "dge",".RData"))
  
  return(list(dge=dge, name1=name1, g.dir.org=g.dir.org))
  
  
}


##### run

sim <- simulate_from_DM(sample.size = 3, s = "1", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 30, nr.genes = 10000, nM = 200, tot="uni", nD=3, out.dir=out.dir)


# sim <- simulate_from_DM(sample.size = 3, s = "1", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 30, nr.genes = 100000, nM = 200, tot="nbinom", nD=3, out.dir=out.dir)


dge <- sim$dge
name1 <- sim$name1
g.dir.org <- sim$g.dir.org

######################################
# run dirmult
######################################
## dirmult estimates gammas!!!

##### function

source(paste0(Rdir, "dirmult_code.R"))


fit_DM <- function(dge, name1, out.dir, mc.cores=20){
  
  split_counts <- split(seq_len(nrow(dge$counts)), dge$genes$gene.id, drop=TRUE)
  gene.id.unique <- names(split_counts)
  group <- dge$samples$group
  group <- as.factor(group)
  
  
  FitDM <- mclapply(seq_len(length(gene.id.unique)), function(g){
    # g=9992
    print(g)
    counts.tmp <- dge$counts[split_counts[[g]], , drop=FALSE]
    
    null.list <- list(loglikh=NULL, Y=NULL, gh=NULL, g0h = NULL ,pih=NULL, th=NULL, df=NULL , gMom=NULL , g0Mom=NULL, piMom=NULL, tMom=NULL, meanY=NULL)
    
    if(nrow(counts.tmp)==1)
      return(null.list)
    
    Y <- t(counts.tmp)
    colnames(Y) <- dge$genes$ete.id[split_counts[[g]]]
    q <- ncol(Y)
    
    if(any(rowSums(Y) == 0))
      return(null.list)
    
    null <- NULL
    try (null <- dirmult( Y , trace=FALSE), silent=TRUE)
    
    if(is.null(null))
      return(null.list)
    
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
  
  
  pih <- grep("pih", params.n) 
  pih <- do.call(rbind, params[pih])
  colnames(pih) <- paste0("e", 1:ncol(pih))
  
  table <- data.frame(gene.id=gene.id.unique, loglikh=loglikh, g0h = g0h, th=th, meanY=meanY,pih, row.names = NULL)
  
  save(FitDM, table, file=paste0(out.dir, "/", name1, "FitDM",".RData"))
  
  return(list(FitDM=FitDM, table=table))
  
}


##### run


fit <- fit_DM(dge, name1=name1, out.dir=out.dir, mc.cores=30)
  




FitDM <- fit$FitDM
tableDM <- table <- fit$table
  


######################################
# plot g0 and pi estimates
######################################
library(beanplot)

name <- "dirmult"


pdf(paste0(out.dir,"/", name1, "ESTIM_g0_",name,".pdf"))

beanplot(table$g0h, col="orangered", log="y", what=c(0,1,1,0), main=paste0("gamma: ", paste(g.dir.org, collapse=", ")), ylim=c(10, 1e+7))
abline(h=sum(g.dir.org), lty=2)

smoothScatter(table$meanY, log(table$g0h), col="orangered", main=paste0("gamma: ", paste(g.dir.org, collapse=", ")), xlab="mean gene expression", ylab="log gamma +", nrpoints = Inf, colramp=colorRampPalette(c("white", blues9)), ylim=c(log(10),log(1e+7)))
abline(h=log(sum(g.dir.org)), lty=2)

beanplot(table[, paste0("e", 1:ncol(pih))], col="cadetblue", what=c(0,1,1,0), main=paste0("pi: ", paste(pi.org, collapse=", ")), ylim=c(0,1), log = "")
abline(h=pi.org, lty=2)


dev.off()
  
  
######################################
# estimate common dispersion with equalTheta
######################################

t.org <- 1/(1+g0.org)
t.st <- mean(tableDM[, "th"])


fit_equalTheta <- function(dge, t, name1, out.dir){
  
  rownames(dge$counts) <- dge$genes$ete.id
  
  split_counts <- lapply(split(as.data.frame(dge$counts), dge$genes$gene.id, drop=TRUE), t) # t - transpose
  
  gene.id.unique <- names(split_counts)
  group <- dge$samples$group
  group <- as.factor(group)
  
  eT <- equalTheta(split_counts, theta=t)
  
  th <- eT$theta[[1]]
  g0h <- (1-th)/th
  
  pih <- do.call(rbind, eT$pi)
  colnames(pih) <- paste0("e", 1:ncol(pih))
  rownames(pih) <- gene.id.unique
  
  loglikh <- eT$loglik
  
  table <- data.frame(gene.id=gene.id.unique, loglikh=loglikh, g0h = g0h, th=th, pih, row.names = NULL)
  
  save(eT, table, file=paste0(out.dir, "/", name1, "eT",".RData"))
  
  return(list(eT = eT, table = table))
  
}


# fit_et <- fit_equalTheta(dge, t=t.st, name1, out.dir)
fit_et <- fit_equalTheta(dge, t=0.03, name1, out.dir)


tableET <- table <- fit_et$table

name <- "equalTheta"


pdf(paste0(out.dir,"/", name1, "ESTIM_g0_",name,".pdf"))

beanplot(table[, paste0("e", 1:ncol(pih))], col="cadetblue", what=c(0,1,1,0), main=paste0("pi: ", paste(pi.org, collapse=", ")), ylim=c(0,1), log = "")
abline(h=pi.org, lty=2)

dev.off()



######################################
# MSE
######################################


mse <- function(est, org){
  colSums((est-org)^2)/nrow(est) 
}




mse(tableDM[, paste0("e", 1:ncol(pih))], pi.org)

mse(tableET[, paste0("e", 1:ncol(pih))], pi.org)



######################################
# MSE(pi) vs. theta
######################################

name <- "dirmult"



pdf(paste0(out.dir,"/", name1, "MSEvsGamma0_",name,".pdf"))



plot(1,1, type="n", xlim=c(log(10),log(1e+7)), ylim=c(min(tableDM[, paste0("e", 1:ncol(pih))], tableET[, paste0("e", 1:ncol(pih))]), max(tableDM[, paste0("e", 1:ncol(pih))], tableET[, paste0("e", 1:ncol(pih))])), xlab="Gamma0", ylab="Pi", main=paste0("pi: ", paste(pi.org, collapse=", "))) 

for(i in 1:ncol(tableDM[, paste0("e", 1:ncol(pih))]))
  points(log(tableDM[, "g0h"]), tableDM[, paste0("e", i)], col=i+1, pch=".:", cex=2)

legend("topright",paste0("e", 1:ncol(pih)), col=c(1:ncol(pih)+1), pch=19)




plot(1,1, type="n", xlim=c(log(10),log(1e+7)), ylim=c(-0.2, 0.2), xlab="Gamma0", ylab="(er_DM = Pi.h - Pi)") 

for(i in 1:ncol(tableDM[, paste0("e", 1:ncol(pih))]))
  points(log(tableDM[, "g0h"]), tableDM[, paste0("e", i)] - pi.org[i], col=i+1, pch=".:", cex=2)

legend("topright",paste0("e", 1:ncol(pih)), col=c(1:ncol(pih)+1), pch=19)




plot(1,1, type="n", xlim=c(log(10),log(1e+7)), ylim=c(-0.02, 0.02), xlab="Gamma0", ylab="er_DM - er_ET") 

for(i in 1:ncol(tableDM[, paste0("e", 1:ncol(pih))]))
  points(log(tableDM[, "g0h"]), tableDM[, paste0("e", i)] - pi.org[i] - (tableET[, paste0("e", i)] - pi.org[i]), col=i+1, pch=".:", cex=2)
abline(v=log(tableET[1, "g0h"]), lty=2)
legend("topright",paste0("e", 1:ncol(pih)), col=c(1:ncol(pih)+1), pch=19)




dev.off()



















































