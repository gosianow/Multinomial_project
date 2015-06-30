
##############################################################################

# BioC 14
# Created 6 Aug 2014
# Updated 12 Aug 2014

# smulate data from known DIR-MULTI 

# test dirmult performance for different n

##############################################################################

setwd("/home/gosia/Multinomial_project/DM_weighted_likelihood")


library(dirmult)
library(edgeR)
library(parallel)

out.dir <- "PLOTS_dirmult_performance/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)

######################################
# simulate data
######################################

##### function
# sample.size = 3; s = "1"; pi.org = c(0.3, 0.3, 0.3, 0.1); g0.org = 30; nr.genes = 100; nM = 200; tot="uni"; nD=3; mc.cores=30

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

sim <- simulate_from_DM(sample.size = 3, s = "1", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 30, nr.genes = 100000, nM = 200, tot="uni", nD=3, out.dir=out.dir)


sim <- simulate_from_DM(sample.size = 3, s = "1", pi.org = c(0.3, 0.3, 0.3, 0.1), g0.org = 30, nr.genes = 100000, nM = 200, tot="nbinom", nD=3, out.dir=out.dir)


dge <- sim$dge
name1 <- sim$name1
g.dir.org <- sim$g.dir.org

######################################
# run dirmult
######################################
## dirmult estimates gammas!!!

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


fit <- fit_DM(dge, name1=name1, out.dir=out.dir, mc.cores=30)
  
FitDM <- fit$FitDM
table <- fit$table
  
######################################
# plot g0 estimates
######################################
library(beanplot)


name <- "dirmult"

pdf(paste0(out.dir,"/", name1, "ESTIM_g0_",name,".pdf"))

beanplot(table$g0h, col="orangered", log="y", what=c(1,1,1,0), main=paste0("gamma: ", paste(g.dir.org, collapse=",")), ylim=c(10, 1e+7))
abline(h=sum(g.dir.org), lty=2)

smoothScatter(table$meanY, log(table$g0h), col="orangered", main=paste0("gamma: ", paste(g.dir.org, collapse=",")), xlab="mean gene expression", ylab="log gamma +", nrpoints = Inf, colramp=colorRampPalette(c("white", blues9)), ylim=c(log(10),log(1e+7)))
abline(h=log(sum(g.dir.org)), lty=2)


dev.off()
  
  















