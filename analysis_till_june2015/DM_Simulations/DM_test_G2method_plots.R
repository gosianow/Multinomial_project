##############################################################################
# smulate data from known DIR-MULTI 
# profile likelihood plot
##############################################################################

# source("/home/gosia/R/R_Multinomial_project/Simulate_DM_test_G2method.R")

##############################################################################

setwd("/home/gosia/Multinomial_project/DM_Simulations")

# install.packages("/home/gosia/R/packages/dirmult_0.1.3-4.tar.gz", lib="/home/gosia/R/libraries/3.1.0/")


library("dirmult")
library("gtools")
library(edgeR)
library(parallel)


######################################
# simulate data
######################################

sample.size <- 30
s <- "1"

#name <- paste0("SIM",sample.size,"_s",s,"_init_NOTexact_")
#name <- paste0("SIM",sample.size,"_s",s,"_initscalar_NOTexact_")

#name <- paste0("SIM",sample.size,"_s",s,"_init_exact_")
#name <- paste0("SIM",sample.size,"_s",s,"_initscalar_exact_")

name1 <- paste0("SIM",sample.size,"_s",s,"_")


# simulate some dirichlets
g.dir.org <- cbind(c(90, 5, 5), c(30, 30, 30) ) # s <- "1"
# g.dir.org <- cbind(c(9, 0.5, 0.5), c(3, 3, 3) ) # s <- "2"
# g.dir.org <- cbind(c(0.9, 0.05, 0.05), c(0.3, 0.3, 0.3) ) # s <- "3"

g.dir <- rbind( rdirichlet( sample.size, g.dir.org[,1] ), rdirichlet( sample.size, g.dir.org[,2] ) )

sim <- NULL
genes <- NULL
nr.genes <- 1000

for(i in 1:nr.genes){
  # i=1
  tot <- rnbinom(2*sample.size, mu=100, size=3)
  # tot <- round(rnorm(2*sample.size, mean=200, sd=20))  
  # tot <- rep(200, 2*sample.size)
  d <- sapply(1:(2*sample.size), function(u) rmultinom(1, prob=g.dir[u,], size=tot[u]))  
  sim <- rbind(sim, d)  
  genes <- c( genes, paste0("g", rep(i, 3))  )  
}

rownames(sim) <- genes

dge <- DGEList( counts=sim, group = c(rep(1, sample.size), rep(2, sample.size)), genes=data.frame(gene.id=genes, ete.id =genes) )

## or add 1
#dge$counts <- dge$counts + 1 
dge$counts[ dge$counts == 0 ] <- 1

save(dge, file=paste0("PLOTS/", name1, "dge_",gsub(":", "-" ,gsub(" ", "_", as.character(date()))),".RData"))

# load(paste0("PLOTS/","SIM3_s1_dge_Wed_Jun_11_10-20-16_2014.RData"))


######################################
# run fitting
######################################



source("/home/gosia/R/R_Multinomial_project/dirmult.R")
source("/home/gosia/R/R_Multinomial_project/glmFitDM.R")


# define initial values
init <- g.dir.org - 1
initscalar <- colSums(g.dir.org) - 5

# run fitting
g2fit.obs <- g2FitDM(dge)

# g2fit.init <- g2FitDM(dge, init=init)
# g2fit.inits <- g2FitDM(dge, initscalar=initscalar)

#g2fit.exp <- g2FitDM(dge, mode="exp") # does not work in this mode



####

plot.trend <- function(g2fit, name){
  
  gh0.tmp <-  mclapply(names(g2fit$FitDM), function(g){
    
    if(!is.null(g2fit$FitDM[[g]])){           
      
      meanY0 <- mean(rowSums(g2fit$FitDM[[g]]$Y[1:sample.size,]))
      meanY1 <- mean(rowSums(g2fit$FitDM[[g]]$Y[(sample.size+1):(2*sample.size),]))
      gh0.0f <- mean(g2fit$FitDM[[g]]$gh0[1])
      gh0.1f <- mean(g2fit$FitDM[[g]]$gh0[sample.size+1])
      mom1 <- g2fit$FitDM[[g]]$mom[1]
      mom2 <- g2fit$FitDM[[g]]$mom[2]

      return(data.frame(gh0.0f=gh0.0f, gh0.1f=gh0.1f,  meanY0=meanY0, meanY1=meanY1, mom1=mom1, mom2=mom2))
    }    
    return(data.frame(gh0.0f=NA, gh0.1f=NA, meanY0=NA, meanY1=NA, mom1=NA, mom2=NA))
    
  }, mc.cores=10)
  
  
  gh0 <- do.call(rbind, gh0.tmp)
  
  gh0 <- data.frame(GeneID=names(g2fit$FitDM), gh0)
  head(gh0)
  table.gh0 <- gh0
  
  
  library(beanplot)
  
  
  pdf(paste0("PLOTS/", name, "ESTIMATION_gh0",gsub(":", "-" ,gsub(" ", "_", as.character(date()))),".pdf"), h=5, w=10)
  par(mfrow=c(1,2 ), cex.main=1.1, cex.axis = 1.1, cex=1.2)
  
  
  plot(table.gh0$meanY0, table.gh0$gh0.0f, col=2, log="xy", main=paste0("gamma: ", paste(g.dir.org[,1], collapse=",")))
  abline(h=sum(g.dir.org[,1]), lty=2)
  plot(table.gh0$meanY1 , table.gh0$gh0.1f, col=3, log="xy", main=paste0("gamma: ", paste(g.dir.org[,2], collapse=",")))
  abline(h=sum(g.dir.org[,2]), lty=2)
  
  plot((1-table.gh0$mom1)/table.gh0$mom1 , table.gh0$gh0.0f, col=2, log="xy", main="log( (1-mom)/mom ) vs log( gh0 )")
  plot((1-table.gh0$mom2)/table.gh0$mom2 , table.gh0$gh0.1f, col=3, log="xy")
  
  beanplot(table.gh0$gh0.0f, col=2, log="y", what=c(1,1,1,0), main=paste0("gamma: ", paste(g.dir.org[,1], collapse=",")))
  abline(h=sum(g.dir.org[,1]), lty=2)
  beanplot(table.gh0$gh0.1f, col=3, log="y", what=c(1,1,1,0), main=paste0("gamma: ", paste(g.dir.org[,2], collapse=",")))
  abline(h=sum(g.dir.org[,2]), lty=2)
  
  dev.off()
  
  
  pdf(paste0("PLOTS/", name, "ESTIMATION_gh0",".pdf"))
  
  plot(table.gh0$meanY0, table.gh0$gh0.0f, col=2, log="xy", xlab="mean gene expression", ylab=expression(paste(gamma, "+")), cex.axis = 1.5, cex.lab=1.7, cex=2 )
  abline(h=sum(g.dir.org[,1]), lty=2, lwd=2)
  
  dev.off()
  
  
  invisible(table.gh0)
  
} 


# g2fit <- g2fit.inits
# name <- paste0(name1, "inits_")
# plot.trend(g2fit, name)
# 
# 
# g2fit <- g2fit.init
# name <- paste0(name1, "init_")
# plot.trend(g2fit, name)


g2fit <- g2fit.obs
name <- paste0(name1, "obs_")
table.gh0 <- plot.trend(g2fit, name)
  


  

######################################
# run equalTheta
######################################

# sample.size <- 3
# s <- "2"
# 
# #name <- paste0("SIM",sample.size,"_s",s,"_init_NOTexact_")
# #name <- paste0("SIM",sample.size,"_s",s,"_initscalar_NOTexact_")
# 
# #name <- paste0("SIM",sample.size,"_s",s,"_init_exact_")
# #name <- paste0("SIM",sample.size,"_s",s,"_initscalar_exact_")
# 
# name <- paste0("SIM",sample.size,"_s",s,"_")
# 
# # simulate some dirichlets
# g.dir.org <- cbind(c(80, 10, 10) )
# g.dir <- rbind( rdirichlet( sample.size, g.dir.org[,1] ) )
# 
# sim <- NULL
# genes <- NULL
# nr.genes <- 50
# 
# for(i in 1:nr.genes){
#   # i=1
#   tot <- rnbinom(2*sample.size, mu=100, size=3)
#   #tot <- round(rnorm(sample.size, mean=200, sd=20))  
#   d <- sapply(1:(sample.size), function(u) rmultinom(1, prob=g.dir[u,], size=tot[u]))  
#   sim <- rbind(sim, d)  
#   genes <- c( genes, paste0("g", rep(i, 3))  )  
# }
# 
# rownames(sim) <- genes
# dge <- DGEList( counts=sim, group = c(rep(1, sample.size)), genes=data.frame(gene.id=genes, ete.id =genes) )


dge <- dge[, 1:3]

### run equalTheta

split_counts <- split(as.data.frame(dge$counts), dge$genes$gene.id)
split_counts2 <- lapply(split_counts, function(g){
  gt <- t(g)
  colnames(gt) <- paste0("I", 1:ncol(gt))
  return(gt)
})



est <- equalTheta(data=split_counts2, theta = 1/(1+sum(g.dir.org)), trace=FALSE)

names(est)

est.gamma <- mean(unlist(lapply(est$gamma, sum)))
est.gamma




################################################################################
# profile likelihood plot
################################################################################


# genes.up <- as.character(table.gh0[table.gh0$gh0.0f > 1e+5, "GeneID"])
# genes.down <- as.character(table.gh0[table.gh0$gh0.0f < 1e+3, "GeneID"])

genes.up <- as.character(table.gh0[table.gh0$gh0.0f > 1e+5, "GeneID"])
genes.down <- as.character(table.gh0[table.gh0$gh0.0f < 4, "GeneID"])

dir.par1 <- g.dir.org[, 1]

genes.plot <- c(genes.up[1:5], genes.down[1:5])

genes.plot <- c(genes.down[1:5])


#gamma.v <- c(seq(2, 100, 2), seq(120, 2000, 20), seq(2200, 100000, 200), seq(100000, 1000000, 10000))

#gamma.v <- c(seq(10, 100, 5), seq(100, 1000, 50), seq(1e3, 1e4, 5e2), seq(1e4, 1e5, 5e3), seq(1e5, 1e6, 1e5))

gamma.v <- c(seq(0.01, 2, 0.01), seq(2, 100, 2), seq(120, 2000, 20))

theta.v <- 1/(1+gamma.v)


pdf(paste0("PLOTS/", name, "PL_",gsub(":", "-" ,gsub(" ", "_", as.character(date()))),".pdf"))
par(mfrow=c(2,2))

for(j in genes.plot){
  # j=genes.plot[1]
  PLL.v <- rep(NA, length(theta.v))
  data <- split_counts2[j][[1]]
  
  for(i in 1:length(theta.v)){
    # i=700
    #cat(paste0(i, "\n"))

    try(PLL <-  estProfLogLik(data, theta=theta.v[i], trace=FALSE, maxit=1e3))
    if(is.null(PLL$loglik)) break
    PLL.v[i] <- PLL$loglik
    
  }
   
  theta.est <- table.gh0[table.gh0$GeneID==j, "gh0.0f"]
  
  plot(theta.v, PLL.v, type="l", main=paste0("gene " , j ,"\n n = ", sample.size, "\n dir par = ", paste(dir.par1, collapse=", ")), xlim=c(0, 1))
  abline(v = 1/(1+sum(theta.est)), col=2 )
  abline(v = 1/(1+sum(dir.par1)) , lty=2)
  
  
  plot(theta.v, PLL.v, type="l", main=paste0("gene " , j ,"\n n = ", sample.size, "\n dir par = ", paste(dir.par1, collapse=", ")), xlim=c(0.4, 0.6))
  abline(v = 1/(1+sum(theta.est)), col=2 )
  abline(v = 1/(1+sum(dir.par1)) , lty=2)
  
  
#   plot(theta.v, PLL.v, type="l", log="x", main=paste0("gene " , j ,"\n n = ", sample.size, "\n dir par = ", paste(dir.par1, collapse=", ")))
#   abline(v = 1/(1+sum(theta.est)), col=2 )
#   abline(v = 1/(1+sum(dir.par1)), lty=2)
#    
  
  #gp.e <- gridProf(data, theta=theta.est, from = - 1e-6, to = 1e-5, len=100)

}

dev.off()




################################################################################
# check what is wrong in estProfLogLik for selected genes 
################################################################################



dir.par1 <- g.dir.org[, 1]
genes.plot <- c("g10")



pdf(paste0("PLOTS/", name, "checkPL_",gsub(":", "-" ,gsub(" ", "_", as.character(date()))),".pdf"))

# gamma.v <- c(seq(1e4, 1e5, 1e3), seq(1e5, 1e6, 1e5))

gamma.v <- c(seq(1e4, 1e5, 1e2), seq(1e5, 1e6, 1e5))
gamma.v <- gamma.v[gamma.v <= 4e4]


theta.v <- 1/(1+gamma.v)

for(j in genes.plot){
  # j=genes.plot[1]
  PLL.v <- rep(NA, length(theta.v))
  data <- split_counts2[j][[1]]
  
  for(i in 1:length(theta.v)){
    # i=700
    cat(paste0(j, " ", theta.v[i]," ", gamma.v[i], "\n"))
    
    try(PLL <-  estProfLogLik(data, theta=theta.v[i], trace=FALSE, maxit=1e4))
    if(is.null(PLL$loglik)) PLL$loglik <- NA #break
    PLL.v[i] <- PLL$loglik
    
    cat(paste0(PLL$loglik), "\n")
    
  }
  
  theta.est <- table.gh0[table.gh0$GeneID==j, "gh0.0f"]
  
  plot(theta.v, PLL.v, type="l", main=paste0("gene " , j ,"\n n = ", sample.size, "\n dir par = ", paste(dir.par1, collapse=", ")))
  abline(v = 1/(1+sum(theta.est)), col=2 )
  abline(v = 1/(1+sum(dir.par1)) , lty=2)
  
  plot(theta.v, PLL.v, type="l", log="x", main=paste0("gene " , j ,"\n n = ", sample.size, "\n dir par = ", paste(dir.par1, collapse=", ")))
  abline(v = 1/(1+sum(theta.est)), col=2 )
  abline(v = 1/(1+sum(dir.par1)), lty=2)
  
  
}

dev.off()


#### check the code


dir.par1 <- g.dir.org[, 1]
genes.plot <- c("g10")
gamma.v <- 23000

theta.v <- 1/(1+gamma.v)


j = genes.plot[1]
PLL.v <- rep(NA, length(theta.v))
data <- split_counts2[j][[1]]

i = 1
cat(paste0(j, " ", theta.v[i]," ", gamma.v[i], "\n"))

PLL <-  estProfLogLik(data, theta=theta.v[i], trace=TRUE, maxit=1e5)


PLL.v[i] <- PLL$loglik

cat(paste0(PLL$loglik), "\n")














































