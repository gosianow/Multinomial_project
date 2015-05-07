##############################################################################
# Creaed 12 Jun 2014
# smulate data from known DIR-MULTI 
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

sample.size <- 3
s <- "1"

#name <- paste0("SIM",sample.size,"_s",s,"_init_NOTexact_")
#name <- paste0("SIM",sample.size,"_s",s,"_initscalar_NOTexact_")

#name <- paste0("SIM",sample.size,"_s",s,"_init_exact_")
#name <- paste0("SIM",sample.size,"_s",s,"_initscalar_exact_")

name1 <- paste0("SIM",sample.size,"_s",s,"_")


# simulate some dirichlets, theta ~ 0.03
g.dir.org <- cbind(c(26, 2, 2), c(10, 10, 10) ) # s <- "1"

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


# load(paste0("PLOTS/","SIM3_s1_dge_Wed_Jun_11_15-51-33_2014.RData"))


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
  




################################################################################
# check what it the difference in dirmult for genes from 2 clouds
################################################################################




genes.up <- as.character(table.gh0[table.gh0$gh0.0f > 1e+5, "GeneID"])
genes.down <- as.character(table.gh0[table.gh0$gh0.0f < 1e+2, "GeneID"])


dir.par1 <- g.dir.org[, 1]




group=dge$samples$group; initscalar=NULL; init=NULL; mode="obs"



split_counts <- split(seq_len(nrow(dge$counts)), dge$genes$gene.id, drop=TRUE)
gene.id.unique <- names(split_counts)


group <- as.factor(group)



levels(group) <- c("g1", "g2")

g <- genes.up[1]
#g <- genes.down[2]


cat("gene ", g, "\n")

counts.tmp <- dge$counts[split_counts[[g]], , drop=FALSE]


Y <- t(counts.tmp)
Y
q <- ncol(Y)


Ygr1 <- Y[group=="g1",]
Ygr1
Ygr2 <- Y[group=="g2",]
Ygr2


g1 <- dirmult( Ygr1 , trace=FALSE, mode=mode)
g2 <- dirmult( Ygr2 , trace=FALSE, mode=mode)

g1


dirmult( Ygr1 , trace=TRUE, mode=mode, epsilon=10^(-6))

dirmult( Ygr1 , trace=TRUE, mode=mode, initscalar=30)



g.dir.org/30



data=Ygr1; epsilon=10^(-4); trace=FALSE; mode=mode


## Estimate parameters in the Dirichlet-Multinomial distribution
dirmult2 <- function(data,init,initscalar,epsilon=10^(-4),trace=TRUE,mode){
  data <- data[rowSums(data)!=0,colSums(data)!=0]
  
  if(missing(initscalar)){
    mom <- weirMoM(data) # theta
    if(mom<=0) mom <- 0.005
    initscalar <- (1-mom)/mom
  }else{
    mom <- 1/(1+initscalar)
  }
  if(missing(init)){
    gamma <- colSums(data)/sum(data)*initscalar
  }else{
    gamma <- init
    initscalar <- sum(gamma)
    mom <- 1/(1+initscalar)
  } 
  
  if(missing(mode)) mode <- "obs"
  if(!is.element(mode,c("obs","exp"))){
    message(paste("Warning: Mode '",mode,"' not valid\n",sep=""))
    mode <- "obs"
  }
  lik1 <- 0
  lik2 <- epsilon*10
  ite <- 1
  gamite <- 0
  conv <- TRUE
  # Iterations
  while(conv){
    
    if(abs(lik2-lik1) < epsilon) conv <- FALSE
    
    
    if(mode=="exp") fim <- expfim(data,gamma)
    else if(mode=="obs") fim <- -1*obsfim(data,gamma)
    
    lik1 <- loglik(data,gamma)
    
    # Updates parameter estimates
    
    update <- solve(fim) %*% u(data,gamma)
    print(update)
    gamma <- gamma + update
    print(gamma)
    
    gamma[gamma<0] <- 0.01 # Negative gamma_j are set to 0.01
    if(any(gamma<0) & (gamite%%10)==0){ print(gamma); gamite <- gamite+1}
    if(trace) message(paste("Iteration ",ite,": Log-likelihood value: ",lik1,sep=""))
    gams <- paste(" Gamma",1:length(gamma),sep="")
    lik2 <- loglik(data,gamma)
    ite <- ite+1
    
  }
  sumgam <- sum(gamma)
  theta <- 1/(sumgam+1)
  pi <- as.numeric(gamma/sumgam)
  names(pi) <- dimnames(data)[[2]]
  list(loglik=lik1,ite=ite-1,gamma=as.numeric(gamma), gamma0=sum(as.numeric(gamma)) ,pi=pi,theta=theta, theta.mom=mom, gamma.mom=(1-mom)/mom)
}



dirmult2( Ygr1 , trace=TRUE, mode=mode)


dirmult2( Ygr1 , trace=TRUE, mode=mode, initscalar=778968)






m <- fim 

im1 <- solve(m)

im2 <- qr.solve(m)

all.equal(im1, im2)































