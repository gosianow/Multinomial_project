##############################################################################
# smulate data from known DIR-MULTI 
# compare g2FitDM dirmulti & MGLM
##############################################################################


library(edgeR)

setwd("/home/gosia/Multinomial_project/DM_dirmult_vs_MGLM/")

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

save(dge, file=paste0(name1, "dge" ,".RData"))



######################################
# fit DM dirmult
######################################

source("/home/gosia/R/R_Multinomial_project/DM_dirmult_vs_MGLM/g2FitDM_dirmult.R")

# define initial values
init <- g.dir.org - 1
initscalar <- colSums(g.dir.org) - 5

# run fitting
g2fit.dirmult <- g2FitDM(dge)

# g2fit.init <- g2FitDM(dge, init=init)
# g2fit.inits <- g2FitDM(dge, initscalar=initscalar)


######################################
# fit DM MGLM
######################################


source("/home/gosia/R/R_Multinomial_project/DM_dirmult_vs_MGLM/g2FitDM_MGLM.R")

# define initial values
init <- g.dir.org - 1

# run fitting
g2fit.mglm <- g2FitDM_MGLM(dge)




######################################
# plot the dispersion estimates
######################################
library(beanplot)

plot.trend <- function(g2fit, name, mc.cores=15){
  
  gh0.tmp <-  mclapply(names(g2fit$FitDM), function(g){
    
    if(!is.null(g2fit$FitDM[[g]])){           
      
      meanY0 <- mean(rowSums(g2fit$FitDM[[g]]$Y[1:sample.size,]))
      meanY1 <- mean(rowSums(g2fit$FitDM[[g]]$Y[(sample.size+1):(2*sample.size),]))
      gh0.0f <- mean(g2fit$FitDM[[g]]$gh0[1])
      gh0.1f <- mean(g2fit$FitDM[[g]]$gh0[sample.size+1])
      
      return(data.frame(gh0.0f=gh0.0f, gh0.1f=gh0.1f,  meanY0=meanY0, meanY1=meanY1))
    }    
    return(data.frame(gh0.0f=NA, gh0.1f=NA, meanY0=NA, meanY1=NA))
    
  }, mc.cores=mc.cores)
  
  
  gh0 <- do.call(rbind, gh0.tmp)
  
  gh0 <- data.frame(GeneID=names(g2fit$FitDM), gh0)
  head(gh0)
  table.gh0 <- gh0
  

  pdf(paste0(name, "ESTIMATION_gh0",".pdf"), h=5, w=10)
  par(mfrow=c(1,2 ), cex.main=1.1, cex.axis = 1.1, cex=1.2)
  
  
  plot(table.gh0$meanY0, table.gh0$gh0.0f, col=2, log="xy", main=paste0("gamma: ", paste(g.dir.org[,1], collapse=",")))
  abline(h=sum(g.dir.org[,1]), lty=2)
  plot(table.gh0$meanY1 , table.gh0$gh0.1f, col=3, log="xy", main=paste0("gamma: ", paste(g.dir.org[,2], collapse=",")))
  abline(h=sum(g.dir.org[,2]), lty=2)
  
  beanplot(table.gh0$gh0.0f, col=2, log="y", what=c(1,1,1,0), main=paste0("gamma: ", paste(g.dir.org[,1], collapse=",")))
  abline(h=sum(g.dir.org[,1]), lty=2)
  beanplot(table.gh0$gh0.1f, col=3, log="y", what=c(1,1,1,0), main=paste0("gamma: ", paste(g.dir.org[,2], collapse=",")))
  abline(h=sum(g.dir.org[,2]), lty=2)
  
  dev.off()
  
  invisible(table.gh0)
  
} 



tgh0d <- plot.trend(g2fit.dirmult, paste0(name1, "dirmult"))


tgh0m <- plot.trend(g2fit.mglm, paste0(name1, "mglm"))


tgh0 <- merge(tgh0d, tgh0m, by="GeneID")




pdf(paste0("Compare_dirmult_mglm",".pdf"))

smoothScatter(log(tgh0$gh0.0f.x), log(tgh0$gh0.0f.y))
abline(a=0, b=1, col=2)
smoothScatter(log(tgh0$gh0.1f.x), log(tgh0$gh0.1f.y))
abline(a=0, b=1, col=2)



dev.off()

























