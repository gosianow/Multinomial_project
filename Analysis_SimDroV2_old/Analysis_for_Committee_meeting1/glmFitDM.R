#######################################################
# glm and g2 methods 
# updated 17 July 2014
#######################################################



library("edgeR")
# library("DEXSeq")
library(parallel)

#######################################################
# Loglik from DirMulti
#######################################################

Loglik <- function( b, Y, X) {
  
  p <- ncol(X)
  q <- ncol(Y)
  b <- matrix(b, q, p)
  
  g <- exp(X %*% t(b))  # n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)
  
  # formula (7) is equivalent to likelihood from dirmult
  res <-   sum(lgamma(gs) - lgamma(ys+gs) + rowSums(lgamma(Y+g) - lgamma(g)))
  # normalization coeff
  # norm <- sum(lgamma(ys + 1) - rowSums(lgamma(Y + 1)))
  # res <- res + norm
  
  -res
  
}


#######################################################
# Loglik from dirmult
#######################################################


ll <- function (x, t)
  log(t + x - 1)


loglik <- function (x, t)
{
  
  
  l <- 0
  ts <- sum(t)
  nc <- ncol(x)
  for (j in 1:nrow(x)) {
    
    l <- l - sum(unlist(lapply(list(1:(rowSums(x)[j])), ll, t = ts)))
    for (i in 1:nc) {
      if (x[j, i] == 0)
        lij <- 0
      else lij <- sum(unlist(lapply(list(1:(x[j, i])),ll, t = t[i])))
      l <- l + lij
    }
    
    
  }
  l
}


Loglik2 <- function( b, Y, X){
  
  p <- ncol(X)
  q <- ncol(Y)
  b <- matrix(b, q, p)
  
  g <- exp(X %*% t(b))  # n * q
  
  res <- sum(unlist(lapply( 1:nrow(Y), function(i){
    
    loglik(x=Y[i, , drop=FALSE], t=g[i, , drop=FALSE])
    
  } )))
  
  
  -res
    
} 




#######################################################
# glmFitDM
#######################################################

glmFitDM <- function(dge, design, optimization="optim"){
  
  
  split_counts <- split(seq_len(nrow(dge$counts)), dge$genes$gene.id, drop=TRUE)
  gene.id.unique <- names(split_counts)
  
  X <- as.matrix(design)
  p <- ncol(X)
  
  FitDM <- mclapply(seq_len(length(gene.id.unique)), function(g){
    # g=1
    
    counts.tmp <- dge$counts[split_counts[[g]], , drop=FALSE]
    
    if(nrow(counts.tmp)==1)
      return(NULL)
    
    Y <- t(counts.tmp)
    q <- ncol(Y)
    b <- rep(0, q*p)
    
    if(optimization=="optim"){
      
      bh_opt <- optim(par=b, fn=Loglik, Y=Y, X=X, method="Nelder-Mead")
      
      loglikh <- -bh_opt$value
      bh <- matrix(bh_opt$par, q, p)
      
    }else if(optimization=="nlm"){
      
      bh_nlm <- nlm(Loglik, b, Y=Y, X=X)
      
      loglikh <- -bh_nlm$minimum
      bh <- matrix(bh_nlm$estimate, q, p)
      
    }
    
    gh <- exp(X %*% t(bh))
    gh0 <- rowSums(gh)
    th <- gh / gh0
    
    return(list(loglikh=loglikh, Y=Y, bh=bh, gh=gh, gh0= gh0 ,th=th, df=p*(q-1) ))
    
  }, mc.cores=20)
  
  names(FitDM) <- gene.id.unique
  
  return(list(counts=dge$counts, samples=dge$samples, genes=dge$genes,  design=design, optimization=optimization, FitDM=FitDM))
  
}


#######################################################
# glmLRTDM
#######################################################



glmLRTDM <- function(glmfit, coef=ncol(glmfit$design)){
  
  fit.full <- glmfit
  
  design0 <- fit.full$design[, -coef, drop = FALSE]
  
  fit.null <- glmFitDM(fit.full, design=design0, optimization=fit.full$optimization)
   
  LRT <- mclapply(names(fit.full$FitDM), function(g){
    
    if(!is.null(fit.full$FitDM[[g]])){           
      LR <-  2*(fit.full$FitDM[[g]]$loglikh - fit.null$FitDM[[g]]$loglikh)
      df <- fit.full$FitDM[[g]]$df - fit.null$FitDM[[g]]$df
      pvalue <- pchisq(LR, df = df , lower.tail = FALSE)
      return(data.frame(LR=LR, df=df, PValue=pvalue, LLfull=fit.full$FitDM[[g]]$loglikh, LLnull=fit.null$FitDM[[g]]$loglikh))
    }    
    return(data.frame(LR=NA, df=NA, PValue=NA, LLfull=NA, LLnull=NA))
    
  }, mc.cores=10)
  

  LRT <- do.call(rbind, LRT)
  FDR <- p.adjust(LRT$PValue, method="BH")
  o <- order(LRT$PValue)
  
  table <- data.frame(GeneID=names(fit.full$FitDM), LR=LRT$LR, df=LRT$df, LLfull=LRT$LLfull, LLnull=LRT$LLnull , PValue=LRT$PValue, FDR=FDR)[o,]
  
  return(list(counts=dge$counts, samples=dge$samples, genes=dge$genes,  design=fit.full$design, design0=design0, optimization=fit.full$optimization, fit.full=fit.full$FitDM, fit.null=fit.null$FitDM,  table=table))
  
  
}




#######################################################
# g2FitDM
#######################################################

library("dirmult")
library("gtools")


library(plyr)
rep.row <- function(r, n){
  r <- as.data.frame(t(as.matrix(r, 1, length(r))))
  colwise(function(x) rep(x, n))(r) 
}

# group=dge$samples$group; initscalar=NULL; init=NULL

g2FitDM <- function(dge, group=dge$samples$group, initscalar=NULL, init=NULL, mode="obs"){
  
  
  split_counts <- split(seq_len(nrow(dge$counts)), dge$genes$gene.id, drop=TRUE)
  gene.id.unique <- names(split_counts)
  

  group <- as.factor(group)

  
  if(length(levels(group))==1){
    
    
    FitDM <- mclapply(seq_len(length(gene.id.unique)), function(g){
      # g=1
      
      counts.tmp <- dge$counts[split_counts[[g]], , drop=FALSE]
      
      if(nrow(counts.tmp)==1)
        return(NULL)
      
      Y <- t(counts.tmp)
      q <- ncol(Y)

      null <- dirmult( Y , trace=FALSE, mode=mode)
      
      loglikh <- null$loglik
      gh <- rep.row(null$gamma, nrow(Y))
      gh0 <- rowSums(gh)
      th <- gh / gh0
      
      mom <- null$mom
      
      return(list(loglikh=loglikh, Y=Y, gh=gh, gh0= gh0 ,th=th, df=(q-1) , mom=mom))
      
    }, mc.cores=20)

  } else {
    
    levels(group) <- c("g1", "g2")
    
    FitDM <- mclapply(seq_len(length(gene.id.unique)), function(g){
      # g=2
      
      #cat("gene ", g, "\n")
      
      counts.tmp <- dge$counts[split_counts[[g]], , drop=FALSE]
      
      if(nrow(counts.tmp)==1)
        return(NULL)
      
      Y <- t(counts.tmp)
      q <- ncol(Y)
      
      Ygr1 <- Y[group=="g1",]
      Ygr2 <- Y[group=="g2",]
      
      if (!is.null(initscalar)) {
        g1 <- dirmult( Ygr1 , trace=FALSE, initscalar=initscalar[1], mode=mode)
        g2 <- dirmult( Ygr2 , trace=FALSE, initscalar=initscalar[2], mode=mode)        
      } else if (!is.null(init)) {
        g1 <- dirmult( Ygr1 , trace=FALSE, init=init[,1], mode=mode)
        g2 <- dirmult( Ygr2 , trace=FALSE, init=init[,2], mode=mode)       
      } else {
        g1 <- dirmult( Ygr1 , trace=FALSE, mode=mode)
        g2 <- dirmult( Ygr2 , trace=FALSE, mode=mode)
      }
      
      
      loglikh <- g1$loglik + g2$loglik
      gh <- rbind(rep.row(g1$gamma, nrow(Ygr1)), rep.row(g2$gamma, nrow(Ygr2)))
      gh0 <- rowSums(gh)
      th <- gh / gh0
      
      mom <- c(g1$mom, g2$mom)      

      return( list( loglikh=loglikh, Y=Y, gh=gh, gh0= gh0 ,th=th, df=2*(q-1), mom=mom ) )
      
    }, mc.cores=10)
    
    
  }
  
  names(FitDM) <- gene.id.unique
  return(list(counts=dge$counts, samples=dge$samples, genes=dge$genes, group=group, FitDM=FitDM))


}



#######################################################
# 2gLRTDM
#######################################################




g2LRTDM <- function(g2fit, group=rep(1, length(dge$samples$group))){
  
  fit.full <- g2fit
  fit.null <- g2FitDM(g2fit, group=group)
  
  LRT <- mclapply(names(fit.full$FitDM), function(g){
    # g = "FBgn0000042"
    
    if( !is.null( fit.full$FitDM[[g]] ) ){           
      LR <-  2*(fit.full$FitDM[[g]]$loglikh - fit.null$FitDM[[g]]$loglikh)
      df <- fit.full$FitDM[[g]]$df - fit.null$FitDM[[g]]$df
      pvalue <- pchisq(LR, df = df , lower.tail = FALSE)
      return(data.frame(LR=LR, df=df, PValue=pvalue, LLfull=fit.full$FitDM[[g]]$loglikh, LLnull=fit.null$FitDM[[g]]$loglikh))
    }    
    return(data.frame(LR=NA, df=NA, PValue=NA, LLfull=NA, LLnull=NA))
    
  }, mc.cores=15)
  
  
  LRT <- do.call(rbind, LRT)
  FDR <- p.adjust(LRT$PValue, method="BH")
  o <- order(LRT$PValue)
  
  table <- data.frame(GeneID=names(fit.full$FitDM), LR=LRT$LR, df=LRT$df, LLfull=LRT$LLfull, LLnull=LRT$LLnull , PValue=LRT$PValue, FDR=FDR)[o,]
  
  return(list(counts=dge$counts, samples=dge$samples, genes=dge$genes, group=fit.full$group, group0=fit.null$group, fit.full=fit.full$FitDM, fit.null=fit.null$FitDM,  table=table))
  
  
}


















