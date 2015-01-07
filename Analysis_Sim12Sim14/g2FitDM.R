#######################################################
# g2 method with estProfLogLik
#######################################################




#######################################################
# 2gLRTDM
#######################################################

library("dirmult")
library("gtools")


library(plyr)
rep.row <- function(r, n){
  r <- as.data.frame(t(as.matrix(r, 1, length(r))))
  colwise(function(x) rep(x, n))(r) 
}



g2FitDM <- function(dge, group=dge$samples$group, initscalar=NULL, init=NULL){
  
  
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
      
      null <- dirmult( Y , trace=FALSE)
      
      loglikh <- null$loglik
      gh <- rep.row(null$gamma, nrow(Y))
      gh0 <- rowSums(gh)
      th <- gh / gh0
      
      mom <- null$mom
      
      return(list(loglikh=loglikh, Y=Y, gh=gh, gh0= gh0 ,th=th, df=(q-1) , mom=mom))
      
    }, mc.cores=20)
    
  } else {
    
    levels(group) <- c("g1", "g2")
    
    
    Y.means <- mclapply(seq_len(length(gene.id.unique)), function(g){
      # g=1
      
      counts.tmp <- dge$counts[split_counts[[g]], , drop=FALSE]
      
      if(nrow(counts.tmp)==1)
        return(NULL)
      
      Y <- t(counts.tmp)
      Ygr1 <- Y[group=="g1",]
      Ygr2 <- Y[group=="g2",]
      
      
    } , mc.cores=20) 
    
    
    
    
    
    
    
    FitDM <- mclapply(seq_len(length(gene.id.unique)), function(g){
      # g=1
      
      counts.tmp <- dge$counts[split_counts[[g]], , drop=FALSE]
      
      if(nrow(counts.tmp)==1)
        return(NULL)
      
      Y <- t(counts.tmp)
      q <- ncol(Y)
      
      Ygr1 <- Y[group=="g1",]
      Ygr2 <- Y[group=="g2",]
      
      if (!is.null(initscalar)) {
        g1 <- dirmult( Ygr1 , trace=FALSE, initscalar=initscalar[1])
        g2 <- dirmult( Ygr2 , trace=FALSE, initscalar=initscalar[2])        
      } else if (!is.null(init)) {
        g1 <- dirmult( Ygr1 , trace=FALSE, init=init[,1])
        g2 <- dirmult( Ygr2 , trace=FALSE, init=init[,2])       
      } else {
        g1 <- dirmult( Ygr1 , trace=FALSE)
        g2 <- dirmult( Ygr2 , trace=FALSE)
      }
      
      
      loglikh <- g1$loglik + g2$loglik
      gh <- rbind(rep.row(g1$gamma, nrow(Ygr1)), rep.row(g2$gamma, nrow(Ygr2)))
      gh0 <- rowSums(gh)
      th <- gh / gh0
      
      mom <- c(g1$mom, g2$mom)
      
      
      # ProfLogLik
      
      
      PPL <- estProfLogLik(data, theta)
      
      
      
      return(list(loglikh=loglikh, Y=Y, gh=gh, gh0= gh0 ,th=th, df=2*(q-1), mom=mom ))
      
    }, mc.cores=20)
    
  }
  
  names(FitDM) <- gene.id.unique
  return(list(counts=dge$counts, samples=dge$samples, genes=dge$genes, group=group, FitDM=FitDM))
  
  
}





g2LRTDM <- function(g2fit, group=rep(1, length(dge$samples$group))){
  
  fit.full <- g2fit
  fit.null <- g2FitDM(g2fit, group=group)
  
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
  
  return(list(counts=dge$counts, samples=dge$samples, genes=dge$genes, group=fit.full$group, group0=fit.null$group, fit.full=fit.full$FitDM, fit.null=fit.null$FitDM,  table=table))
  
  
}


















