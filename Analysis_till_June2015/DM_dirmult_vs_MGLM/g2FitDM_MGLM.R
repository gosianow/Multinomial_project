#######################################################
# g2 DM methods based on MGLM package 
# created 22July 2014 / updated 22 July 2014
#######################################################
# BioC 14


library(MGLM)
library(gtools)
library(edgeR)
library(parallel)


library(plyr)
rep.row <- function(r, n){
  r <- as.data.frame(t(as.matrix(r, 1, length(r))))
  colwise(function(x) rep(x, n))(r) 
}




#######################################################
# g2FitDM_MGLM
#######################################################


# group=dge$samples$group; init=NULL

g2FitDM_MGLM <- function(dge, group=dge$samples$group, init=NULL, mc.cores=15){
  
  
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

      null <- MGLMfit( Ygr1, dist="DM")
      
      loglikh <- null$logL
      gh <- rep.row(null$estimate, nrow(Y))
      gh0 <- rowSums(gh)
      th <- gh / gh0
          
      return(list(loglikh=loglikh, Y=Y, gh=gh, gh0= gh0 ,th=th, df=(q-1)))
      
    }, mc.cores=mc.cores)

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
      
      if (!is.null(init)) {
        g1 <- MGLMfit( Ygr1 , dist="DM", init=init)
        g2 <- MGLMfit( Ygr2 , dist="DM", init=init)        
      } else {
        g1 <- MGLMfit( Ygr1, dist="DM")
        g2 <- MGLMfit( Ygr2, dist="DM")
      }
      
      
      loglikh <- g1$logL + g2$logL
      gh <- rbind(rep.row(g1$estimate, nrow(Ygr1)), rep.row(g2$estimate, nrow(Ygr2)))
      gh0 <- rowSums(gh)
      th <- gh / gh0

      return( list( loglikh=loglikh, Y=Y, gh=gh, gh0= gh0 ,th=th, df=2*(q-1)) )
      
    }, mc.cores=mc.cores)
    
    
  }
  
  names(FitDM) <- gene.id.unique
  return(list(counts=dge$counts, samples=dge$samples, genes=dge$genes, group=group, FitDM=FitDM))


}



#######################################################
# 2gLRTDM
#######################################################




g2LRTDM <- function(g2fit, group=rep(1, length(dge$samples$group)), mc.cores=15){
  
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
    
  }, mc.cores=mc.cores)
  
  
  LRT <- do.call(rbind, LRT)
  FDR <- p.adjust(LRT$PValue, method="BH")
  o <- order(LRT$PValue)
  
  table <- data.frame(GeneID=names(fit.full$FitDM), LR=LRT$LR, df=LRT$df, LLfull=LRT$LLfull, LLnull=LRT$LLnull , PValue=LRT$PValue, FDR=FDR)[o,]
  
  return(list(counts=dge$counts, samples=dge$samples, genes=dge$genes, group=fit.full$group, group0=fit.null$group, fit.full=fit.full$FitDM, fit.null=fit.null$FitDM,  table=table))
  
  
}


















