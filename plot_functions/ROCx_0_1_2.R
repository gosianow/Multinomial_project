#######################################################
# ROCx plot  --- my function
#######################################################

## with cumsum 


ROCx <- function(pvs, apvs, status, add=FALSE, ...){
  # status==1, means true DS, status!=1, means not DS
  
  status.org <- status
  P.org <- sum(status.org==1)
  N.org <- sum(status.org==0)
  
  NAs <- is.na(pvs)
  pvs <- pvs[!NAs]
  status <- status[!NAs]
  
  P <- sum(status==1)
  N <- sum(status==0)
  
  ord <- order(pvs, decreasing = FALSE, na.last = TRUE)
  status <- status[ord]
  pvs <- pvs[ord]
  
  TPRv <- cumsum(status) / P
  FPRv <- cumsum(!status) / N
  
  
  if(add==FALSE){
    plot(FPRv, TPRv, xlab="False positive rate", ylab="True positive rate", type="l", ...)
  }else{   
    lines(FPRv, TPRv, type="l", ...)
  }
  
  TPR <- sum(apvs[status.org==1] < 0.05 , na.rm = T) / P.org
  points((approx(TPRv, FPRv, xout=TPR)$y), TPR, ...)  
  
}




ROCx_notNorm <- function(pvs, apvs, status, add=FALSE, ...){
  # status==1, means true DS, status!=1, means not DS
  
  status.org <- status
  P <- sum(status.org==1)
  N <- sum(status.org==0)
  
  NAs <- is.na(pvs)
  pvs <- pvs[!NAs]
  status <- status[!NAs]
  
  ord <- order(pvs, decreasing = FALSE, na.last = TRUE)
  status <- status[ord]
  pvs <- pvs[ord]
  
  TPRv <- cumsum(status) / P
  FPRv <- cumsum(!status) / N
  
  
  if(add==FALSE){
    plot(FPRv, TPRv, xlab="False positive rate", ylab="True positive rate", type="l", ...)
  }else{   
    lines(FPRv, TPRv, type="l", ...)
  }
  
  TPR <- sum(apvs[status.org==1] < 0.05 , na.rm = T) / P
  points((approx(TPRv, FPRv, xout=TPR)$y), TPR, ...)  
  
}

library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)



ggROCx_notNorm <- function(scores, colours){

  methods <- names(scores)
  Xlist <- list()
  ROClist <- list()
  
  for(m in seq(length(scores))){
    # m = 1
    
    sc <- scores[[m]]
    
    NAs <- apply(sc[, "status", drop = FALSE], 1, function(r){any(is.na(r)) } )   
    sc <- sc[!NAs, ]
    
    pvs <- sc$pvs
    apvs <- sc$apvs
    status <- sc$status
    
    table(status, useNA = "always")
    
    status.org <- status
    P <- sum(status.org==1)
    N <- sum(status.org==0)
    
    NAs <- is.na(pvs)
    pvs <- pvs[!NAs]
    status <- status[!NAs]
    
    ord <- order(pvs, decreasing = FALSE, na.last = TRUE)
    status <- status[ord]
    pvs <- pvs[ord]
    
    TPRv <- cumsum(status) / P
    FPRv <- cumsum(!status) / N
    ROClist[[m]] <- data.frame(FPR = FPRv, TPR = TPRv, Method = methods[m])
    
    ### what is the TPR if the error was controlled 
    TPR <- sum(apvs[status.org==1] < 0.05 , na.rm = T) / P    
    Xlist[[m]] <- data.frame(FPR = approx(TPRv, FPRv, xout=TPR)$y, TPR = TPR, Method = methods[m])    
  
  }
  
  X <- do.call(rbind, Xlist)
  ROC <- do.call(rbind, ROClist)
  
  ggp <- ggplot(data = ROC, aes(x = FPR, y = TPR, group = Method, colour = Method)) +
    theme_bw() +
    geom_line(size=2) +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.justification=c(1,0), legend.position=c(1,0)) +
    guides(colour = guide_legend(override.aes = list(size = 2, shape = NA))) +
    scale_color_manual(values = colours[levels(ROC$Method)]) +
    geom_point(data = X, aes(x = FPR, y = TPR, group = Method, colour = Method), size = 10, shape = "X") 

  
  # ggsave(filename = paste0(out.dir, "/test2.pdf"), width = 10, height = 7)
  
print(ggp)
}

































