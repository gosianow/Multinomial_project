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
