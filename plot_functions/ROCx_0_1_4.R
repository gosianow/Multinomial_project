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



ROCx_notNorm <- function(scores, colours, plotPath, xlim=c(0,1), ylim=c(0,1)){
	
	methods <- names(scores)
	
  pdf(plotPath)
	opar <- par()      # make a copy of current settings
	par(mar = c(5, 5, 4, 2) + 0.1, mgp = c(3.5, 1, 0)) # c(5, 4, 4, 2) + 0.1 # c(bottom, left, top, right)
		
	pvs <- scores[[1]][, "pvs"]
  apvs <- scores[[1]][, "apvs"]
  status <- scores[[1]][, "status"]
  
  NAs <- is.na(status)
  apvs <- apvs[!NAs]
  status <- status[!NAs]
  
  
  ROCx_notNorm_perMethod(pvs, apvs, status, col=colours[methods[1]], xlim=xlim, ylim=ylim, lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in 2:length(methods)){
		pvs <- scores[[i]][, "pvs"]
	  apvs <- scores[[i]][, "apvs"]
	  status <- scores[[i]][, "status"]
    NAs <- is.na(status)
    apvs <- apvs[!NAs]
    status <- status[!NAs]
    
    ROCx_notNorm_perMethod(pvs, apvs, status, col = colours[methods[i]], add=TRUE, lwd=4, cex=2.5)
  }
 
  
  legend("bottomright", methods, col = colours[methods], lty=1, lwd=4, cex=1, bty="o", bg="white", box.col="white")
  dev.off()
       	par(mar = opar$mar, mgp = opar$mgp)          # restore original settings
	
	
}


ROCx_notNorm_perMethod <- function(pvs, apvs, status, add=FALSE, ...){
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
    plot(FPRv, TPRv, xlab="False positive rate", ylab="True positive rate", type="l", las = 1, ...)
  }else{   
    lines(FPRv, TPRv, type="l", ...)
  }

  TPR <- sum(apvs[status.org==1] < 0.05 , na.rm = T) / P
  points((approx(TPRv, FPRv, xout=TPR)$y), TPR, pch = "x", ...)  

}








library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)



ggROCx_notNorm <- function(scores, colours, plotPath){

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
  pdf(plotPath, width = 7, height = 7)
 
print(ggp)

  dev.off()
}

































