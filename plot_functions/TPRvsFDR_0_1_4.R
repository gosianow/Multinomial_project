#######################################################
# TPR vs acieved FDR --- my functions 
#######################################################

TPRvsFDR <- function(scores, colours, plotPath, xlim=c(0,0.5), ylim=c(0,1)){
	
	methods <- names(scores)
	
  pdf(plotPath)
	opar <- par()      # make a copy of current settings
	par(mar = c(5, 5, 4, 2) + 0.1, mgp = c(3.5, 1, 0)) # c(5, 4, 4, 2) + 0.1 # c(bottom, left, top, right)
	
  apvs <- scores[[1]][, "apvs"]
  status <- scores[[1]][, "status"]
  
  NAs <- is.na(status)
  apvs <- apvs[!NAs]
  status <- status[!NAs]
  
  
  TPRvsFDR_perMethod(status, apvs, col=colours[methods[1]], xlim=xlim, ylim=ylim, lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in 2:length(methods)){
	  apvs <- scores[[i]][, "apvs"]
	  status <- scores[[i]][, "status"]
    NAs <- is.na(status)
    apvs <- apvs[!NAs]
    status <- status[!NAs]
    
    TPRvsFDR_perMethod(status, apvs, col = colours[methods[i]], add=TRUE, lwd=4, cex=2.5)
  }
 
  
  legend("bottomright", methods, col = colours[methods], lty=1, lwd=4, cex=1, bty="o", bg="white", box.col="white")
  dev.off()
   	par(mar = opar$mar, mgp = opar$mgp)          # restore original settings
	
}







TPRvsFDR_perMethod <- function(status, apvs, FDR.cut.off=c(0.01, 0.05, 0.1), pch=c( 22, 23, 24), col="red", cex.axis=1.5, cex=2.5, add=FALSE, ...){
  
  apvs[is.na(apvs)] <- 1
  
  n <- length(status)
  q <- length(FDR.cut.off)
  TPR <- rep(0, q)
  FDR <- rep(0, q)
  
  for(i in 1:q){
    # i=1
    status.est <- as.numeric(apvs < FDR.cut.off[i])
    
    TP <- sum(status==1 & status.est==1)
    FP <- sum(status==0 & status.est==1)
    FN <- sum(status==1 & status.est==0)
    
    TPR[i] <- TP/(TP+FN)
    FDR[i] <- FP/(FP+TP)
    
  }  
  
  
  tf <- cbind( FDR , TPR)
  
  bg <- rep(col,q)
  bg[FDR > FDR.cut.off] <- "white"
  
  if(add==FALSE){
    
    plot(1:2, type="n", xlab="Achieved FDR", ylab="TPR", xaxt="n", col=col, cex.axis=cex.axis, las = 1, ...)
    axis(side=1, at=(2:10)/10, labels=(2:10)/10, las=1, col.ticks="grey", col.axis="grey", cex.axis=cex.axis)
    axis(side=1, at=FDR.cut.off, labels=FDR.cut.off, las=1, cex.axis=cex.axis)
    
    for(i in 1:q)
      lines(rep(FDR.cut.off[i], 50), seq(-0.1,1.1,length.out=50), type="b", pch=pch[i], cex=0.5, bg=1)
    
    lines(tf, type="l", col=col, ...)
    points(tf, pch=pch, bg=bg, col=col, cex=cex, ...)
    
  }else{    
    lines(tf, type="l", col=col, ...)
    points(tf, pch=pch, bg=bg, col=col, cex=cex,...)
  }
  
  
  
}



TPRvsFDRextend_perMethod <- function(status, apvs, FDR.cut.off=c(0.01, 0.05, 0.1), pch=c( 22, 23, 24), col="red", cex.axis=1.5, cex=2.5, add=FALSE, ...){
  
  apvs[is.na(apvs)] <- 1
  
  FDR.cut.off.org <- FDR.cut.off
  FDR.cut.off <- c(FDR.cut.off.org, seq(max(FDR.cut.off.org), 1, 0.01))
  
  n <- length(status)
  q.org <- length(FDR.cut.off.org)
  q <- length(FDR.cut.off)
  TPR <- rep(0, q)
  FDR <- rep(0, q)
  
  for(i in 1:q){
    # i=1
    status.est <- as.numeric(apvs < FDR.cut.off[i])
    
    TP <- sum(status==1 & status.est==1)
    FP <- sum(status==0 & status.est==1)
    FN <- sum(status==1 & status.est==0)
    
    TPR[i] <- TP/(TP+FN)
    FDR[i] <- FP/(FP+TP)
    
  }  
  tf <- cbind( FDR , TPR)
  
  bg <- rep(col,q.org)
  bg[FDR[1:q.org] > FDR.cut.off.org] <- "white"
  
  if(add==FALSE){
    
    plot(1:2, type="n", xlab="Achieved FDR", ylab="TPR", xaxt="n", col=col, cex.axis=cex.axis, ...)
    axis(side=1, at=(2:10)/10, labels=(2:10)/10, las=1, col.ticks="grey", col.axis="grey", cex.axis=cex.axis)
    axis(side=1, at=FDR.cut.off.org, labels=FDR.cut.off.org, las=1, cex.axis=cex.axis)
    
    for(i in 1:q.org)
      lines(rep(FDR.cut.off[i], 50), seq(-0.1,1.1,length.out=50), type="b", pch=pch[i], cex=0.5, bg=1)
    
    lines(tf, type="l", col=col, ...)
    points(tf[1:q.org,], pch=pch, bg=bg, col=col, cex=cex, ...)
    
  }else{    
    lines(tf, type="l", col=col, ...)
    points(tf[1:q.org,], pch=pch, bg=bg, col=col, cex=cex,...)
  }
  
  
  
}

