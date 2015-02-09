##############################################################################

# BioC 3.0
# Created 4 Feb 2015:

# Investigate problem in _SimDroV1_run_DM_0.1.2.R with TG-grid common dispersion when the filtering is relaxed

#

##############################################################################


table(dgeDM$tagwiseDispersion, useNA = "always")


colMeans(dgeDM$plotLoglik0)

colMedians(dgeDM$plotLoglik0)

apply(dgeDM$plotLoglik0, 2, median)



pdf(paste0(out.dir, "/",name1,"_DispGrid_CommonHist.pdf"))

df <- data.frame(genes = rownames(dge$plotLoglik0), dge$plotLoglik0)
df.m <- melt(df, id.vars = "genes", variable.name = "SplinePts", value.name = "Loglikelihood")
require(ggplot2)
ggp <- ggplot(data = df.m, aes(x = SplinePts, y = Loglikelihood, fill = SplinePts)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=18, size=4) +
  guides(fill=FALSE)

print(ggp)

dev.off()



df[which.min(df$X1), ]

dgeDM$counts[["FBgn0053196"]]


#### if plot == TRUE

plotGenes <- names(dgeDM$meanExpr)[which(dgeDM$meanExpr > 10 & dgeDM$meanExpr < 100 & dgeDMnone$tagwiseDispersion > 1e5)]
subset <- 10


pdf(paste0(out.dir, "/",name1,"_DispGrid_Likelihoods.pdf"))

for(i in plotGenes[seq(subset)]){
  
  
  
  ylim <- c(min(c(dgeDM$plotLoglik[plotGenes[i],], dgeDM$plotLoglik0[plotGenes[i],], dgeDM$plotModeration[plotGenes[i],])) ,max(c(dgeDM$plotLoglik[plotGenes[i],], dgeDM$plotLoglik0[plotGenes[i],], dgeDM$plotModeration[plotGenes[i],])))
  
  #   plot(log10(dgeDM$plotSplineDisp), log10(dgeDM$plotLoglik[plotGenes[i],] - ylim[1] + 1), type="l", lwd=3, col = "black", ylim = log10(ylim - ylim[1] +1) )
  #   lines( log10(dgeDM$plotSplineDisp), log10(dgeDM$plotLoglik0[plotGenes[i],]- ylim[1] + 1), type="l", lwd=3, col = "black", lty=2)
  #   lines( log10(dgeDM$plotSplineDisp), log10(dgeDM$plotModeration[plotGenes[i],]- ylim[1] + 1), type="l", lwd=3, col = "grey", lty=2)
  #   abline(v = log10(dgeDM$plotOutDisp[plotGenes[i]]), col = "pink")
  #   abline(v = log10(dgeDM$commonDispersion), col="grey", lty=2)
  
  r <- 1:15
  
  plot(dgeDM$plotSplineDisp[r], dgeDM$plotLoglik[plotGenes[i],r], type="l", lwd=3, col = "black", ylim=ylim )
  lines( dgeDM$plotSplineDisp[r], dgeDM$plotLoglik0[plotGenes[i],r], type="l", lwd=3, col = "black", lty=2)
  lines( dgeDM$plotSplineDisp[r], dgeDM$plotModeration[plotGenes[i],r], type="l", lwd=3, col = "grey", lty=2)
  abline(v = dgeDM$plotOutDisp[plotGenes[i]], col = "pink")
  abline(v = dgeDM$commonDispersion, col="grey", lty=2)
  abline( v = dgeDNnone$tagwiseDispersion[] )
  
  plot(dgeDM$plotSplineDisp[r], dgeDM$plotLoglik[plotGenes[i],r], type="l", lwd=3, col = "black", )
  abline(v = dgeDM$plotOutDisp[plotGenes[i]], col = "pink")
  abline(v = dgeDM$commonDispersion, col="grey", lty=2)
  
  plot( dgeDM$plotSplineDisp[r], dgeDM$plotLoglik0[plotGenes[i],r], type="l", lwd=3, col = "black", lty=2)
  abline(v = dgeDM$plotOutDisp[plotGenes[i]], col = "pink")
  abline(v = dgeDM$commonDispersion, col="grey", lty=2)
  
  plot( dgeDM$plotSplineDisp[r], dgeDM$plotModeration[plotGenes[i],r], type="l", lwd=3, col = "grey", lty=2)
  abline(v = dgeDM$plotOutDisp[plotGenes[i]], col = "pink")
  abline(v = dgeDM$commonDispersion, col="grey", lty=2)
}

dev.off()







