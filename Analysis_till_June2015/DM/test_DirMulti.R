#######################################################
# run example
#######################################################


# install.packages("dirmult", dependencies=TRUE, lib="/home/gosia/R/libraries/3.0.2")


source("/home/gosia/R/R_Multinomial_project/DirMulti.R")


sim.obj <- DirmultSim()

# n=100 - samples, q=40 - spieces
Yex <- sim.obj$Y
head(Yex)

# n=100 - samples, p=101 - covariates
# first column of 1s
Xex <- sim.obj$X
head(Xex)


dmG.obj <- DirmultGrpGrid(Yex, Xex, model="dirmult")

names(dmG.obj[[1]])


dir.obj <- DirmultGrpGrid(Yex, Xex, model="dir")
mult.obj <- DirmultGrpGrid(Yex, Xex, model="mult")



#######################################################
# run on simulated data
#######################################################


setwd("/home/gosia/Multinomial_project/Simulations")

# load information about simulation 
simu_info <- read.table("Simu_info/simu12_reads_final.txt", header=T)
# write.table(simu_info, "Simu_info/simu12_reads_final.xls", sep="\t", quote=F, row.names=F, col.names=T)
head(simu_info)


# create metadata file
metadata <- data.frame(SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "R",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))


library("edgeR")
library("DEXSeq")

# run read.HTSeqCounts to create ExonCountSet object
metadata$condition <- as.factor(metadata$condition)


ecs <- read.HTSeqCounts(countfiles = paste0("MISO_counts/", metadata$SampleName, "_miso_sim12.txt"), design = metadata[,c("condition")])
sampleNames(ecs) = metadata$SampleName


counts <- counts(ecs)
gene.id <- fData(ecs)$geneID
ete.id <- fData(ecs)$exonID


gene.id.unique <- unique(gene.id)

#g = gene.id.unique[1]
g="FBgn0000054"

counts.tmp <- counts[gene.id == g, ] + 1
dim(counts.tmp)

design <- model.matrix(~condition, data=metadata)
design



Y <- t(counts.tmp)
Y

Y <- t(d)
Y



X <- as.matrix(design)



Loglik <- function( b, Y, X) {
  
  p <- ncol(X)
  q <- ncol(Y)
  b <- matrix(b, q, p)
  
  g <- exp(X %*% t(b))  # n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)
  
  # formula (7) is equivalent to likelihood from dirmult
  res <- 	sum(lgamma(gs) - lgamma(ys+gs) + rowSums(lgamma(Y+g) - lgamma(g)))
  # normalization coeff
  # norm <- sum(lgamma(ys + 1) - rowSums(lgamma(Y + 1)))
  # res <- res + norm
  
  -res
    
}


p <- ncol(X)
q <- ncol(Y)

b <- rep(0, q*p)

# optim
bh_opt <- optim(par=b, fn=Loglik, Y=Y, X=X, method="Nelder-Mead")

loglikh <- bh_opt$value
loglikh

bh <- matrix(bh_opt$par, q, p)
bh

gh <- exp(X %*% t(bh))
gh
th <- gh / rowSums(gh)
th



# nlm
bh_nlm <- nlm(Loglik, b, Y=Y, X=X)

loglikh <- bh_nlm$minimum
loglikh


bh <- matrix(bh_nlm$estimate, q, p)
bh

gh <- exp(X %*% t(bh))
gh
th <- gh / rowSums(gh)
th



#################################
# testing
#################################

coef <- 2


X0 <- X[,-coef, drop=FALSE]

#################################
# rocX
#################################


#ROC curve for score of prediction between SVM and NN
data(rocX.score)

head(rocX.score)

r <- rocX(rocX.score$score, rocX.score$label)

pdf("roc.pdf")
plot(r)
dev.off()


plot(r, cutoff.fpr = 0.8, col = c("red", "green"), cex.X = 3, pch.X = 2, lwd = 5, lty = 3)

# ROC cure of p_value and threshold.X is set at FDR = 0.05
# score = 1- p_value
# threshold.X = 1 -0.05
data(rocX.pvals)
score <- 1-do.call("cbind", rocX.pvals$pval)
score.X <- 1-do.call("cbind", rocX.pvals$padj)
label <- rep(FALSE, nrow(score))
label[rocX.pvals$indDE] <- TRUE
r1 <- rocX(score, label, score.X = score.X, threshold.X = 0.95, label.ordering = c(FALSE, TRUE))
plot(r1)
plot(r1, cutoff.fpr = 0.8, col = c("red", "green"), cex.X = 3, pch.X = 2, lwd = 5, lty = 3)

















