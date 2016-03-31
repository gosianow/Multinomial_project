
# install.packages("dirmult", dependencies=TRUE, lib="/home/gosia/R/libraries/3.0.2")

# install.packages("dirmult", dependencies=TRUE)

source("/home/gosia/R/multinomial_project/DirMulti/DirMulti.R")


#######################################################
# run example from DirMult
#######################################################


sim.obj <- DirmultSim()

Y <- sim.obj$Y
X <- sim.obj$X


dm.obj <- DirmultGrpGrid(Y, X, model="dirmult")

dir.obj <- DirmultGrpGrid(Y, X, model="dir")

mult.obj <- DirmultGrpGrid(Y, X, model="mult")


names(dm.obj[[1]])
names(dm.obj[[2]])



#######################################################
# run example from DirMult on data simulated from DM 
#######################################################


ysim <- t(simulate_from_dm_table(m = 1, n = 10, pi = c(1/3, 1/3, 1/3), g0 = 100, nm = 100, tot = "uni", nd = 3, BPPARAM = BiocParallel::MulticoreParam(workers = 1)))


Loglik <- function(Y, X, b, model = "dirmult") {
  # Compute the log likelihood, constant part discarded
  # The likelihood is scaled. Be careful when computing AIC BIC
  
  g <- exp(X %*% t(b))  # n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)
  if (model == "dirmult") {
    res <-  sum(lgamma(gs) - lgamma(ys+gs) + rowSums(lgamma(Y+g) - lgamma(g))) # formula (7)
  }
  if (model == "mult") {
    res <-  sum(rowSums(Y * log(g)) - ys * log(gs))
  }
  if (model == "dir") {
    res <-  sum(lgamma(gs) + rowSums((g-1) * log(Y) - lgamma(g)))
  }

  res / nrow(X)
}


S <- function(Y, X, b, model = "dirmult") {
  # Compute the Score function at b
  
  S <- 0
  g <- exp(X %*% t(b))  # n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)
  
  if (model == "dirmult") {
    S <-  t((digamma(gs) - digamma(ys+gs) + digamma(Y+g) - digamma(g)) * g) %*% X
  }
  
  if (model == "mult") {
    S <- t((Y / g - ys / gs) * g) %*% X
  }
  
  if (model == "dir") {
    S <-  t((digamma(gs)  - digamma(g) + log(Y)) * g) %*% X
  }
  
  S / nrow(X)
}


H <- function(Y, X, b, model = "dirmult"){
  # Compute the diagonal of the hessian matrix at b

  H <- 0
  g <- exp(X %*% t(b))  # n * q
  gs <- rowSums(g)
  ys <- rowSums(Y)
  if (model == "dirmult") {
      H <- t((trigamma(gs) - trigamma(ys+gs) + trigamma(Y+g) - trigamma(g)) * g^2 + (digamma(gs) - digamma(ys+gs) + digamma(Y+g) - digamma(g)) * g) %*% X^2
  }
  if (model == "mult") {
      H <- t((-Y / g^2 + ys / gs^2) * g^2  +
          (Y / g - ys / gs) * g) %*% X^2
  }
  if (model == "dir") {
      H <- t((trigamma(gs) - trigamma(g)) * g^2  +
          (digamma(gs)  - digamma(g) + log(Y)) * g) %*% X^2
  }
  # Divided by the sample size
  H / nrow(X)
}

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
  res <-  sum(lgamma(gs) - lgamma(ys+gs) + rowSums(lgamma(Y+g) - lgamma(g)))
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



















