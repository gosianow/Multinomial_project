##############################################################################

# BioC 14
# Created 06 Oct 2014: update of dmFunctions version 4
# + change the estimations of common dispersion - pooled gamma+, but pi estimated per group

# Update 06 Oct 2014:
# + add sQTL analysis
# 08 Oct 2014
# + in dmAdj remove one group adjCROneGeneGroup and keep adjCROneGeneManyGroups
# + dge$counts must be a split list!


##############################################################################
# help functions
##############################################################################

repRow<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
repCol<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


##############################################################################
# log-likelihood for k parameters (pi)
##############################################################################


## Computes the log-likelihood
dmLogLik <- function(pi, gamma0, y){
  
  k <- ncol(y)
  N <- nrow(y)
  S <- rowSums(y)  
  l <- 0
  
  for(i in 1:N){  
    # i=1
    l <- l - sum(sapply(1:S[i], function(r){log(gamma0 + r - 1)}))    
    for(j in 1:k){   
      # j=3
      if(y[i,j] == 0) lij <- 0
      else lij <- sum(sapply(1:y[i,j], function(r){log(pi[j] * gamma0 + r - 1)}))    
      l <- l + lij      
    }
  }
  
  return(l)
  
}


## Computes the log-likelihood -- with gamma functions -- 
dmLogLikG <- function(pi, gamma0, y){
  
  N <- nrow(y)
  S <- rowSums(y)
  
  #l <- sum(sapply(1:N, function(i){ sum(lgamma(y[i,] + pi * gamma0) - lgamma(pi * gamma0)) }))
  l <- sum( rowSums( lgamma(y + repRow(pi, nrow(y)) * gamma0) - lgamma(repRow(pi, nrow(y)) * gamma0) ) )

  l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + l
  
  return(l)
  
}



##############################################################################
# log-likelihood for k-1 parameters (pi)
##############################################################################


## Computes the log-likelihood
dmLogLikkm1 <- function(pi, gamma0, y){
  ## pi has length of k-1
  
  k <- ncol(y)
  N <- nrow(y)
  S <- rowSums(y)  
  l <- 0
  
  pi <- c(pi, 1 - sum(pi))
  names(pi) <- colnames(y)
  
  for(i in 1:N){  
    # i=1
    l <- l - sum(log(gamma0 + 1:S[i] - 1))  # sum(sapply(1:S[i], function(r){log(gamma0 + r - 1)}))    
    for(j in 1:k){   
      # j=3
      if(y[i,j] == 0) lij <- 0
      else lij <- sum(log(pi[j] * gamma0 + 1:y[i,j] - 1)) # sum(sapply(1:y[i,j], function(r){log(pi[j] * gamma0 + r - 1)}))    
      l <- l + lij      
    }
  }
  
  return(l)
  
}


## Computes the log-likelihood -- with gamma functions -- for k-1 parameters
dmLogLikGkm1 <- function(pi, gamma0, y){
  ## pi has length of k-1
  
  N <- nrow(y)
  S <- rowSums(y)
  
  pi <- c(pi, 1 - sum(pi))
  names(pi) <- colnames(y)
  
  l <- sum( rowSums( lgamma(y + repRow(pi, nrow(y)) * gamma0) - lgamma(repRow(pi, nrow(y)) * gamma0) ) )
  
  l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + l
  
  return(l)
  
}



##############################################################################
# Score-function for k-1 parameters (pi) 
##############################################################################

## in version 4 there was mistake in calculation of dmScoreFunkm1
dmScoreFunkm1 <- function(pi, gamma0, y){  
  
  k <- ncol(y)
  N <- nrow(y) 
  ykm1 <- y[, -k, drop=FALSE]
  yk <- y[, k, drop=TRUE]
  pik <- 1-sum(pi)
  
  S <- rep(0, k-1)
  
  for(i in 1:N){ 
    # i=1
    if(yk[i] == 0) Sik <- 0
    else Sik <- sum(gamma0 / (pik * gamma0 + 1:yk[i] - 1))
			# else Sik <- sum(sapply(1:yk[i], function(r){gamma0 / (pik * gamma0 + r - 1) }))
    for(j in 1:(k-1)){
      # j=1
      if(y[i,j] == 0) Sij <- 0
      else Sij <- sum(gamma0 / (pi[j] * gamma0 + 1:ykm1[i,j] - 1)) 
				# else Sij <- sum(sapply(1:ykm1[i,j], function(r){gamma0 / (pi[j] * gamma0 + r - 1) })) 
        S[j] <- S[j] + Sij
    }
		 S <- S - Skj
  }
  return(S)
  
}



## Score-function for k-1 parameters
dmScoreFunGkm1 <- function(pi, gamma0, y){  
  ## pi has length of k-1
  
  k <- ncol(y)
  N <- nrow(y)
  ykm1 <- y[, -k, drop=FALSE]
  yk <- y[, k]
  pik <- 1-sum(pi)
  
  S <- gamma0 * colSums( digamma(ykm1 + repRow(pi, N) * gamma0) - digamma(repRow(pi, N) * gamma0) - (digamma(yk + gamma0 * pik) - digamma(gamma0 * pik)) ) 
  
  return(S)
  
}

# Error in digamma(ykm1 + repRow(pi, N) * gamma0) - digamma(repRow(pi, N) *  :
  # non-conformable arrays



##############################################################################
# observed FIM for k-1 parameters (pi) /// plenty of mistakes / WRONG FUNCTION
##############################################################################

# pi <- pi[-length(pi)]

#### Computing the inverce of observed Fisher Information Matrix
dmInvObsFIMkm1 <- function(pi, gamma0, y, inv = TRUE){  
  ## pi has length of k-1

  k <- ncol(y)
  N <- nrow(y)
  D <- rep(0, k-1)
  Djl <- 0
  pik <- 1-sum(pi)
  
  for(i in 1:N){
    # i=1
    if(y[i, k] == 0) Djl <- 0
    else Djl <- Djl +  sum(-gamma0^2 / (pik * gamma0 + 1:y[i,k] - 1) ^2)    
    for(j in 1:(k-1)){
      # j=1
      if(y[i,j] == 0) Dij <- 0
      else Dij <- sum(-gamma0^2 / (pi[j] * gamma0 + 1:y[i,j] - 1) ^2) 
      #print(Dij)
      D[j] <- D[j] + Dij
    }
  }
  
  H <- matrix(Djl, k-1, k-1)
  diag(H) <- diag(H) + D
  
  if(!inv)
    return(-H)
  
  invFIM <- solve(-H)
  return(invFIM)
}


##############################################################################
# expected FIM for k-1 parameters (pi) /// plenty of mistakes / WRONG FUNCTION
##############################################################################


#### Beta-Binomial
## does not work for big gamma0 !!!
dBetaBin <- function(x, n, a, b){
  beta(x+a, n-x+b) / beta(a,b) * choose(n,x)
}

pBetaBin <- function(x, n, a, b) {
  ## x must be >=1
  sapply(x, function(xx){ 
    1 - sum(dBetaBin(0:(xx-1), n, a, b))
    })
}

## works for big gamma0
dBetaBin2 <- function(x,n,a,b){ 
  sapply(x, function(xx){
    res <- lchoose(n,xx)
    if(xx==0) res <- res else res <- res + sum(log(a + 1:xx - 1))
    if(xx==n) res <- res else res <- res + sum(log(b + 1:(n-xx) -1))
    res <- res - sum(log(a + b + 1:n - 1))
    names(res) <- NULL
    exp(res)
  })
}

pBetaBin2 <- function(x, n, a, b) {
  ## x must be >=1
  sapply(x, function(xx){ 
    1 - sum(dBetaBin2(0:(xx-1), n, a, b))
  })
}


# x <- 1:5
# n <- 10
# a <- b <- 5
# 
# x = 1:y[i,k]; n = S[i]; a = pik * gamma0; b = (1 - pik) * gamma0

# dBetaBin(x, n, a, b)
# pBetaBin(x, n, a, b)
# 
# dBetaBin2(x,n,a,b)
# pBetaBin2(x,n,a,b)



  
#### Computing the inverce of expected Fisher Information Matrix
dmInvExpFIMkm1 <- function(pi, gamma0, y, inv = TRUE){  
  
  k <- ncol(y)
  N <- nrow(y)
  D <- rep(0, k-1)
  Djl <- 0
  pik <- 1-sum(pi)
  S <- rowSums(y)
   
  for(i in 1:N){
    # i=1
    if(y[i, k] == 0) Djl <- 0
    else Djl <- Djl +  sum(-gamma0^2 * pBetaBin2(1:y[i,k], n = S[i], a = pik * gamma0, b = (1 - pik) * gamma0) / (pik * gamma0 + 1:y[i,k] - 1) ^2)  
#     cat("Djl", Djl, fill = T)
    for(j in 1:(k-1)){
      # j=1
      if(y[i,j] == 0) Dij <- 0
      else Dij <- sum(-gamma0^2  * pBetaBin2(1:y[i,j], n = S[i], a = pi[j]*gamma0, b = (1 - pi[j]) * gamma0) / (pi[j] * gamma0 + 1:y[i,j] - 1) ^2) 
#       cat(Dij, fill = T)
      D[j] <- D[j] + Dij
    }
  }
  
  H <- matrix(Djl, k-1, k-1)
  diag(H) <- diag(H) + D
  
if(!inv)
  return(-H)

  invFIM <- solve(-H)
  return(invFIM)
  
}






##############################################################################
# estimates of pi -- Fisher scoring or constrOptim
##############################################################################
## estimate pi for given dispersion  with Fisher scoring -- one gene, one group

# y <- y[["g1"]]; gamma0 <- 1000;  mode = "obs"; epsilon=1e-05; maxIte = 1000; verbose = TRUE; plot = FALSE


dmOneGeneGroup <- function(y, gamma0, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose = FALSE, plot = FALSE){
  # NULL for filtered genes or genes with one exonNULL for filtered genes
  if(dim(y)[1] <= 1) return(NULL)
  
  ### y must be samples x exons
  y <- t(y)
  
  ### check for 0s in cols or rows
  keepCol <- colSums(y) > 0
  if(sum(keepCol) < 2) return(NULL)
  y <- y[, keepCol , drop=FALSE]
  df <- sum(keepCol)
  
  keepRow <- rowSums(y) > 0
  if(sum(keepRow) < 2) return(NULL) 
  y <- y[keepRow, , drop=FALSE]
  
  piInit <- colSums(y)/sum(y)
  k <- length(piInit)
  N <- nrow(y)
  
  if(plot){
    piX <- seq(0, 1, by = 0.01)
    piX <- piX[c(-1, -length(piX))]
    loglikY <- rep(0, length(piX))
    for(i in 1:length(loglikY))
      loglikY[i] <- dmLogLikkm1(piX[i], gamma0, y)
    
    plot(piX, loglikY, type="l", col="deeppink", lwd=4, main = gamma0)
    dev.off()
  }
  
  switch(mode, 
         
         constrOptim={
           epsilon <- .Machine$double.eps # *10
           
           ### must have constraint for SUM pi = 1
           ui <- rbind(diag(rep(1, k), k), diag(rep(-1, k), k), rep(1, k), rep(-1, k))
           ci <- c(rep(0, k), rep(-1, k), 1 - epsilon, -1 - epsilon)
           
           co <- constrOptim(piInit, f=dmLogLikG, g=NULL, ui=ui, ci=ci, control = list(fnscale = -1), gamma0=gamma0, y=y )
           
           piH <- co$par
           lik1 <- co$value
           
           piH
           sum(piH) == 1
         }, 
         
         
         constrOptim2={
           if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           
           epsilon <- .Machine$double.eps 
           
           ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
           ci <- c(rep(0, k-1), rep(-1, k-1), -1 + epsilon) # used to be with epsilon = 1e-5  but "initial value is not in the interior of the feasible region" because it is too big
           #         ui <- rbind(diag(rep(1, k-1)), diag(rep(-1, k-1)))
           #         ci <- c(rep(0, k-1), rep(-1, k-1))
           
           co <- constrOptim(piInit[-k], f=dmLogLikkm1, g=dmScoreFunkm1, ui=ui, ci=ci, control=list(fnscale = -1), gamma0=gamma0, y=y )
           #     co <- constrOptim(piInit[-k], f=dmLogLikkm1, g=NULL, ui=ui, ci=ci, control=list(fnscale = -1), gamma0=gamma0, y=y )
           
           piH <- co$par
           lik1 <- co$value
           
           if(verbose) cat("piH:",c(piH, 1-sum(piH)), fill = T)
           if(verbose) cat("lik1:", lik1, fill = T)
           
           piH <- c(piH, 1-sum(piH))
           names(piH) <- names(piInit)
         }, 
         
         
         optim2={
           if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           
           o <- optim(par = piInit[-k], fn = dmLogLikkm1, gr = dmScoreFunkm1, 
                      gamma0=gamma0, y=y,
                      method = "L-BFGS-B", lower = 0 + epsilon, upper = 1 - epsilon, control=list(fnscale = -1))
           
           piH <- o$par
           lik1 <- o$value
           piH <- c(piH, 1-sum(piH))
           names(piH) <- names(piInit)
           if(verbose) cat("piH:", piH, fill = T)
         }, 
         
         
         optim2NM={
           if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
           o <- optim(par = piInit[-k], fn = dmLogLikkm1, gr = NULL, 
                      gamma0=gamma0, y=y,
                      method = c("Nelder-Mead", "BFGS", "CG", "SANN")[1], control=list(fnscale = -1))
           
           piH <- o$par
           lik1 <- o$value
           piH <- c(piH, 1-sum(piH))
           names(piH) <- names(piInit)
           if(verbose) cat("piH:", piH, fill = T)
         }, 
         
         
         FisherScoring={
           # k-1 parameters for Fisher scoring
           piInitOrg <- piInit <- colSums(y)/sum(y)
           piMAX <- piH <- piInit[-k]
           lik1 <- 0
           likMAX <- lik2 <- dmLogLikkm1(piH, gamma0, y)
           ite <- 1
           conv <- TRUE 
           
           if(plot){
             piX <- seq(0, 1, by = 0.01)
             piX <- piX[c(-1, -length(piX))]
             loglikY <- rep(0, length(piX))
             for(i in 1:length(loglikY))
               loglikY[i] <- dmLogLikkm1(piX[i], gamma0, y)
             
             plot(piX, loglikY, type="l", col="deeppink", lwd=4, main = gamma0)
             abline(v = piH, lty =  3)
             plot(piX, loglikY, type="l", col="deeppink", lwd=6, main = gamma0, xlim=c(piH - 0.3, piH + 0.3))
             abline(v = piH, lty =  3)
             points(piH, lik2, pch="0")
           }
           
           # Iterations
           while(conv & ite <= maxIte){
             if(verbose) cat("\n gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
             if(verbose) cat("ite:", ite, fill = T)
             if(verbose) cat("lik2-lik1:", lik2 - lik1, fill = T)
             
             if(abs(lik2 - lik1) < epsilon) conv <- FALSE
             
             scoreFun <- dmScoreFunkm1(piH, gamma0, y)
             if(verbose) cat("scoreFun:", scoreFun, fill = T)
             
             if(mode=="exp") invFIM <- dmInvExpFIMkm1(piH, gamma0, y)
             else if(mode=="obs") invFIM <- dmInvObsFIMkm1(piH, gamma0, y)
             if(verbose) cat("invFIM:", invFIM, fill = T)
             
             # Updates parameter estimates
             piH <- piH + invFIM %*% scoreFun
             
             ## check if pi is negative, then restart with new init params
             if(any(c(piH, 1-sum(piH)) <= 0 )){
               cat("* Negative piH:",c(piH, 1-sum(piH)), fill = T)
               ### generate new starting params
               #         randInit <- runif(k)
               #         piInit <- randInit/sum(randInit)
               piInit <- rdirichlet(1, piInitOrg*10)
               while(any(piInit < 1e-10))
                 piInit <- rdirichlet(1, piInitOrg*10)
               cat("piInit:", piInit, fill = T)
               piH <- piInit[-k]
             }
             
             if(verbose) cat("piH:",c(piH, 1-sum(piH)), fill = T)
             
             lik1 <- lik2
             lik2 <- dmLogLikkm1(piH, gamma0, y)
             
             if(lik2 > likMAX){
               likMAX  <- lik2
               piMAX <- piH
             }
             
             
             if(plot)
               points(piH, lik2, pch="|")
             if(verbose) cat("lik2:", lik2, fill = T)
             ite <- ite + 1
           }
           
           
           
           if(ite > maxIte){
             piH <- piMAX
             lik1 <- likMAX
           }
           
           if(plot){
             points(piH, lik1, pch="*", cex=2, col="darkturquoise")
             dev.off()
           }
           
           piH <- c(piH, 1-sum(piH))
           names(piH) <- names(piInit) 
         }) ## End of FisherScoring
  

  keepCol[keepCol] <- piH
  piH <- keepCol
  
#   ## check if the sum of pi equals 1
#   if(sum(piH) != 1){
#     cat("gene:", colnames(y)[1], "gamma0:", gamma0, fill = T)
#     cat("**** Not 1. sum(piH) = ", sum(piH), fill = T)
#     cat("**** piH:",piH, fill = T)
#   }
  
  ## check if pi is negative
  if(any(piH < 0 | piH >1 )){
    cat("**** Negative piH:",piH, fill = T)
  }
  
  return(list(piH = as.matrix(piH), gamma0 = gamma0, logLik = lik1, df=df))
  
  
}


### -- one gene, many groups

# y <- y[["ENSG00000135778"]]; mode = "constrOptim2"; epsilon = 1e-05; maxIte = 1000; verbose=FALSE

dmOneGeneManyGroups <- function(y, group, gamma0, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose = FALSE){
  # NULL for filtered genes or genes with one exon
  k <- dim(y)[1]
  if(k <= 1) return(NULL)
  
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  piH = matrix(0, nrow = k, ncol = ngroups)
  logLik = rep(NA, ngroups)
  df = rep(0, ngroups)

  for(gr in 1:ngroups){
    # gr=2
    grIndx <- which(group == lgroups[gr])
    
    fitGr <- dmOneGeneGroup(y = y[, grIndx, drop = FALSE], gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose=verbose)

    if(is.null(fitGr)) return(NULL)
    
    piH[,gr] <- fitGr$piH
    logLik[gr] <- fitGr$logLik
    df[gr] <- fitGr$df
    
  }
    
  colnames(piH) <- names(df) <- lgroups
  logLik <- sum(logLik)
   
  return(list(piH = piH, gamma0 = gamma0, logLik = logLik, df=df))
  
  
}



##############################################################################
# adjustements to profile likelihood
##############################################################################

# g=1
# piH <- fit[[g]]$piH
# y <- y[[g]]


adjCROneGeneGroup <- function(y, gamma0, piH){
  # NULL for filtered genes or genes with one exon
  if(dim(y)[1] <= 1) return(NULL)
  
  ### y must be samples x exons
  y <- t(y)
  
  ### check for 0s in cols or rows
  keepCol <- colSums(y) > 0
  if(sum(keepCol) < 2) return(NULL)
  y <- y[, keepCol , drop=FALSE]
  
  keepRow <- rowSums(y) > 0
  if(sum(keepRow) < 2) return(NULL) 
  y <- y[keepRow, , drop=FALSE]
  
  N <- nrow(y) 
  piH <- piH[keepCol]
  
  ## 1/2 * log(det(N* FIM))
  adjCR <- log(det(N * dmInvObsFIMkm1(pi = piH[-length(piH)], gamma0, y, inv = FALSE)))/2

  if(adjCR == Inf)
    return(NULL)
  
  return(adjCR)
  
}



adjCROneGeneManyGroups <- function(y, group, gamma0, piH){  
  # NULL for filtered genes or genes with one exon
  k <- dim(y)[1]
  if(k <= 1) return(NULL)
  
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  adjCR = rep(0, ngroups)

  for(gr in 1:ngroups){
    # gr=1
    grIndx <- which(group == lgroups[gr])
    a <- adjCROneGeneGroup(y = y[, grIndx, drop = FALSE], gamma0, piH = piH[, lgroups[gr]])
    
    if(is.null(a)) return(NULL)
    adjCR[gr] <- a
    
  }
  
   adjCR <- sum(adjCR)
  
  return(adjCR)
  
}



# dge = dgeFit

dmAdj <- function(gamma0, dge, group=NULL, mcCores=20){
  
  if(is.null(group)) group <- dge$samples$group
  
  group <- as.factor(group)
  y <- dge$counts
  #   genes <- names(y)
  
  adj <- unlist(mclapply(seq(length(y)), function(g){  
    # g = "FBgn0013733"
    a <- adjCROneGeneManyGroups(y = y[[g]], group, gamma0 = gamma0, piH = dge$fit[[g]]$piH) 
    
    #       if(!is.null(a))
    #       names(a) <- genes[g]
    return(a)
    
  }, mc.cores=mcCores))
  
  return(adj)
  
}



dmSQTLAdj <- function(gamma0, dgeSQTL, mcCores=20){
  
  y <- dgeSQTL$counts
  
  adj <- unlist(mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){  
    # snp = 1
    
    g <- dgeSQTL$SNPs[snp, "gene_id"]
    gt <- dgeSQTL$genotypes[snp,]
    
    a <- adjCROneGeneManyGroups(y = y[[g]][, !is.na(gt)], group = as.factor(na.omit(gt)), gamma0 = gamma0, piH = dgeSQTL$fit[[snp]]$piH) 
    
    return(a)
    
  }, mc.cores=mcCores))
  
  return(adj)
  
}




# ## TODO 
# adjBNGene <- function(piH, gamma0, y){
#   
#   y <- t(y)  
#   ### check for 0s in cols or rows
#   keepRow <- rowSums(y) > 0
#   y <- y[keepRow, , drop=FALSE] 
#   N <- nrow(y)  
#   keepCol <- colSums(y) > 0
#   y <- y[, keepCol , drop=FALSE]  
#   piH <- piH[keepCol]
#   
#   ## 1/2 * log(det(N* FIM))
#   log(det(N * dmInvObsFIMkm1(pi = piH[-length(piH)], gamma0, y, inv = FALSE)))/2
#   
# }



##############################################################################
# calculate profile likelihood + adjustements
##############################################################################
# returns common likelihood = sum of all gene likelihoods

# gamma0=38196.6; dge <- dge; group=NULL; common = TRUE; adjust = TRUE; mode = "constrOptim2"; epsilon = 1e-05; maxIte = 1000; mcCores=40; verbose = FALSE


dmAdjustedProfileLik <- function(gamma0, dge, group=NULL, common = FALSE, adjust = FALSE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, mcCores=20, verbose = FALSE){
  
  cat("Gamma in optimize:", gamma0, fill = TRUE)
  
  dgeFit <- dmFit(dge, group=group, dispersion=gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose=verbose, mcCores = mcCores)

  loglik <- unlist(lapply(dgeFit$fit, function(g){g$logLik})) 
  
  
  if(common)
    loglik <- sum(loglik)
  
  cat("loglik:", loglik, fill = TRUE)
  
  if(!adjust)
    return(loglik)

  
  ## Cox-Reid adjustement
  adj <- dmAdj(gamma0, dge = dgeFit, group = group, mcCores=mcCores)
  
  if(common)
    adj <- sum(adj)
  
  adjloglik <- loglik - adj
  
  cat("adjloglik:", adjloglik, fill = TRUE)
  
  
  return(adjloglik)
  
}




dmSQTLAdjustedProfileLik <- function(gamma0, dgeSQTL, common = FALSE, adjust = FALSE, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, mcCores=20, verbose = FALSE){
  
  cat("Gamma in optimize:", gamma0, fill = TRUE)
  
  dgeSQTLFit <- dmSQTLFit(dgeSQTL, model = "full", dispersion=gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose=verbose, mcCores = mcCores)
  
  loglik <- unlist(lapply(dgeSQTLFit$fit, function(g){g$logLik})) 
  
  
  if(common)
    loglik <- sum(loglik)
  
  cat("loglik:", loglik, fill = TRUE)
  
  if(!adjust)
    return(loglik)
  
  
  ## Cox-Reid adjustement
  adj <- dmSQTLAdj(gamma0, dgeSQTL = dgeSQTLFit, mcCores=mcCores)
  
  if(common)
    adj <- sum(adj)
  
  adjloglik <- loglik - adj
  
  cat("adjloglik:", adjloglik, fill = TRUE)
  
  
  return(adjloglik)
  
}




##############################################################################
# calculate common dispersion 
##############################################################################

# group=NULL; adjust = FALSE; mode = "constrOptim2"; epsilon = 1e-05; maxIte = 1000; interval = c(0, 1e+5); tol = 1e-03; mcCores=1; verbose=FALSE


dmEstimateCommonDisp <- function(dge, group=NULL, adjust = FALSE, subset=Inf, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=20, verbose=FALSE){
  
  if(is.null(group)) group <- dge$samples$group

  if(subset != Inf){
    
    allGenes <- unique(as.character(dge$genes$gene_id))
    nGenes <- length(allGenes)
    
    ## take random genes
    genesSubset <- sample.int(nGenes, min(subset, nGenes))
    
    dge <- dge[dge$genes$gene_id %in% allGenes[genesSubset], ]
    dge$genes$gene_id <- as.factor(as.character(dge$genes$gene_id)) # otherwise number of levels astays as before subsetting      
    
  }

  out <- optimize(f = dmAdjustedProfileLik, interval = interval,
                  dge = dge, group = group, common = TRUE, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, mcCores = mcCores, verbose = verbose,
                  maximum = TRUE, tol = tol) 
  
  
  dge$commonDispersion <- out$maximum
  cat("** Connom dispersion: ", out$maximum, fill = TRUE)
  return(dge)
  
}



dmSQTLEstimateCommonDisp <- function(dgeSQTL, adjust = FALSE, subset=Inf, mode = "constrOptim2", epsilon = 1e-05, maxIte = 1000, interval = c(0, 1e+5), tol = 1e-03, mcCores=20, verbose=FALSE){
  

  if(subset != Inf){    
    nSNPs <- nrow(dgeSQTL$SNPs)
    
    ## take random SNP
    SNPSubset <- sample.int(nSNPs, min(subset, nSNPs))
  
    dgeSQTL$SNPs <- dgeSQTL$SNPs[SNPSubset, ]
    dgeSQTL$genotypes <- dgeSQTL$genotypes[SNPSubset, ]
    
  }
  
  out <- optimize(f = dmSQTLAdjustedProfileLik, interval = interval,
                  dgeSQTL = dgeSQTL, common = TRUE, adjust = adjust, mode = mode, epsilon = epsilon, maxIte = maxIte, mcCores = mcCores, verbose = verbose,
                  maximum = TRUE, tol = tol) 
  
  
  dgeSQTL$commonDispersion <- out$maximum
  cat("** Connom dispersion: ", out$maximum, fill = TRUE)
  return(dgeSQTL)
  
}



##############################################################################
# multiple group fitting 
##############################################################################


# dge <- dgeDM; group=NULL; dispersion=NULL; mode="constrOptim2"; epsilon = 1e-05; maxIte = 1000; verbose=FALSE; mcCores = 20


dmFit <- function(dge, group=NULL, dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20){

  if(is.null(group)) group <- dge$samples$group
  if(is.null(dispersion)) dispersion <- dge$commonDispersion
  gamma0 <- dge$dispersion <- dispersion

  group <- as.factor(group)
  ngroups <- nlevels(group)

y <- dge$counts
  genes <-  names(y)
  

  if( nlevels(group)==1 )
    model <- "null"
 else 
   model <- "full" 
     
 
 switch(model, 
        full={
          
          fit <- mclapply(seq(length(y)), function(g){  
            # g = "ENSG00000135778"
            # cat("Gene:", genes[g], fill = TRUE)
            f <- dmOneGeneManyGroups(y[[g]], group = group, gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  
            return(f)
          }, mc.cores=mcCores)
          
          names(fit) <- genes
          
        }, 
        null={
          
          fit <- mclapply(seq(length(y)), function(g){  
            dmOneGeneGroup(y[[g]], gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  
          }, mc.cores=mcCores)
          
          names(fit) <- genes
          
        })



# switch(model, ## use y instead of seq(length(y))
#        full={         
#          fit <- mclapply(y, function(yy){  
# 
#            f <- dmOneGeneManyGroups(yy, group = group, gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  
#            return(f)
#          }, mc.cores=mcCores)
#          
# #          names(fit) <- genes
#          
#        }, 
#        null={
#          
#          fit <- mclapply(y, function(yy){  
#            dmOneGeneGroup(yy, gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  
#          }, mc.cores=mcCores)
#          
# #          names(fit) <- genes
#          
#        })
    
  dge$fit <- fit
  
  return(dge)
   
}



# model="full"; dispersion=3000; mode="constrOptim2"; epsilon = 1e-05; maxIte = 1000; verbose=FALSE; mcCores = 20


dmSQTLFit <- function(dgeSQTL, model="full", dispersion=NULL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20){
  
  if(is.null(dispersion)) dispersion <- dgeSQTL$commonDispersion
  gamma0 <- dgeSQTL$dispersion <- dispersion

  y <- dgeSQTL$counts

  switch(model, 
         full={
                   
           fit <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){  
             # snp = 1
             # cat("SNP:", dgeSQTL$SNPs$SNP_id[snp], fill = TRUE)

             f <- dmOneGeneManyGroups(y = y[[dgeSQTL$SNPs[snp, "gene_id"]]][, !is.na(dgeSQTL$genotypes[snp,])], group = as.factor(na.omit(dgeSQTL$genotypes[snp,])), gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  
             return(f)
           }, mc.cores=mcCores)
           
           names(fit) <- dgeSQTL$SNPs$SNP_id
           
         }, 
         null={
           
           fit <- mclapply(seq(nrow(dgeSQTL$SNPs)), function(snp){  
             dmOneGeneGroup(y[[dgeSQTL$SNPs[snp, "gene_id"]]][, !is.na(dgeSQTL$genotypes[snp,])], gamma0 = gamma0, mode = mode, epsilon = epsilon, maxIte = maxIte, verbose = verbose)  
           }, mc.cores=mcCores)
           
           names(fit) <- dgeSQTL$SNPs$SNP_id
           
         })
  
  dgeSQTL$fit <- fit

  return(dgeSQTL)
  
}




#######################################################
#  group testing
#######################################################


# dge = dgeDM; mode="constrOptim2"; mcCores=30; verbose=FALSE

dmTest <- function(dge, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20){
  
  fit.full <- dge$fit
  nrgroups <- nlevels(as.factor(dge$samples$group))
  group <- as.factor(rep(1, length(dge$samples$group)))
  dge$genes$gene_id <- as.factor(dge$genes$gene_id)
  
  ## fit null model
  fit.null <- dmFit(dge=dge, group=group, dispersion=dge$dispersion, mode=mode, epsilon = epsilon, maxIte = maxIte, verbose= verbose, mcCores = mcCores)$fit
  
  LRT <- mclapply(unique(levels(dge$genes$gene_id)), function(g){
    # g = "FBgn0000014"
    if(verbose) cat("testing gene: ",g, fill = TRUE)
    
    if(is.null(fit.null[[g]]) || is.null(fit.full[[g]])) 
      return(data.frame(LR=NA, df=NA, PValue=NA, LLfull=NA, LLnull=NA))
    
      LLnull <- fit.null[[g]]$logLik

       LLfull <- sum(fit.full[[g]]$logLik)

      LR <-  2*(LLfull - LLnull)
      
      DFnull <- fit.null[[g]]$df - 1
      DFfull <- sum(fit.full[[g]]$df) - nrgroups
      
      df <- DFfull - DFnull
      
      pValue <- pchisq(LR, df = df , lower.tail = FALSE)
      
      return(data.frame(LR=LR, df=df, PValue=pValue, LLfull=LLfull, LLnull=LLnull))
    
  }, mc.cores=mcCores)
  
  
  LRT <- do.call(rbind, LRT)
  FDR <- p.adjust(LRT$PValue, method="BH")
  o <- order(LRT$PValue)
  
  table <- data.frame(GeneID=unique(levels(dge$genes$gene_id)), LR=LRT$LR, df=LRT$df, LLfull=LRT$LLfull, LLnull=LRT$LLnull , PValue=LRT$PValue, FDR=FDR)[o,]
  
  fit <- dge
  fit$fit.null <- fit.null
  fit$table <- table
  
  return(fit)
  
  
}




dmSQTLTest <- function(dgeSQTL, mode="constrOptim2", epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 20){
  
  fit.full <- dgeSQTL$fit

  ## fit null model
  fit.null <- dmSQTLFit(dgeSQTL=dgeSQTL, model = "null", dispersion = dgeSQTL$dispersion, mode=mode, epsilon = epsilon, maxIte = maxIte, verbose= verbose, mcCores = mcCores)$fit
  
  
  LRT <- mclapply(dgeSQTL$SNPs$SNP_id, function(snp){
    # snp = "snp_19_502623"
    if(verbose) cat("testing SNP: ", snp, fill = TRUE)
    
    if(is.null(fit.null[[snp]]) || is.null(fit.full[[snp]])) 
      return(data.frame(LR=NA, df=NA, PValue=NA, LLfull=NA, LLnull=NA))
    
    LLnull <- fit.null[[snp]]$logLik
    
    LLfull <- sum(fit.full[[snp]]$logLik)
    
    LR <-  2*(LLfull - LLnull)
    
    DFnull <- fit.null[[snp]]$df - 1
    DFfull <- sum(fit.full[[snp]]$df) - length(fit.full[[snp]]$df)
    
    df <- DFfull - DFnull
    
    pValue <- pchisq(LR, df = df , lower.tail = FALSE)
    
    return(data.frame(LR=LR, df=df, PValue=pValue, LLfull=LLfull, LLnull=LLnull))
    
  }, mc.cores=mcCores)
  
  
  LRT <- do.call(rbind, LRT)
  FDR <- p.adjust(LRT$PValue, method="BH")
  o <- order(LRT$PValue)
  
  table <- data.frame(SNP_id = dgeSQTL$SNPs$SNP_id, gene_id = dgeSQTL$SNPs$gene_id, LR=LRT$LR, df=LRT$df, LLfull=LRT$LLfull, LLnull=LRT$LLnull , PValue=LRT$PValue, FDR=FDR, stringsAsFactors = FALSE)[o,]
  
  fit <- dgeSQTL 
  fit$fit.null <- fit.null
  fit$table <- table
  
  return(fit)
  
  
}





#######################################################
# changing printing options from limma :D
#######################################################


# install.packages("/home/gosia/R/packages/limma", lib="/home/gosia/R/libraries/3.1.0/")


my.printHead <- function(x)
  #  Print leading 5 elements or rows of atomic object
  #  Gordon Smyth
  #  May 2003.  Last modified 14 April 2009.
{
  if(is.atomic(x)) {
    d <- dim(x)
    if(length(d)<2) which <- "OneD"
    if(length(d)==2) which <- "TwoD"
    if(length(d)>2) which <- "Array"
  } else {
    if(inherits(x,"data.frame")) {
      d <- dim(x)
      which <- "TwoD"
    } else {
      if(is.call(x))
        which <- "Call"
      else {
        if(is.recursive(x))
          which <- "Recursive"
        else
          which <- "Other"
      }
    }
  }
  switch(which,
         OneD={
           n <- length(x)
           if(n > 20) {
             print(x[1:5])
             cat(n-5,"more elements ...\n")
           } else
             print(x)
         },
         TwoD={
           n <- d[1]
           if(n > 10) {
             print(x[1:5,])
             cat(n-5,"more rows ...\n")
           } else
             print(x)
         },
         Array={
           n <- d[1]
           if(n > 10) {
             dn <- dimnames(x)
             dim(x) <- c(d[1],prod(d[-1]))
             x <- x[1:5,]
             dim(x) <- c(5,d[-1])
             if(!is.null(dn[[1]])) dn[[1]] <- dn[[1]][1:5]
             dimnames(x) <- dn
             print(x)
             cat(n-5,"more rows ...\n")
           } else
             print(x)
         },
         Recursive={
           n <- length(x)
           if(n) {
             i <- names(x)
             if(is.null(i)) i <- seq_len(n)
             if(length(i) >= 2) i <- i[1:2]
             for (what in i) {
               y <- x[[what]]
               cat("$",what, "\n",sep="")
              print(y) #Recall(y)
         
             }
             cat(n-2,"more elements ...\n")
           }
         },
         Call=,
         Other=print(x)
  )
}



unlockBinding("printHead", as.environment("package:limma"))
assignInNamespace("printHead", my.printHead, ns="limma", envir=as.environment("package:limma"))
assign("printHead", my.printHead, as.environment("package:limma"))
lockBinding("printHead", as.environment("package:limma"))










































