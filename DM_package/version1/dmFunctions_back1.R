##############################################################################

# BioC 14
# Created 19 Aug 2014

# Update 23 Aug 2014:

# all loglik, score, and FIM functions for k and k-1 parameters

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
# Score-function for k parameters (pi)
# -- ERROR: does not cut zero 
##############################################################################


## Score-function
dmScoreFun <- function(pi, gamma0, y){  
  ## pi has length of k-1
  k <- ncol(y)
  N <- nrow(y)
  S <- rep(0, k)
  
  for(i in 1:N){
    for(j in 1:k){
      if(y[i,j] == 0) Sij <- 0
      else Sij <- sum(sapply(1:y[i,j], function(r){gamma0 / (pi[j] * gamma0 + r - 1)}))   
      S[j] <- S[j] + Sij
    }
  }
  return(S)
  
}


## Score-function
dmScoreFunG <- function(pi, gamma0, y){  
  
  S <- gamma0 * (colSums(digamma(y + repRow(pi, nrow(y)) * gamma0) - digamma(repRow(pi, nrow(y)) * gamma0) ) )
  
  return(S)
  
}



##############################################################################
# Score-function for k-1 parameters (pi) 
##############################################################################


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
    else Sik <- sum(gamma0 / (pik * gamma0 + 1:yk[i] - 1)) # sum(sapply(1:yk[i], function(r){gamma0 / (pik * gamma0 + r - 1) }))   
    for(j in 1:(k-1)){
      # j=1
      if(y[i,j] == 0) Sij <- 0
      else Sij <- sum(gamma0 / (pi[j] * gamma0 + 1:ykm1[i,j] - 1)) - Sik # sum(sapply(1:ykm1[i,j], function(r){gamma0 / (pi[j] * gamma0 + r - 1) })) - Sik
        S[j] <- S[j] + Sij
    }
  }
  return(S)
  
}



## Score-function for k-1 parameters
dmScoreFunGkm1 <- function(pi, gamma0, y){  
  ## pi has length of k-1
  
  k <- ncol(y)
  N <- nrow(y)
  ykm1 <- y[, -k, drop=FALSE]
  yk <- y[, k, drop=FALSE]
  pik <- 1-sum(pi)
  
  S <- gamma0 * colSums( digamma(ykm1 + repRow(pi, N) * gamma0) - digamma(repRow(pi, N) * gamma0) - (digamma(yk + gamma0 * pik) - digamma(gamma0 * pik)) ) 
  
  return(S)
  
}



##############################################################################
# FIM for k parameters (pi)
# -- ERROR
##############################################################################


## Computing the inverce of observed Fisher Information Matrix
dmInvObsFMI <- function(pi, gamma0, y){  

  k <- ncol(y)
  N <- nrow(y)
  D <- rep(0, k)
  
  for(i in 1:N){
    for(j in 1:k){
      if(y[i,j] == 0) Dij <- 0
      else Dij <- sum(sapply(1:y[i,j], function(r){gamma0^2 / (pi[j] * gamma0 + r - 1) ^2}))   
      D[j] <- D[j] + Dij
    }
  }
  
  if(!any(D==0)){
    invD <- 1/D
    return(diag(invD)) ## should be -diag(invD)
  }else{
    return(NULL)
  }
  
}


## Beta-Binomial
dBetaBin <- function(x, n, a, b) {
  beta(x+a, n-x+b)/beta(a,b)*choose(n,x)
}

pBetaBin <- function(x, n, a, b) {
    ## x must be >=1
 1 - sum(dBetaBin(0:(x-1), n, a, b))
}

dBetaBin2 <- function(x,n,a,b){
  ## does not work for x - vector
  res <- lchoose(n,x)
  if(x==0) res <- res else res <- res + sum(sapply(1:x, function(xx){log(a + xx - 1)}))
  if(x==n) res <- res else res <- res + sum(sapply(1:(n-x), function(xx){log(b + xx -1)}))
  res <- res - sum(sapply(1:n, function(xx){log(a + b + xx - 1)}))
  exp(res)
}



## Computing the inverce of expected Fisher Information Matrix
dmInvExpFMI <- function(pi, gamma0, y){  
  
  S <- rowSums(y)
  k <- ncol(y)
  N <- nrow(y)
  D <- rep(0, k)
  
  for(i in 1:N){
    for(j in 1:k){
      if(y[i,j] == 0) Dij <- 0
      else Dij <- sum(sapply(1:y[i,j], function(r){ pBetaBin(r, n = S[i], a = pi[j]*gamma0, b = (1 - pi[j]) * gamma0) * gamma0^2 / (pi[j] * gamma0 + r - 1) ^2}))   
      D[j] <- D[j] + Dij
    }
  }
  
  if(!any(D==0)){
    invD <- 1/D
    return(diag(invD)) ## should be -diag(invD)
  }else{
    return(NULL)
  }
  
}



##############################################################################
# FIM for k-1 parameters (pi)
##############################################################################

# pi <- pi[1]

## Computing the inverce of observed Fisher Information Matrix
dmInvObsFMIkm1 <- function(pi, gamma0, y){  
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
      D[j] <- D[j] + Dij
    }
  }
  
  FIM <- matrix(Djl, k-1, k-1)
  diag(FIM) <- diag(FIM) + D
  
  invFIM <- solve(-FIM)
  return(invFIM)
}








##############################################################################
# estimates of pi -- Fisher scoring or constrOptim
##############################################################################
## estimate pi for given dispersion  with Fisher scoring -- per gene

# mode = "constrOptim"; epsilon=1e-05; mcCores=20
# y <- y[["g36"]]
# gamma0 <- gamma0[1]

dmFitGene <- function(y, gamma0, mode = "constrOptim", epsilon=1e-05){
  # print(rownames(y)[1])
  ### y must be samples x exons
  y <- t(y)

  ### check for 0s in cols or rows
  keepRow <- rowSums(y) > 0
  y <- y[keepRow, , drop=FALSE]
  
  keepCol <- colSums(y) > 0
  y <- y[, keepCol , drop=FALSE]
  
  
  piInit <- colSums(y)/sum(y)
  k <- length(piInit)
  
  if(mode == "constrOptim"){
    
    epsilon <- .Machine$double.eps * 10
    
    ui <- rbind(diag(rep(1, k)), diag(rep(-1, k)), rep(1, k), rep(-1, k))
    ci <- c(rep(0, k), rep(-1, k), 1 - epsilon, -1 - epsilon)
    
#     ui <- rbind(diag(rep(1, k)), diag(rep(-1, k)))
#     ci <- c(rep(0, k), rep(-1, k))
    
    co <- constrOptim(piInit, f=dmLogLikG, g=NULL, ui=ui, ci=ci, control=list(fnscale = -1), gamma0=gamma0, y=y )

    piH <- co$par
    lik1 <- co$value

  }else{
    # k-1 parameters for Fisher scoring
    piH <- piInit[-k]
    lik1 <- 0
    lik2 <- epsilon*10
    ite <- 1
    conv <- TRUE 
    
    # Iterations
    while(conv){
      print(ite)
      if(abs(lik2 - lik1) < epsilon) conv <- FALSE
      
      scoreFun <- dmScoreFunkm1(piH, gamma0, y)
      #print(scoreFun)
      
      if(mode=="exp") invFMI <- dmInvExpFMI(piH, gamma0, y)
      else if(mode=="obs") invFMI <- dmInvObsFMIkm1(piH, gamma0, y)
      #print(invFMI)
      
      lik1 <- dmLogLikkm1(piH, gamma0, y)
      
      # Updates parameter estimates
      piH <- piH + invFMI %*% scoreFun
      print(piH)
      
      lik2 <- dmLogLikkm1(piH, gamma0, y)
      ite <- ite + 1
    }
    
    piH <- c(piH, 1-sum(piH))
    names(piH) <- names(piInit)
  }
  
  keepCol[keepCol] <- piH
  piH <- keepCol

  return(list(piH = piH, gamma0 = gamma0, logLik = lik1))
  
  
}




### like glmFit in edgeR -- estimate pi for given dispersion for all genes

# gamma0 <- sum(g.dir.org)
# y <- split(as.data.frame(dge$counts), dge$genes$gene.id)


dmFit <- function(y, gamma0, mode = "constrOptim", mcCores=20){

  if(length(gamma0==1))
    gamma0 <- rep(gamma0, length(y))
  
  fit <- mclapply(seq(length(y)), function(g){ dmFitGene(y[[g]], gamma0[g], mode = mode)  }, mc.cores=mcCores)
  
  names(fit) <- names(y)
  
  return(fit)
  
}



##############################################################################
# calculate profile likelihood + adjustements
##############################################################################

# gamma0 <- sum(g.dir.org)
# y <- split(as.data.frame(dge$counts), dge$genes$gene.id)


# returns common likelihood = sum of all genes likelihoods
dmAdjustedProfileLik <- function(gamma0, y, adjust = FALSE, mode = "constrOptim", mcCores=20, common = FALSE){
  
  
  fit <- dmFit(y, gamma0, mode = mode, mcCores=mcCores)
  

  loglik <- unlist(lapply(fit, function(g){g$logLik})) 
  
  if(common)
    loglik <- sum(loglik)
  
  if(!adjust)
    return(loglik)

  ### add adjustements
  
  
  
  
  
  
  
  
  
  
  
}



##############################################################################
# calculate common dispersion 
##############################################################################


# y <- split(as.data.frame(dge$counts), dge$genes$gene.id)
# adjust = FALSE; mode = "constrOptim"; mcCores=20; interval = c(0, 1e+2); tol = 1e-05


dmEstimateCommonDisp <- function(y, adjust = FALSE, mode = "constrOptim", mcCores=20, interval = c(0, 1e+5), tol = 1e-05){
  
  
  ### function to optimize -- common loglikelihood
  fun <- function(gamma0, y, adjust, mode, mcCores){
    dmAdjustedProfileLik(gamma0, y, adjust, mode, mcCores, common = TRUE) 
  }
  
  out <- optimize(f = fun, interval = interval,
                  y = y, adjust = adjust, mode = mode, mcCores = mcCores,
                  maximum = TRUE, tol = tol) 
  
  return(out$maximum)
  
}




##############################################################################
# 2 group fitting 
##############################################################################























































