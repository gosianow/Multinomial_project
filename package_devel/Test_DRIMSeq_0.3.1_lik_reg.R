library(gtools)
library(BiocParallel)

# m = 10; n = 5; pi = c(1/3, 1/3, 1/3); g0 = 100; nm = 100; tot = "uni"; nd = 3; BPPARAM = BiocParallel::MulticoreParam(workers = 4)


simulate_from_dm_table <- function(m = 10, n = 5, pi = c(1/3, 1/3, 1/3), g0 = 100, nm = 100, tot = "nbinom", nd = 3, BPPARAM = BiocParallel::MulticoreParam(workers = 4)){
  # m - number of genes
  # n - total number of samples/replicates
  # pi - proportions of features 
  # g0 - dispersion gamma0 (can be a vector of length m)
  # nm - total number of counts per gene (can be a vector of length m)
  # tot - method to generate total number of counts per gene: nbinom, norm, uni
  # nd - dispersion of nm if simulated from nbinom or norm
  
  if(length(g0 == 1))
    g0 <- rep(g0, m)
  if(length(nm == 1))
    nm <- rep(nm, m)
  
  sim <- bplapply(1:m, function(i){
    # i=1
    
    g_dir_org <- pi * g0[i] # gamma
    
    # simulate dirichlets
    g_dir <- rdirichlet(n, g_dir_org)
    
    # simulate total counts
    if(tot == "nbinom")
      t <- rnbinom(n, mu = nm[i], size = nd)
    if(tot == "norm")
      t <- round(rnorm(n, mean = nm[i], sd = nd))  
    if(tot == "uni")
      t <- rep(nm[i], n)
    
    # simulate multinomial
    dm <- sapply(1:n, function(j) rmultinom(1, prob = g_dir[j, ], size = t[j]))    
   
    rownames(dm) <- paste0("g", rep(i, length(g_dir_org)), ":f", 1:length(g_dir_org))
    
    return(dm)
    
  }, BPPARAM = BPPARAM)
  
sim <- do.call(rbind, sim)

sim

}



dm_logit <- function(x){
  n <- length(x)
  log(x[-n]/x[n]) 
}

n <- 10
pi_org <- c(0.2, 0.6, 0.2)
b_org <- dm_logit(pi_org)
b_org

x <- matrix(1, nrow = n) # 10 samples, 1 covariate
b <- matrix(b_org, nrow = 1, ncol = length(b_org)) # 2 transcripts
b


z <- exp(x %*% b)
pi <- z/(1 + rowSums(z))
pi <- cbind(pi, 1 - rowSums(pi))
pi


ysim <- simulate_from_dm_table(m = 1, n = n, pi = pi[1, ], g0 = 100, nm = 1000, tot = "uni", nd = 3, BPPARAM = BiocParallel::MulticoreParam(workers = 1))



dm_lik_regG(b = b, x = x, gamma0 = 100, y = ysim)

dm_likG(pi = pi_org[-length(pi_org)], gamma0 = 100, y = ysim)


dm_lik_regG_neg <- function(b, x, gamma0, y){
  - dm_lik_regG(b, x, gamma0, y)
}

b_init <- c(1, 1)

# nlm 
bh_nlm <- nlm(dm_lik_regG_neg, p = b_init, x = x, gamma0 = 100, y = ysim)
bh_nlm$estimate



# optim
bh_opt <- optim(par = b_init, fn = dm_lik_regG, x = x, gamma0 = 100, y = ysim, method="Nelder-Mead", control = list(fnscale = -1))

bh_opt$par


# optim
bh_opt <- optim(par = b_init, fn = dm_lik_regG, x = x, gamma0 = 100, y = ysim, method = "BFGS", control = list(fnscale = -1))

bh_opt$par


b_est <- matrix(bh_opt$par, nrow = 1, ncol = length(b_org))

z <- exp(x %*% b_est)
pi <- z/(1 + rowSums(z))
pi <- cbind(pi, 1 - rowSums(pi))
pi[1, ]










