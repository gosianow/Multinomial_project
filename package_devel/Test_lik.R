
### matrix with 0s
gamma0 = 10
y = matrix(c(0, 0, 35, 70), nrow = 2)
pi = 0.3

DM:::dm_lik(pi, gamma0, y)
DM:::dm_likG(pi, gamma0, y)

DM:::dm_score(pi, gamma0, y)
DM:::dm_scoreG(pi, gamma0, y)




gamma0 = 10
y = matrix(c(35, 70), nrow = 2)
pi = 0.3

DM:::dm_lik(pi, gamma0, y)
DM:::dm_likG(pi, gamma0, y)

DM:::dm_score(pi, gamma0, y)
DM:::dm_scoreG(pi, gamma0, y)




gamma0 = 10
y = matrix(c(35, 70, 100, 100), nrow = 2)
pi = 0.3

DM:::dm_lik(pi, gamma0, y)
DM:::dm_likG(pi, gamma0, y)




gamma0 = 10
y = matrix(c(0, 0, 0, 70), nrow = 2)
pi = 0.3

DM:::dm_lik(pi, gamma0, y)
DM:::dm_likG(pi, gamma0, y)

dm_score(pi, gamma0, y)
dm_scoreG(pi, gamma0, y)



pi_init <- rowSums(y)/sum(y)
pi <- pi_init[-length(pi_init)]

dm_lik(pi, gamma0, y)
dm_likG(pi, gamma0, y)

dm_score(pi, gamma0, y)
dm_scoreG(pi, gamma0, y)




gamma0 = 10
y = matrix(c(1, 100, 1, 200), nrow = 2)
pi = 0.0001

DM:::dm_lik(pi, gamma0, y)
DM:::dm_likG(pi, gamma0, y)

dm_score(pi, gamma0, y)
dm_scoreG(pi, gamma0, y)





################################################################################


dfe <- width(dt@counts) - 1
dfo <- dt@table[names(dfe), "df"]

table(dfe == dfo)

zero_groups <- which(!dfe == dfo)

zero_counts <- dt@counts[names(zero_groups)]

zero_gamma0 <- dt@tagwise_dispersion[names(zero_groups)]

likA <- numeric(length(zero_gamma0))
likB <- numeric(length(zero_gamma0))

for(i in 1:length(zero_gamma0)){
  # i = 20
  zero_counts[[i]]
  
  y <- zero_counts[[i]][, 4:6]
  y
  gamma0 <- zero_gamma0[i]


  ### keep 0 features
  y_init <- y
  y_init[y_init == 0] <- 1

  pi_init <- rowSums(y_init)/sum(y_init)

  k <- length(pi_init) ## k - number of exons


  ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
  ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) 


  co <- constrOptim(pi_init[-k], f = DM:::dm_likG, grad = DM:::dm_scoreG, ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol), gamma0 = gamma0, y = y)

  pi <- co$par
  pi <- c(pi, 1-sum(pi))
  lik <- co$value

  pi
  lik

  likA[i] <- lik

  ### remove 0 features

  keep_row <- rowSums(y) > 0
  y <- y[keep_row, , drop=FALSE]

  pi_init <- rowSums(y)/sum(y)

  k <- length(pi_init) ## k - number of exons


  ui <- rbind(diag(rep(1, k-1), k-1), diag(rep(-1, k-1), k-1), rep(-1, k-1))
  ci <- c(rep(0, k-1), rep(-1, k-1), -1 + .Machine$double.eps) 


  co <- constrOptim(pi_init[-k], f = DM:::dm_likG, grad = DM:::dm_scoreG, ui = ui, ci = ci, control = list(fnscale = -1, reltol = prop_tol), gamma0 = gamma0, y = y)

  pi <- co$par
  pi <- c(pi, 1-sum(pi))
  lik <- co$value

  pi
  lik

  likB[i] <- lik

  
}


likA - likB



################################################################################


y <- dt@counts[["ENSG00000010704"]][, 1:3]
gamma0 <- dt@tagwise_dispersion["ENSG00000010704"]


dm_fitOneGeneOneGroup(y, gamma0, prop_mode = c("constrOptim", "constrOptimG")[2], prop_tol = 1e-12, verbose = FALSE)


y <- dt@counts[["ENSG00000010704"]][, 1:3] 
y[y == 0] <- 1
gamma0 <- dt@tagwise_dispersion["ENSG00000010704"]


dm_fitOneGeneOneGroup(y, gamma0, prop_mode = c("constrOptim", "constrOptimG")[2], prop_tol = 1e-12, verbose = FALSE)



y <- dt@counts[["ENSG00000010704"]]
gamma0 <- dt@tagwise_dispersion["ENSG00000010704"]


dm_fitOneGeneOneGroup(y, gamma0, prop_mode = c("constrOptim", "constrOptimG")[2], prop_tol = 1e-12, verbose = FALSE)


y <- dt@counts[["ENSG00000010704"]] + 1
gamma0 <- dt@tagwise_dispersion["ENSG00000010704"]


dm_fitOneGeneOneGroup(y, gamma0, prop_mode = c("constrOptim", "constrOptimG")[2], prop_tol = 1e-12, verbose = FALSE)















