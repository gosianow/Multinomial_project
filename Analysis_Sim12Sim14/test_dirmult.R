library("dirmult")
library("gtools")

setwd("/home/gosia/Multinomial_project")


################################################################################
# profile likelihood plot
################################################################################


#dir.par1 <- c(85, 5, 5, 5)
dir.par1 <- c(10, 10, 10, 10)
n <- 3

# simulate some dirichlets
g.dir <- rbind( rdirichlet(n, dir.par1 ))
g.dir

# give some variation in total counts
tot <- rnbinom(n,mu=100,size=3)
tot

# simulate transcript counts
d <- sapply(1:n,function(u) rmultinom(1,prob=g.dir[u,],size=tot[u]))
d <- t(d) + 1

colnames(d) <- paste0("I", 1:4)
rownames(d) <- paste0("S", 1:n)

d


source("/home/gosia/R/R_Multinomial_project/dirmult.R")

data <- d[1:n, ]
data 

g1 <- dirmult( data, trace=FALSE, initscalar=100)
g1
dirmult.summary( data, g1)


gamma <- colSums(data)/sum(data)


gamma.v <- c(seq(2, 100, 2), seq(120, 2000, 20), seq(2200, 100000, 200), seq(100000, 1000000, 10000))
#gamma.v <- c(seq(2, 100, 2), seq(120, 1000, 20), seq(1200, 20000, 200))

theta.v <- 1/(1+gamma.v)
PLL.v <- rep(NA, length(theta.v))

for(i in 1:length(theta.v)){
  # i=151
  #cat(paste0(i, "\n"))
  
  PLL <-  estProfLogLik(data, theta=theta.v[i], trace=FALSE)
  PLL.v[i] <- PLL$loglik
  
  
}

pdf(paste0("Profile_lokelihood/PLL_", gsub(":", "-" ,gsub(" ", "_", as.character(date()))) , ".pdf") )

plot(theta.v[1:150], PLL.v[1:150], type="l", main=paste0("n = ", n, "\n dir par = ", paste(dir.par1, collapse=", ")))
abline(v=g1$theta, col=2 )
abline(v=1/(1+sum(dir.par1)), lty=2)


plot(theta.v[1:150], PLL.v[1:150], type="l", main=paste0("n = ", n, "\n dir par = ", paste(dir.par1, collapse=", ")), xlim=c(0, 0.1))
abline(v=g1$theta, col=2 )
abline(v=1/(1+sum(dir.par1)), lty=2)


plot(theta.v, PLL.v, type="l", log="x", main=paste0("n = ", n, "\n dir par = ", paste(dir.par1, collapse=", ")))
abline(v=g1$theta, col=2 )
abline(v=1/(1+sum(dir.par1)), lty=2)


dev.off()



################################################################################
# profile likelihood check
################################################################################


dir.par1 <- c(85, 5, 5, 5)
# dir.par1 <- c(10, 10, 10, 10)
n <- 3

# simulate some dirichlets
g.dir <- rbind( rdirichlet(n, dir.par1 ))
g.dir

# give some variation in total counts
tot <- rnbinom(n,mu=100,size=3)
tot

# simulate transcript counts
d <- sapply(1:n,function(u) rmultinom(1,prob=g.dir[u,],size=tot[u]))
d <- t(d)

colnames(d) <- paste0("I", 1:4)
rownames(d) <- paste0("S", 1:n)


d[ d == 0 ] <- 1

### run dirmult
source("/home/gosia/R/R_Multinomial_project/dirmult.R")

data <- d
 
data

rowSums(data)

g1 <- dirmult( data, trace=FALSE)
g1
dirmult.summary( data, g1)

## MOM
gamma <- colSums(data)/sum(data)
gamma



gamma.v <- 100

theta.v <- 1/(1+gamma.v)
theta.v

PL <- estProfLogLik(data, theta=theta.v, trace=FALSE)



gp.t <- gridProf(data, theta=theta.v, from=-0.01, to=0.01, len=20)

gp.e <- gridProf(data, theta=g1$theta, from = - 1e-6, to = 1e-5, len=10)
gp.e

PL <- estProfLogLik(data, theta=g1$theta, initPi=g1$pi, epsilon=1e-10, trace=FALSE, maxit=1e5)

theta=g1$theta; trace=FALSE



plot(gp.t$theta, gp.t$loglik, type="l")
abline(v=theta.v, col=2)
dev.off()



plot(gp.e$theta, gp.e$loglik, type="l")
abline(v=g1$theta, col=2)
dev.off()




adapGridProf(data, delta = qchisq(0.95,df=1)*2, stepsize=50)







################################################################################
# testing between 2 groups
################################################################################


dir.par1 <- c(85, 5, 5, 5)
dir.par2 <- c(5, 5, 5, 85)

# simulate some dirichlets
g.dir <- rbind( rdirichlet(3, dir.par1 ),  rdirichlet(3, dir.par2 ) )

# g.dir <- rbind( rdirichlet(3, c(17, 1, 1, 1) ), rdirichlet(3, c(1, 1, 1, 17) ) )

g.dir

# give some variation in total counts
tot <- rnbinom(6,mu=100,size=3)
tot

# simulate transcript counts
d <- sapply(1:6,function(u) rmultinom(1,prob=g.dir[u,],size=tot[u]))
d <- t(d) + 1

colnames(d) <- paste0("I", 1:4)
rownames(d) <- paste0("S", 1:6)

d


source("/home/gosia/R/R_Multinomial_project/dirmult.R")

data <- d[1:3, ]
data <- d

g1 <- dirmult( data, trace=FALSE, initscalar=100)
g1
dirmult.summary( data, g1)


g2 <- dirmult( d[4:6, ], trace=FALSE)
g2


null <- dirmult( d, trace=FALSE)
null
dirmult.summary( d, null )




lr <- -2*( null$loglik - (g1$loglik+g2$loglik) )

pchisq(lr,df=3,lower.tail=FALSE)


# 

rmultinom(5, size=100, prob=c(0.1, 0.4, 0.5))


#######



data(us)
fit <- dirmult(us[[1]], epsilon=10^(-4), trace=FALSE)
fit
dirmult.summary(us[[1]], fit)





###################################
# loglik
###################################




ll <- function (x, t)
  log(t + x - 1)



x <- t(d)
t <- full$gamma

loglik <- function (x, t)
{
  
  
  l <- 0
  ts <- sum(t)
  nc <- ncol(x)
  for (j in 1:nrow(x)) {
    
    l <- l - sum(unlist(lapply(list(1:(rowSums(x)[j])), ll, t = ts)))
    for (i in 1:nc) {
      if (x[j, i] == 0)
        lij <- 0
      else lij <- sum(unlist(lapply(list(1:(x[j, i])),ll, t = t[i])))
      l <- l + lij
    }
    
    
  }
  l
}


loglik(x, t)




###################################
# dmultinom The Multinomial Distribution
###################################



dmultinom <- function (x, size = NULL, prob, log = FALSE)
{
  K <- length(prob)
  if (length(x) != K)
    stop("x[] and prob[] must be equal length vectors.")
  if (any(prob < 0) || (s <- sum(prob)) == 0)
    stop("probabilities cannot be negative nor all 0")
  prob <- prob/s
  x <- as.integer(x + 0.5)
  if (any(x < 0))
    stop("'x' must be non-negative")
  N <- sum(x)
  if (is.null(size))
    size <- N
  else if (size != N)
    stop("size != sum(x), i.e. one is wrong")
  i0 <- prob == 0
  if (any(i0)) {
    if (any(x[i0] != 0))
      return(if (log) -Inf else 0)
    if (all(i0))
      return(if (log) 0 else 1)
    x <- x[!i0]
    prob <- prob[!i0]
  }
  
  
  r <- lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
  
  
  if (log)
    r
  else exp(r)
}






dDmultinom <- function (x, size , gamma, log = FALSE)
{

  gammaP <- sum(gamma)
  
  r <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(gammaP +1) - lgamma(size + gammaP +1) + sum(lgamma(x + gamma +1) - lgamma(gamma)) 

  
  if (log)
    r
  else exp(r)
}





###################################
# GLM
###################################




## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())

summary(glm.D93)







################################


pdf("rocX.pdf")
#ROC curve for score of prediction between SVM and NN
data(rocX.score)
r <- rocX(rocX.score$score, rocX.score$label)
plot(r)
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


dev.off()





















