
# BioC 14


# install.packages("/home/gosia/R/packages/MGLM_0.0.6.tar.gz", lib="/home/gosia/R/libraries/3.1.0/")

library(MGLM)


n <- 200
d <- 4
alpha <- rep(1, d)
m <- 50
Y <- rdirm(m, alpha, n)
dmFit <- MGLMfit(Y, dist="DM")
print(dmFit)



dmFit$estimate

dmFit$logL













