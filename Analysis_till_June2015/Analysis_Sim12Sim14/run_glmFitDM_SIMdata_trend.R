
# continuation of run_glmFitDM_SIMdata.R 
# check the trend and two clouds



#######################################################
# trend: overall gene expression vs dispersion
#######################################################

rowSums(glmlrt$fit.full[[1]]$th)

glmlrt$fit.full[[1]]$Y

gh0.tmp <-  mclapply(names(glmlrt$fit.full), function(g){
  
  if(!is.null(glmlrt$fit.full[[g]])){           
    
    meanY <- mean(rowSums(glmlrt$fit.full[[g]]$Y))
    meanY0 <- mean(rowSums(glmlrt$fit.full[[g]]$Y[1:3,]))
    meanY1 <- mean(rowSums(glmlrt$fit.full[[g]]$Y[4:6,]))
    gh0.0f <- mean(glmlrt$fit.full[[g]]$gh0[1])
    gh0.1f <- mean(glmlrt$fit.full[[g]]$gh0[4])
    gh0.0n <- mean(glmlrt$fit.null[[g]]$gh0[1])
    
    return(data.frame(gh0.0f=gh0.0f, gh0.1f=gh0.1f, gh0.0n=gh0.0n, meanY=meanY, meanY0=meanY0, meanY1=meanY1))
  }    
  return(data.frame(gh0.0f=NA, gh0.1f=NA, gh0.0n=NA, meanY=NA, meanY0=NA, meanY1=NA))
  
}, mc.cores=10)


gh0 <- do.call(rbind, gh0.tmp)

gh0 <- data.frame(GeneID=names(glmlrt$fit.full), gh0)

head(gh0)



table.gh0 <- merge(gh0, table, by=1, all.x=TRUE, sort=FALSE)

table.gh0[!complete.cases(table.gh0), ]

table.gh0 <- table.gh0[complete.cases(table.gh0), ]
head(table.gh0)


table.gh0$status.est <- ifelse(table.gh0$FDR <= 0.05, 1, 0)




t.gh0.up <- table.gh0[table.gh0$gh0.0f > 1e+04, ]
t.gh0.dn <- table.gh0[table.gh0$gh0.0f <= 1e+04, ]



frac.up <- sum(t.gh0.up$status.est == 1 & t.gh0.up$status == 0) / sum(t.gh0.up$status.est == 1)

frac.dn <- sum(t.gh0.dn$status.est == 1 & t.gh0.dn$status == 0) / sum(t.gh0.dn$status.est == 1)





pdf(paste0("PLOTS/", name, "TREND_meanY01_vs_gh0.pdf"))

plot(table.gh0$meanY0, table.gh0$gh0.0f, col=2, log="xy", main="log(meanY) vs log(gh0)")
plot(table.gh0$meanY1 , table.gh0$gh0.1f, col=3, log="xy")
plot(table.gh0$meanY , table.gh0$gh0.0n, col=4, log="xy")


dev.off()




pdf(paste0("PLOTS/", name, "TREND_rpkm_vs_gh0.pdf"))

plot(table.gh0$rpkm, table.gh0$gh0.0f, col=2, log="xy", main="log(rpkm) vs log(gh0)")
plot(table.gh0$rpkm , table.gh0$gh0.1f, col=3, log="xy")
plot(table.gh0$rpkm , table.gh0$gh0.0n, col=4, log="xy")

plot(table.gh0$rpkm, table.gh0$gh0.0f, col=ifelse(table.gh0$status==1, 1, 2), log="xy", main="log(rpkm) vs log(gh0)")
plot(table.gh0$rpkm , table.gh0$gh0.1f, col=ifelse(table.gh0$status==1, 1, 3), log="xy")
plot(table.gh0$rpkm , table.gh0$gh0.0n, col=ifelse(table.gh0$status==1, 1, 4), log="xy")


dev.off()




pdf(paste0("PLOTS/", name, "TREND_meanY_vs_gh0.pdf"))

plot(table.gh0$meanY, table.gh0$gh0.0f, col=2, log="xy", main="log(meanY) vs log(gh0)")
plot(table.gh0$meanY , table.gh0$gh0.1f, col=3, log="xy")
plot(table.gh0$meanY , table.gh0$gh0.0n, col=4, log="xy")


plot(table.gh0$meanY, table.gh0$gh0.0f, col=ifelse(table.gh0$status==1, 1, 2), log="xy", main="log(meanY) vs log(gh0) \n black: status = 1")
plot(table.gh0$meanY , table.gh0$gh0.1f, col=ifelse(table.gh0$status==1, 1, 3), log="xy")
plot(table.gh0$meanY , table.gh0$gh0.0n, col=ifelse(table.gh0$status==1, 1, 4), log="xy")


plot(table.gh0$meanY, table.gh0$gh0.0f, col=ifelse(table.gh0$FDR <= 0.05, 1, 2), log="xy", main="log(meanY) vs log(gh0) \n black: FDR <= 0.05")
plot(table.gh0$meanY , table.gh0$gh0.1f, col=ifelse(table.gh0$FDR <= 0.05, 1, 3), log="xy")
plot(table.gh0$meanY , table.gh0$gh0.0n, col=ifelse(table.gh0$FDR <= 0.05, 1, 4), log="xy")


plot(table.gh0$meanY, table.gh0$gh0.0f, col=2, log="xy", pch=paste0(table.gh0$df), cex=0.5 ,main="log(meanY) vs log(gh0)")
plot(table.gh0$meanY , table.gh0$gh0.1f, col=3, log="xy", pch=paste0(table.gh0$df), cex=0.5)
plot(table.gh0$meanY , table.gh0$gh0.0n, col=4, log="xy", pch=paste0(table.gh0$df), cex=0.5)


plot(table.gh0$meanY, table.gh0$gh0.0f, col=table.gh0$df, log="xy", cex=0.5, main="log(meanY) vs log(gh0)")
plot(table.gh0$meanY , table.gh0$gh0.1f, col=table.gh0$df, log="xy", cex=0.5)
plot(table.gh0$meanY , table.gh0$gh0.0n, col=table.gh0$df, log="xy", cex=0.5)



dev.off()


pdf(paste0("PLOTS/", name, "HIST_gh0.pdf"))

hist(log10(table.gh0$gh0.0f), breaks=50)
hist(log10(table.gh0$gh0.1f), breaks=50)
hist(log10(table.gh0$gh0.0n), breaks=50)

dev.off()

#######################################################
# check the slice from two clouds plots
#######################################################


head(table.gh0)

slice <- table.gh0[table.gh0$meanY0 >= 100 & table.gh0$meanY0 <= 130, ]


slice.up <- slice[slice$gh0.0f > 1e+04, ]
slice.down <- slice[slice$gh0.0f < 1e+04, ]

table(slice.up$df)
table(slice.down$df)


gene.up <- slice.up[slice.up$df==2, "GeneID"][1]
gene.down <- slice.down[slice.down$df==2, "GeneID"][1]


gene.up <- slice.up[slice.up$df==2, "GeneID"][3]
gene.down <- slice.down[slice.down$df==2, "GeneID"][3]


glmlrt$fit.full[[gene.up]]
glmlrt$fit.full[[gene.down]]



save.image(paste0("PLOTS/", name, ".Rdata"))



pdf(paste0("PLOTS/", name, "CHECK_slice.pdf"))

dev.off()






#######################################################
# check
#######################################################

# optim + nlm

table.all.optim <- read.table("PLOTS/optim/table_all_MISO.xls", header=T)
head(table.all.optim)

table.all.nlm <- read.table("PLOTS/nlm/table_all_MISO.xls", header=T)
head(table.all.nlm)


table.all.all <- merge(table.all.optim, table.all.nlm, by=1)
write.table(table.all.all, paste0("PLOTS/",name,"table_all_on.xls"),quote=F, sep="\t", row.names=F, col.names=T)


# optim + g2


table.all.optim <- read.table("PLOTS/optim_p1_table_all_MISO.xls", header=T)
head(table.all.optim)

table.all.g2 <- read.table("PLOTS/G2_p1_table_all_MISO.xls", header=T)
head(table.all.g2)


table.all.all <- merge(table.all.optim, table.all.g2, by=1)
write.table(table.all.all, paste0("PLOTS/","merge_table_all_og2.xls"),quote=F, sep="\t", row.names=F, col.names=T)


#######################################################
# check for g2 & optim
#######################################################

gene <- "FBgn0038369"


glmlrt$fit.full[[gene]]
glmlrt$fit.null[[gene]]


gene <- "FBgn0035422"

g2lrt$fit.full[[gene]]
g2lrt$fit.null[[gene]]


#######################################################
# check for optim
#######################################################


# FP genes with very big LR
gene <- "FBgn0002283"
gene <- "FBgn0040284"


# FN
gene <- "FBgn0033232"
gene <- "FBgn0002626"


# LR < 0

# TN
gene <- "FBgn0025725"
gene <- "FBgn0032987"


# FN
gene <- "FBgn0031228"
gene <- "FBgn0259937"

gene <- "FBgn0031228"
gene <- "FBgn0036286"


glmlrt$fit.full[[gene]]
glmlrt$fit.null[[gene]]



#######################################################
# LR versus bh estimates
#######################################################

LRbh <-  mclapply(names(glmlrt$FitDM), function(g){
  
  if(!is.null(glmlrt$FitDM[[g]])){           
    LR <-  2*(glmlrt$FitDM[[g]]$loglikh - glmlrt$fit.null[[g]]$loglikh)
    
    bh0f <- mean(glmlrt$FitDM[[g]]$bh[,1])
    bh1f <- mean(glmlrt$FitDM[[g]]$bh[,2])
    bh0n <- mean(glmlrt$fit.null[[g]]$bh[,1])
    
    return(data.frame(LR=LR, bh0f=bh0f, bh1f=bh1f, bh0n=bh0n))
  }    
  return(data.frame(LR=NA, bh0f=NA, bh1f=NA, bh0n=NA))
  
}, mc.cores=10)


LRbh <- do.call(rbind, LRbh)
head(LRbh)


pdf(paste0("PLOTS/", name, "LR_vs_bh_MISO.pdf"))
plot(LRbh[,1], LRbh[,2], col=1)
# points(LRbh[,1], LRbh[,3], col=2)
# points(LRbh[,1], LRbh[,4], col=3)
dev.off()


#######################################################
# LR versus gh estimates
#######################################################

LRgh0 <-  mclapply(names(glmlrt$FitDM), function(g){
  
  if(!is.null(glmlrt$FitDM[[g]])){           
    LR <-  2*(glmlrt$FitDM[[g]]$loglikh - glmlrt$fit.null[[g]]$loglikh)
    
    gh0.0f <- mean(glmlrt$FitDM[[g]]$gh0[1])
    gh0.1f <- mean(glmlrt$FitDM[[g]]$gh0[4])
    gh0.0n <- mean(glmlrt$fit.null[[g]]$gh0[1])
    
    return(data.frame(LR=LR, gh0.0f=gh0.0f, gh0.1f=gh0.1f, gh0.0n=gh0.0n))
  }    
  return(data.frame(LR=NA, gh0.0f=NA, gh0.1f=NA, gh0.0n=NA))
  
}, mc.cores=10)


LRgh0 <- do.call(rbind, LRgh0)
head(LRgh0)

LRgh0$grp <- LRgh0$LR > 0


pdf(paste0("PLOTS/", name, "BOXPLOT_LR_vs_gh0_MISO.pdf"))

boxplot(gh0.0f ~ grp , LRgh0, col=2, log="y")
boxplot(gh0.1f ~ grp , LRgh0, col=3, log="y")
boxplot(gh0.0n ~ grp , LRgh0, col=4, log="y")

dev.off()





########################
# LogLik
########################

gene <- "FBgn0001291"

b <- glmlrt$FitDM[[gene]]$bh
Y <- glmlrt$FitDM[[gene]]$Y
X <- design

p <- ncol(X)
q <- ncol(Y)


g <- exp(X %*% t(b))  # n * q
gs <- rowSums(g)
ys <- rowSums(Y)

# formula (7) is equivalent to likelihood from dirmult
res <-   sum( lgamma(gs) - lgamma(ys+gs) + rowSums(lgamma(Y+g) - lgamma(g)) )
# normalization coeff
norm <- sum(lgamma(ys + 1) - rowSums(lgamma(Y + 1)))
#res <- res + norm

res


#########


ll <- function (x, t)
  log(t + x - 1)


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

#########


sum(unlist(lapply( 1:nrow(Y), function(i){
  
  loglik(x=Y[i, , drop=FALSE], t=g[i, , drop=FALSE])
  
} )))



sum(unlist(lapply( 1:nrow(glmlrt$FitDM[[gene]]$Y), function(i){
  
  loglik(x=glmlrt$FitDM[[gene]]$Y[i, , drop=FALSE], t=glmlrt$FitDM[[gene]]$gh[i, , drop=FALSE])
  
} )))



sum(unlist(lapply( 1:nrow(glmlrt$fit.null[[gene]]$Y), function(i){
  
  loglik(x=glmlrt$fit.null[[gene]]$Y[i, , drop=FALSE], t=glmlrt$fit.null[[gene]]$gh[i, , drop=FALSE])
  
} )))


loglik(x=glmlrt$fit.null[[gene]]$Y, t=glmlrt$fit.null[[gene]]$gh[1,, drop=FALSE])





