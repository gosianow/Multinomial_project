# BioC 3.1

# Created 3 June 2015
# Compare the speed of bplapply and bpmapply 
# Updated 5 June 2015
# Use constrOptim for testing the speed


##############################################################################################################
# use BiocParallel
##############################################################################################################

library(BiocParallel)

registered()

bpparam()

BPPARAM <- SnowParam(workers = 5, type = "SOCK")

BPPARAM <- MulticoreParam(workers = 5)


funMean <- function(x){mean(x, na.rm = TRUE)}


bpmapply(function(y, funMean){
  yy <- funMean(y)
}, list(a = 1:10, b = 20:27), MoreArgs = list(funMean), BPPARAM = BPPARAM)



bpmapply(function(y){
  yy <- funMean(y)
}, list(a = 1:10, b = 20:27), BPPARAM = BPPARAM)




library(DM)


fun_bpmapply <- function(){
bpmapply(function(y){
  yy <- DM::rdirichlet(n = y, alpha = c(5, 5))
}, rep(10, 1000), BPPARAM = BPPARAM)
  return(NULL)
}


fun_bplapply <- function(){
bplapply(rep(10, 1000), function(y){
  yy <- DM::rdirichlet(n = y, alpha = c(5, 5))
},BPPARAM = BPPARAM)
return(NULL)
}


microbenchmark(fun_bpmapply(), fun_bplapply(), times=5)


##############################################################################################################
# load dgeSQTL data 
##############################################################################################################


setwd("/home/Shared/data/seq/GEUVADIS/")

library(DM)
library(BiocParallel)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)


BPPARAM <- MulticoreParam(workers=5)


######### run on DM_0_1_5 data 
out.dir <- "DM_0_1_5_sQTL_analysis/Results_Data_DM_0_1_5_TagwiseDisp_gridNone/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

load(paste0("DM_0_1_5_Data/dgeSQTL.RData"))

dgeSQTL_org <- dgeSQTL

### subset 100 genes
geneList <- names(dgeSQTL_org$counts)[1:1000]
o <- order(sapply(dgeSQTL_org$genotypes[geneList], nrow), decreasing = TRUE)
geneList <- geneList[o]
dgeSQTL$counts <- dgeSQTL_org$counts[geneList]
dgeSQTL$genotypes <- dgeSQTL_org$genotypes[geneList]




##############################################################################################################
# run elements of dmSQTLEstimateTagwiseDisp
##############################################################################################################

library(microbenchmark)


adjustDisp = TRUE; modeDisp = c("optimize", "optim", "constrOptim", "grid")[2]; intervalDisp = c(0, 1e+5); tolDisp = 1e-00;  initDisp = 10; initWeirMoMDisp = TRUE; gridLengthDisp = 15; gridRangeDisp = c(-7, 7); trendDisp = c("none", "commonDispersion", "trendedDispersion")[1]; priorDfDisp = 10; spanDisp = 0.3; modeProp = c( "constrOptim2", "constrOptim2G", "FisherScoring")[2]; tolProp = 1e-12; verbose = FALSE; plot = FALSE

### bpmapply

fun_bpmapply <- function(){
  
  dispList <- bpmapply(function(y, snps){
    # g = 1; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
    
    snps <- snps[1:min(nrow(snps), 100), , drop = FALSE]
    
    disp <- rep(NA, nrow(snps))
    names(disp) <- rownames(snps)
    
    for(i in 1:nrow(snps)){
      # i = 1
      
      NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
      yg <- y[, !NAs]             
      group <- snps[i, !NAs]
      group <- factor(group)
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      names(lgroups) <- lgroups
      igroups <- lapply(lgroups, function(x){which(group == x)})
      
      ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
      gamma0 <- initDisp
      if(is.null(dmAdjustedProfileLikTG(gamma0 = gamma0, y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose))){
        next
      }
      
      if(initWeirMoMDisp)
        initDisp <- weirMoM(data = yg, se=FALSE)
      
      optimum <- constrOptim(theta = initDisp, dmAdjustedProfileLikTG, grad = NULL, method = "Nelder-Mead",
                             ui=1, ci=1e-8, control = list(fnscale = -1, reltol = tolDisp), 
                             y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                             adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose )
      
      disp[i] <- optimum$par				 
      
    }
    
    return(disp)  
    
  }, dgeSQTL$counts,  dgeSQTL$genotypes, SIMPLIFY=FALSE, BPPARAM = BPPARAM)
  
  return(dispList)
}


### bplapply


fun_bplapply <- function(){
  
  dispList <- bplapply(geneList, function(g){
    # g = "ENSG00000052723.7"; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
    
    g = 1

    y = dgeSQTL$counts[[g]]
    snps = dgeSQTL$genotypes[[g]]
    
    snps <- snps[1:min(nrow(snps), 100), , drop = FALSE]
    
    disp <- rep(NA, nrow(snps))
    names(disp) <- rownames(snps)
    
    for(i in 1:nrow(snps)){
      # i = 1
      
      NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
      yg <- y[, !NAs]             
      group <- snps[i, !NAs]
      group <- factor(group)
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      names(lgroups) <- lgroups
      igroups <- lapply(lgroups, function(x){which(group == x)})
      
      
      ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
      gamma0 <- initDisp
      if(is.null(dmAdjustedProfileLikTG(gamma0 = gamma0, y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose))){
        next
      }
      
      if(initWeirMoMDisp)
        initDisp <- weirMoM(data = yg, se=FALSE)
      
      optimum <- constrOptim(theta = initDisp, dmAdjustedProfileLikTG, grad = NULL, method = "Nelder-Mead",
                             ui=1, ci=1e-8, control = list(fnscale = -1, reltol = tolDisp), 
                             y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                             adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose )
      
      disp[i] <- optimum$par					 
      
    }
    
    return(disp)  
    
  }, BPPARAM = BPPARAM)
  names(dispList) <- geneList
  return(dispList)
}


fun_bplapply2 <- function(){
  
  dispList <- bplapply(geneList, function(g, counts, genotypes, initDisp, adjustDisp, modeProp, tolProp, verbose, initWeirMoMDisp, tolDisp){
    # g = "ENSG00000052723.7"; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
    
    y = counts[[g]]
    snps = genotypes[[g]]
    
    snps <- snps[1:min(nrow(snps), 100), , drop = FALSE]
    
    disp <- rep(NA, nrow(snps))
    names(disp) <- rownames(snps)
    
    for(i in 1:nrow(snps)){
      # i = 1
      
      NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
      yg <- y[, !NAs]             
      group <- snps[i, !NAs]
      group <- factor(group)
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      names(lgroups) <- lgroups
      igroups <- lapply(lgroups, function(x){which(group == x)})
      
      ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
      gamma0 <- initDisp
      if(is.null(dmAdjustedProfileLikTG(gamma0 = gamma0, y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose))){
        next
      }
      
      if(initWeirMoMDisp)
        initDisp <- weirMoM(data = yg, se=FALSE)
      
      optimum <- constrOptim(theta = initDisp, dmAdjustedProfileLikTG, grad = NULL, method = "Nelder-Mead",
                             ui=1, ci=1e-8, control = list(fnscale = -1, reltol = tolDisp), 
                             y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                             adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose )
      
      disp[i] <- optimum$par					 
      
    }
    
    return(disp)  
    
  }, BPPARAM = BPPARAM, counts = dgeSQTL$counts, genotypes = dgeSQTL$genotypes, initDisp = initDisp, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose, initWeirMoMDisp = initWeirMoMDisp, tolDisp = tolDisp)
  names(dispList) <- geneList
  return(dispList)
}



### bpmapply with snow

fun_bpmapply_snow <- function(){
  
  dispList <- bpmapply(function(y, snps, initDisp, adjustDisp, modeProp, tolProp, verbose, initWeirMoMDisp, tolDisp){
    # g = 1; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
    
    snps <- snps[1:min(nrow(snps), 100), , drop = FALSE]
    
    disp <- rep(0, nrow(snps))
    names(disp) <- rownames(snps)
    
    for(i in 1:nrow(snps)){
      # i = 1
      
      NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
      yg <- y[, !NAs]             
      group <- snps[i, !NAs]
      group <- factor(group)
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      names(lgroups) <- lgroups
      igroups <- lapply(lgroups, function(x){which(group == x)})
      
      ### return NA if gene has 1 exon or observations in one sample in group (anyway this gene would not be fitted by dmFit)
      gamma0 <- initDisp
      if(is.null(DM::dmAdjustedProfileLikTG(gamma0 = gamma0, y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose))){
        disp[i] <- NA
        next
        
      }
      
      if(initWeirMoMDisp)
        initDisp <- DM::weirMoM(data = yg, se=FALSE)
      
      optimum <- constrOptim(theta = initDisp, DM::dmAdjustedProfileLikTG, grad = NULL, method = "Nelder-Mead",
                             ui=1, ci=1e-8, control = list(fnscale = -1, reltol = tolDisp), 
                             y = yg, ngroups=ngroups, lgroups=lgroups, igroups=igroups, 
                             adjustDisp = adjustDisp, modeProp = modeProp, tolProp = tolProp, verbose = verbose )
      
      disp[i] <- optimum$par				 
      
    }
    
    return(disp)  
    
  }, dgeSQTL$counts,  dgeSQTL$genotypes, MoreArgs = list(initDisp, adjustDisp, modeProp, tolProp, verbose, initWeirMoMDisp, tolDisp), SIMPLIFY=FALSE, BPPARAM = BPPARAM_snow)
  
  return(dispList)
}




timeM <- system.time(outM <- fun_bpmapply())
timeM

timeL <- system.time(outL <- fun_bplapply())
timeL

timeMS <- system.time(outM <- fun_bpmapply_snow())
timeMS


BPPARAM <- MulticoreParam(workers = 5)

microbenchmark(fun_bpmapply(), fun_bplapply(), times=1)



BPPARAM <- MulticoreParam(workers = 5)


microbenchmark(fun_bplapply(), fun_bplapply2(), times=1)



BPPARAM <- MulticoreParam(workers = 5)
BPPARAM_snow <- SnowParam(workers = 5, type = "SOCK")

microbenchmark(fun_bpmapply(), fun_bpmapply_snow(), times=1)





##############################################################################################################
# run elements of dmSQTLFit
##############################################################################################################


model = c("full", "null")[1]; dispersion = c("commonDispersion", "tagwiseDispersion")[1]; modeProp=c("constrOptim2", "constrOptim2G", "FisherScoring")[2]; tolProp = 1e-12; verbose=FALSE



### subset 100 genes
geneList <- names(dgeSQTL_org$counts)[1:1000]
o <- order(sapply(dgeSQTL_org$genotypes[geneList], nrow), decreasing = TRUE)
geneList <- geneList[o]
dgeSQTL$counts <- dgeSQTL_org$counts[geneList]
dgeSQTL$genotypes <- dgeSQTL_org$genotypes[geneList]


dgeSQTL$commonDispersion <- 4

gamma0 <- dgeSQTL$genotypes
gamma0 <- lapply(gamma0, function(g){
  # g = gamma0[[1]]
  disp <- rep(dgeSQTL$commonDispersion, nrow(g))
  names(disp) <- rownames(g)
  return(disp)
})




### bplapply


fun_bplapply <- function(){
  
  
  fit <- bplapply(geneList, function(g){
    # g = geneList[1]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
    
    y = dgeSQTL$counts[[g]]
    snps = dgeSQTL$genotypes[[g]]
    snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
    
    f <- list()
    
    for(i in 1:nrow(snps)){
      # i = 1
      if(is.na(gamma0[[g]][i]))
        f[[i]] <- NULL
      
      NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
      yg <- y[, !NAs]             
      group <- snps[i, !NAs]
      group <- factor(group)
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      names(lgroups) <- lgroups
      igroups <- lapply(lgroups, function(x){which(group == x)})
      
      
      f[[i]] <- dmOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                    gamma0 = gamma0[[g]][i], modeProp = modeProp, tolProp = tolProp, verbose = verbose)  
      
    }
    names(f) <- rownames(snps)
    
    return(f)
  }, BPPARAM = BPPARAM)
  
  names(fit) <- geneList
  
  return(fit)
}


fun_bplapply_numeric <- function(){

  
  fit <- bplapply(1:length(geneList), function(g){
    # g = geneList[1]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
    
    y = dgeSQTL$counts[[g]]
    snps = dgeSQTL$genotypes[[g]]
    snps <- snps[1:min(nrow(snps), 5), , drop = FALSE]
    
    f <- list()
    
    for(i in 1:nrow(snps)){
      # i = 1
      if(is.na(gamma0[[g]][i]))
        f[[i]] <- NULL
      
      NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
      yg <- y[, !NAs]             
      group <- snps[i, !NAs]
      group <- factor(group)
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      names(lgroups) <- lgroups
      igroups <- lapply(lgroups, function(x){which(group == x)})
      
      
      f[[i]] <- dmOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                    gamma0 = gamma0[[g]][i], modeProp = modeProp, tolProp = tolProp, verbose = verbose)  
      
    }
    names(f) <- rownames(snps)
    
    return(f)
  }, BPPARAM = BPPARAM)
  
  names(fit) <- geneList
  
  return(fit)
}



timeL <- system.time(outL <- fun_bplapply())
timeL



library(microbenchmark)


BPPARAM <- MulticoreParam(workers = 5)

microbenchmark(fun_bplapply(), fun_bplapply_numeric(), times=1)







##############################################################################################################
# check the size of objects: counts, genotypes
##############################################################################################################

c1 <- dgeSQTL_org$counts

c2 <- do.call(rbind, dgeSQTL_org$counts)

object.size(c1)
object.size(c2)



g0 <- dgeSQTL_org$genotypes

g1 <- lapply(names(dgeSQTL_org$genotypes), function(g){
  snps <- dgeSQTL_org$genotypes[[g]]
  rownames(snps) <- paste0(g, "-", rownames(snps))
  return(snps)
  })


g2 <- do.call(rbind, g1)


object.size(g0)
object.size(g1)
object.size(g2)




##############################################################################################################
# check the size of objects: fit list and array
##############################################################################################################



fun_fit_list <- function(){
  
  fit <- bplapply(geneList, function(g){
    # g = geneList[1]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
    
    # g = geneList[1]
    y = dgeSQTL$counts[[g]]
    snps = dgeSQTL$genotypes[[g]]
    snps <- snps[1:min(nrow(snps), 10), , drop = FALSE]
    
    f <- list()
    
    for(i in 1:nrow(snps)){
      # i = 1
      if(is.na(gamma0[[g]][i])){
        f[[i]] <- NULL
        next
      }

      NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
      yg <- y[, !NAs]             
      group <- snps[i, !NAs]
      group <- factor(group)
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      names(lgroups) <- lgroups
      igroups <- lapply(lgroups, function(x){which(group == x)})
      
      
      f[[i]] <- dmOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                    gamma0 = gamma0[[g]][i], modeProp = modeProp, tolProp = tolProp, verbose = verbose)  
      
    }
    names(f) <- rownames(snps)
    
    return(f)
  }, BPPARAM = BPPARAM)
  
  names(fit) <- geneList
  
  return(fit)
  
}


fun_fit_array <- function(){
  
  fit <- bplapply(geneList, function(g){
    # g = geneList[1]; y = dgeSQTL$counts[[g]]; snps = dgeSQTL$genotypes[[g]]
    
    # g = geneList[1]
    y = dgeSQTL$counts[[g]]
    snps = dgeSQTL$genotypes[[g]]
    snps <- snps[1:min(nrow(snps), 10), , drop = FALSE]
    
    ny <- nrow(y)
    ly <- rownames(y)
    nsnps <- nrow(snps)
    lsnps <- rownames(snps)
    
    LpiH <- array(NA, dim = c(ny, 3, nsnps), dimnames = list( transcripts = ly, groups = 0:2, snps = lsnps))
    Lgamma0 <- rep(NA, nsnps); names(Lgamma0) <- lsnps
    LlogLik <- matrix(NA, nsnps, 3, dimnames = list(lsnps, 0:2))
    Ldf <- matrix(NA, nsnps, 3, dimnames = list(lsnps, 0:2))
    
    for(i in 1:nrow(snps)){
      # i = 1
      if(is.na(gamma0[[g]][i]))
        next
      
      NAs <- is.na(snps[i, ]) | is.na(y[1, ])            
      yg <- y[, !NAs]             
      group <- snps[i, !NAs]
      group <- factor(group)
      ngroups <- nlevels(group)
      lgroups <- levels(group)
      names(lgroups) <- lgroups
      igroups <- lapply(lgroups, function(x){which(group == x)})
      
      
      fit <- dmOneGeneManyGroups(y = yg, ngroups = ngroups, lgroups = lgroups, igroups = igroups, 
                                 gamma0 = gamma0[[g]][i], modeProp = modeProp, tolProp = tolProp, verbose = verbose)
      
      if(is.null(fit))
        next
      
      LpiH[, lgroups , i] <- fit$piH
      Lgamma0[i] <- fit$gamma0
      LlogLik[i, lgroups] <- fit$logLik
      Ldf[i, lgroups] <- fit$df
      
    }
    
    return(list(piH = LpiH, gamma0 = Lgamma0, logLik = LlogLik, df = Ldf))
    
  }, BPPARAM = BPPARAM)
  
  
  names(fit) <- geneList
  
  return(fit)
  
}




BPPARAM <- MulticoreParam(workers = 10)


timeL <- system.time(fit_list <- fun_fit_list())
timeL

timeA <- system.time(fit_array <- fun_fit_array())
timeA


object.size(fit_list)

object.size(fit_array)


microbenchmark(fun_bplapply(), fun_bplapply_numeric(), times=1)




















piH <- array(0, dim = c(nrow(snps), nrow(y), 3))

piH_l <- lapply(1:nrow(snps), function(x){matrix(0, nrow(y), 3)})

object.size(piH)

object.size(piH_l)















































