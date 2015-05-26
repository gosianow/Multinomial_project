##############################################################################

# BioC 3.0
# Created 7 May 2015

# Check how changes the estimates of gamma+ and pi when reducing number of low expressed transcripts
# Take the real data expression from GEUVADIS analysis
# Produce PLL and MSE of piH plots

# Update 11 May 2015
# Simulate data with parameters from GEUVADIS real example and check the PLL and MSE of piH plots

# Update 17 May 2015 
# Plot likelihood versus PiH for genes with 2 transcripts

# Update 19 May 2015 
# Plot likelihood versus PiH for genes with 3 transcripts

# Update 20 May 2015 
# Plot score versus PiH for genes with 3 transcripts


##############################################################################

setwd("/home/gosia/Multinomial_project/Simulations_DM/")

out.dir <- "NumberOfTranscripts/"
dir.create(out.dir, showWarnings=F, recursive=T)


library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)

#### function to simulate data from DM
source("/home/gosia/R/R_Multinomial_project/Analysis_SimDM/simulate_from_DM.R")




### Source all R files in DM package
Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])


library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)



##################################################################################
# Get the real data from GEUVADIS results
##################################################################################


### the SNP that have the most negative LR - snp_5_179056159-ENSG00000169045.13
load("/home/Shared/data/seq/GEUVADIS/DM_0_1_2_sQTL_analysis/Results_Data_DM_0_1_2_TagwiseDisp_gridCommonDispersion/dgeSQTL_chr5.RData")

gene <- "ENSG00000169045.13"
snp <- "snp_5_179056159"



### the SNP that have the most negative LR  with 2 transcripts - ENSG00000112893.5 − snp_5_109101871
load("/home/Shared/data/seq/GEUVADIS/DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone/dgeSQTL_chr5.RData")

gene <- "ENSG00000112893.5"
snp <- "snp_5_109101871"


gene <- "ENSG00000164574.11"
snp <- "snp_5_153569041"



### the SNP that have the most negative LR  with 3 transcripts -
load("/home/Shared/data/seq/GEUVADIS/DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone/dgeSQTL_chr5.RData")

gene <- "ENSG00000170571.6"
snp <- "snp_5_49699894"


# null proportions are more weard
gene <- "ENSG00000181904.8"
snp <- "snp_5_134183242"

out.dir.s <- paste0(out.dir,"/GEUVADIS_", snp, "_3trans/")
dir.create(out.dir.s, showWarnings=F, recursive=T)
out.dir.s <- paste0(out.dir.s, snp, "_")



### the SNP that have the most negative LR
load("/home/Shared/data/seq/GEUVADIS/DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tolEps/dgeSQTL_chr5.RData")

gene <- "ENSG00000145730.15"
snp <- "snp_5_102094652"


out.dir.s <- paste0(out.dir,"/GEUVADIS_", snp, "_mostNegative/")
dir.create(out.dir.s, showWarnings=F, recursive=T)
out.dir.s <- paste0(out.dir.s, snp, "_")



##################################################################################
# Prepare dge object
##################################################################################




y <- dgeSQTL$counts[[gene]]

### filter out the transcripts with low proportions
keep <- nrow(y)

prop <- sort(rowSums(y)/sum(y), decreasing = TRUE)
counts <- y[rownames(y) %in% names(prop[1:keep]), ]

rownames(dgeSQTL$genes) <- dgeSQTL$genes$ete_id
genes <- dgeSQTL$genes[rownames(counts), ]

genotype <- dgeSQTL$genotypes[dgeSQTL$SNPs$SNP_id == snp & dgeSQTL$SNPs$gene_id == gene, ]
table(genotype, useNA = "always")

nas <- is.na(genotype) | is.na(counts[1, ])

counts <- counts[, !nas]
genotype <- genotype[!nas]



### create the DGE object
dge <- DGEList( counts=counts , group = genotype , genes =  genes)
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels=unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR


dge$commonDispersion <- dgeSQTL$commonDispersion


dgeSQTL$commonDispersion
dgeSQTL$tagwiseDispersion[paste(snp, gene, sep="-")]







##################################################################################
################################# run DM  tagwiseDispersion
##################################################################################

plotName <- ""
# plotName <- "_tol10"
plotName <- "_tolEps"

dgeDMList <- list()



modeList = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2G", "optim2NM", "FisherScoring")
mode <- modeList[2]


modeDispList = c("optimize", "optim", "constrOptim", "grid")
modeDisp <- modeDispList[4]

adjust = TRUE

dgeDM <- dmEstimateTagwiseDisp(dge, group = NULL, adjust = adjust, mode = mode, epsilon = 1e-03, maxIte = 1000, modeDisp = modeDisp, interval = c(0, 1e+5), tol = .Machine$double.eps,  initDisp = 2, initWeirMoM = TRUE, gridLength = 10, gridRange = c(-6, 6), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 10, span = 0.3, mcCores = 1, verbose = TRUE, plot = FALSE)
  

dgeDM$tagwiseDispersion



##################################################################################
################################# plot PLL and MSE of piH
##################################################################################

calculatePLL_MSE <- function(dgeDM, gamma0, mode){
  
  g = 1
  
  tagwiseDispersion <- dgeDM$tagwiseDispersion[g]
  
  y <- dgeDM$counts[[g]]
  trans <- rownames(y)
  ngenes <- length(y)
  
  group <- dgeDM$samples$group
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  nlibs <- length(group)
  
  igroups <- list()
  for(gr in 1:ngroups){
    # gr=2
    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
  }
  
  
  PLL <- rep(0, length(gamma0))
  PLLnull <- rep(0, length(gamma0))
  PLLadj <- rep(0, length(gamma0))
  MSE <- matrix(0, nrow(y), length(gamma0))
  rownames(MSE) <- trans
  colnames(MSE) <- gamma0
  MSSE <- matrix(0, nrow(y), length(gamma0))
  rownames(MSSE) <- trans
  colnames(MSSE) <- gamma0
  
  PLLconditions <-  matrix(0, ngroups, length(gamma0))
  rownames(PLLconditions) <- lgroups
  colnames(PLLconditions) <- gamma0
  
  for(i in 1:length(gamma0)){
    # i = 1
    
    f <- dmOneGeneManyGroups(y = y, ngroups=ngroups, lgroups=lgroups, igroups=igroups, gamma0=gamma0[i], mode = mode, epsilon = 1e-05, maxIte = 1000, verbose = FALSE)
    
    PLLconditions[, i] <- f$logLik 
    
#     PLL[i] <- dmAdjustedProfileLikTG(gamma0 = gamma0[i], y = y, ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = FALSE, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose = FALSE)
    PLL[i] <- sum(f$logLik)
    
    
    igroups.null <- list()
    igroups.null[[lgroups[1]]] <- sort(unlist(igroups))
    
    PLLnull[i] <- dmAdjustedProfileLikTG(gamma0 = gamma0[i], y = y, ngroups=1, lgroups=lgroups[1], igroups=igroups.null, adjust = FALSE, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose = FALSE)
    

adj <- dmAdjCROneGeneManyGroups(y = y , ngroups = ngroups, lgroups = lgroups, igroups = igroups, gamma0 = gamma0[i], piH = f$piH)
if(adj == Inf)
  PLLadj[i] <- NA

    PLLadj[i] <- sum(f$logLik) - adj
  


    expr <- y 
    prop.smp <- prop.table(expr, 2)
    SE <- matrix(0, nrow(prop.smp), ncol(prop.smp))
    SSE <- matrix(0, nrow(prop.smp), ncol(prop.smp))
    
    for(gr in 1:ngroups){
      # gr = 1
      ps = prop.smp[, igroups[[gr]], drop = FALSE]
      
      ps.est <- f$piH[, lgroups[gr]]
      
      ### only difference
      SE[, igroups[[gr]]] <- abs(sweep(ps, 1, ps.est, FUN = "-"))
      
      ### standardized
      SSE[, igroups[[gr]]] <- sweep(abs(sweep(ps, 1, ps.est, FUN = "-")), 1, ps.est, FUN = "/")    
      
    }
    
    MSE[, i] <- rowMeans(SE)

    MSSE[, i] <- rowMeans(SSE)

    
  }
  
  
  return(list(gamma0 = gamma0, PLLconditions = PLLconditions, PLL = PLL, PLLnull = PLLnull, PLLadj = PLLadj, MSE = MSE, MSSE = MSSE, tagwiseDispersion = tagwiseDispersion))
  
}


plotPLL_MSE <- function(PLL_MSE, plotPath){
  
  gamma0 <- PLL_MSE$gamma0
  PLL <- PLL_MSE$PLL
  PLLconditions <- PLL_MSE$PLLconditions
  PLLadj <- PLL_MSE$PLLadj
  PLLnull <- PLL_MSE$PLLnull
  tagwiseDispersion <- PLL_MSE$tagwiseDispersion
  MSE <- PLL_MSE$MSE
  MSSE <- PLL_MSE$MSSE
  
  
  pdf(plotPath, width = 8, height = 5)

  
  df <- data.frame(gamma0=gamma0, PLL = PLL, PLLadj = PLLadj, PLLnull = PLLnull)
  df <- melt(df, id.vars = "gamma0", variable.name = "Method", value.name = "PLL")
  
  pll <- ggplot(df, aes(gamma0, PLL, group = Method, colour = Method)) + geom_line() + geom_vline(xintercept = tagwiseDispersion, linetype = 2)
  print(pll)
  
  dfs <- subset(df, gamma0 < (tagwiseDispersion + 1) & gamma0 > (tagwiseDispersion - 1) & Method != "PLLadj")
  pll <- ggplot(dfs, aes(gamma0, PLL, group = Method, colour = Method)) + geom_line() + geom_vline(xintercept = tagwiseDispersion, linetype = 2)
  print(pll)
  
  dfs <- subset(df, gamma0 < (tagwiseDispersion + 1) & gamma0 > (tagwiseDispersion - 1) & Method == "PLLadj")
  pll <- ggplot(dfs, aes(gamma0, PLL, group = Method, colour = Method)) + geom_line() + geom_vline(xintercept = tagwiseDispersion, linetype = 2)
  print(pll)
  
  
  df <- melt(PLLconditions, varnames = c("Condition", "gamma0"), value.name = "PLL")
  df$Condition <- factor(df$Condition)
  
  pll <- ggplot(df, aes(gamma0, PLL, group = Condition, colour = Condition)) + geom_line() + geom_vline(xintercept = tagwiseDispersion, linetype = 2) + facet_grid(Condition ~ ., scales = "free_y")
  print(pll)
  
  
  df <- melt(MSE, varnames = c("ete_id", "gamma0"), value.name = "MSE")
  mse <- ggplot(df, aes(gamma0, MSE, colour = ete_id)) + geom_line() + geom_vline(xintercept = tagwiseDispersion, linetype = 2) 
  print(mse)
  
  
  df <- melt(MSSE, varnames = c("ete_id", "gamma0"), value.name = "MSSE")
  mse <- ggplot(df, aes(gamma0, MSSE, colour = ete_id)) + geom_line() + geom_vline(xintercept = tagwiseDispersion, linetype = 2) + stat_summary(fun.y = mean, colour="black", linetype = 2, geom="line", size = 1)
  print(mse)
  
  
  dev.off() 
  
  
  
}



gridLength = 10
gridRange = c(-6, 6)
splinePts <- seq(from = gridRange[1], to = gridRange[2], length = gridLength)
splineDisp <- dgeDM$commonDispersion * 2^splinePts


gamma0 <- sort(c(seq(from = 6, to = 8, by = 0.02)[-1], dgeDM$tagwiseDispersion), decreasing = FALSE)
plotPath <- paste0(out.dir.s, "PLL_MSE_", mode, "_", modeDisp,"_keep", keep, plotName, ".pdf")


PLL_MSE <- calculatePLL_MSE(dgeDM, gamma0, mode)

# PLL_MSE$tagwiseDispersion <- dgeDM$tagwiseDispersion
plotPLL_MSE(PLL_MSE, plotPath)




##################################################################################
################################# run DM
##################################################################################
# dgeDM$tagwiseDispersion[1] <- 1000


# modeList = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2G", "optim2NM", "FisherScoring")
# mode <- modeList[3]


dgeDM <- dmFit(dgeDM, group=NULL, dispersion = "tagwiseDispersion", mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 1)

dgeDM <- dmTest(dgeDM, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 1)

dgeDM$table



dgeDMList[[mode]] <- dgeDM

names(dgeDMList)

##################################################################################
################################# plot LL of piH for genes with 2 transcripts
##################################################################################





plotLL <- function(dgeDM, piX, gamma0, mode, plotPath){
  
  g = 1
  y <- dgeDM$counts[[g]]
  transcripts <- rownames(y)
  genes <- names(y)
  ngenes <- length(y)
  
  group <- dgeDM$samples$group
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  igroups <- list()
  for(gr in 1:ngroups){
    # gr=2
    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
  }
   
  LL <- matrix(0, length(piX), ngroups)
  colnames(LL) <- lgroups
  LLnull <- rep(0, length(piX))
  
  for(i in 1:length(piX)){
    # i = 1    
    switch(mode, 
           
           constrOptim2 = {
             for(gr in 1:ngroups){                 
               LL[i, gr] <- dmLogLikkm1(piX[i], gamma0, y[, igroups[[gr]], drop = FALSE])               
             }             
             LLnull[i] <- dmLogLikkm1(piX[i], gamma0, y)            
           },
           
           constrOptim2G = {
             for(gr in 1:ngroups){                
               LL[i, gr] <- dmLogLikGkm1(piX[i], gamma0, y[, igroups[[gr]], drop = FALSE])               
             }             
             LLnull[i] <- dmLogLikGkm1(piX[i], gamma0, y)             
           })
       
  }
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  colour <- c(gg_color_hue(ngroups), "orange")
  
  
  df.est <- data.frame(pi = c(dgeDM$fit[[g]]$piH[1, ], dgeDM$fit.null[[g]]$piH[1, ]), LL = c(dgeDM$fit[[g]]$logLik, dgeDM$fit.null[[g]]$logLik), Condition = factor(c(paste0("X", lgroups), "LLnull"), levels = c(paste0("X", lgroups), "LLnull") ))
  
  df <- data.frame(pi = piX, LL, LLnull = LLnull)
  df <- melt(df, id.vars = "pi", variable.name = "Condition", value.name = "LL")
  
  pll <- ggplot(df, aes(pi, LL, group = Condition, colour = Condition)) + 
    geom_line() +
    ggtitle(paste0(transcripts[1])) +
    geom_point(data = df.est, aes(pi, LL, colour = Condition), size = 3, alpha = 0.6, shape = 18) +
    facet_grid(Condition ~ ., scales = "free_y") +
    scale_colour_manual(values = colour)

  pdf(plotPath, width = 8, height = 10)
  print(pll) 
  dev.off() 
  
}






piX <- seq(0.2, 0.4, by = 0.001)
# piX <- seq(0, 1, by = 0.1)
piX <- piX[c(-1, -length(piX))]
gamma0 <- dgeDM$tagwiseDispersion

plotPath <- paste0(out.dir.s, "LLpiH_", mode, "_", modeDisp,"_keep", keep, plotName, ".pdf")

plotLL(dgeDM, piX, gamma0, mode, plotPath)




##################################################################################
################################# plot LL of piH for genes with 3 transcripts
##################################################################################



library(plotly)


calculateLL <- function(dgeDM, piXY, gamma0, mode){
  
  g = 1
  y <- dgeDM$counts[[g]]
  transcripts <- rownames(y)
  genes <- names(y)
  ngenes <- length(y)
  
  group <- dgeDM$samples$group
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  igroups <- list()
  for(gr in 1:ngroups){
    # gr=2
    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
  }
  
  LL <- matrix(0, nrow(piXY), ngroups)
  colnames(LL) <- lgroups
  LLnull <- rep(0, nrow(piXY))
  
  for(i in 1:nrow(piXY)){
    # i = 1    
    
    if(sum(piXY[i, ]) >= 1 - sqrt(.Machine$double.eps)){
      LL[i, ] <- NA
      LLnull[i] <- NA
    }else{
      switch(mode, 
             
             constrOptim2 = {
               for(gr in 1:ngroups){                 
                 LL[i, gr] <- dmLogLikkm1(as.numeric(piXY[i, ]), gamma0, y[, igroups[[gr]], drop = FALSE])               
               }             
               LLnull[i] <- dmLogLikkm1(as.numeric(piXY[i, ]), gamma0, y)            
             },
             
             constrOptim2G = {
               for(gr in 1:ngroups){                
                 LL[i, gr] <- dmLogLikGkm1(as.numeric(piXY[i, ]), gamma0, y[, igroups[[gr]], drop = FALSE])               
               }             
               LLnull[i] <- dmLogLikGkm1(as.numeric(piXY[i, ]), gamma0, y)             
             })
    }
    
  }
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  colour <- c(gg_color_hue(ngroups), "orange")
  
  
  df.est <- data.frame(x = c(dgeDM$fit[[g]]$piH[1, ], dgeDM$fit.null[[g]]$piH[1, ]), y = c(dgeDM$fit[[g]]$piH[2, ], dgeDM$fit.null[[g]]$piH[2, ]), LL = c(dgeDM$fit[[g]]$logLik, dgeDM$fit.null[[g]]$logLik), Condition = factor(c(paste0("X", lgroups), "LLnull"), levels = c(paste0("X", lgroups), "LLnull") ))
  
  df <- data.frame(piXY, LL, LLnull = LLnull)
  
  df <- melt(df, id.vars = c("x", "y"), variable.name = "Condition", value.name = "LL")
  
  
  return(list(df = df, df.est = df.est, colour = colour, gamma0 = gamma0))
  
}


plotLL <- function(LL, plotPath, range = 0.75){
  
 df <- LL$df
 df.est <- LL$df.est
 colour <- LL$colour
 gamma0 <- LL$gamma0
  
# py <- plotly(username="gosianow", key="yy902bfs7g")  # open plotly connection
# plotly_url <- rep(NA, nrow(df.est))

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


pdf(plotPath, width = 8, height = 6)

for(i in 1:nrow(df.est)){
  # i = 1
  Condition <- df.est$Condition[i]
  
  df.est.tmp <- df.est[df.est$Condition == Condition,]
  df.tmp <- df[df$Condition == Condition,] 
  
  
  
  pll <- ggplot(df.tmp, aes(x=x, y=y, fill = LL)) + 
    geom_tile() +
    xlab(paste0(transcripts[1])) +
    ylab(paste0(transcripts[2])) +
    ggtitle(paste0("Condition ", Condition, "\n gamma0 = ", gamma0)) +
    geom_point(data = df.est.tmp, aes(x = x, y = y), fill = colour[i], size = 4, alpha = 1, shape = 23) +
    scale_fill_gradientn(colours = myPalette(200), limits = c(min(df.tmp$LL, na.rm = TRUE) + diff(range(df.tmp$LL, na.rm = TRUE))*range, max(df.tmp$LL, na.rm = TRUE)), na.value = myPalette(200)[1])
    # scale_fill_gradientn(colours = myPalette(200), limits = c(min(df.tmp$LL, na.rm = TRUE) + diff(range(df.tmp$LL, na.rm = TRUE))/2, max(df.tmp$LL, na.rm = TRUE)))
    # scale_fill_gradientn(colours = myPalette(200))

  print(pll)
  
  
  
  ### for plotly can not have NAs 
#   df.tmp[is.na(df.tmp$LL), "LL"] <- min(df.tmp$LL, na.rm = TRUE)
#   
#   pll <- ggplot(df.tmp, aes(x=x, y=y, fill = LL)) + 
#     geom_tile() +
#     xlab(paste0(transcripts[1])) +
#     ylab(paste0(transcripts[2])) +
#     ggtitle(paste0("Condition ", Condition))
#   
#   out <- py$ggplotly(pll, kwargs=list(filename=paste0("LLpiH_", i), fileopt="overwrite"))
#   plotly_url[i] <- out$response$url
  
}

dev.off() 
  


}




piX <- seq(0, 1, by = 0.1)
piX <- piX[c(-1, -length(piX))]
piY <- seq(0, 1, by = 0.1)
piY <- piY[c(-1, -length(piY))]

piXY <- expand.grid(x=piX, y=piY)  
# piXY <- piXY[rowSums(piXY) < 1, ]


gamma0 <- dgeDM$tagwiseDispersion

plotPath <- paste0(out.dir.s, "LLpiH_", mode, "_", modeDisp,"_keep", keep, plotName, ".pdf")

LL <- calculateLL(dgeDM, piXY, gamma0, mode)

plotLL(LL, plotPath)





#### plot LL for "constrOptim2"  and "constrOptim2G"



piX <- seq(0.1, 0.3, by = 0.001)
piX <- piX[c(-1, -length(piX))]
piY <- seq(0.6, 0.8, by = 0.001)
piY <- piY[c(-1, -length(piY))]
piXY <- expand.grid(x=piX, y=piY)  

gamma0 <- dgeDMList[["constrOptim2"]]$tagwiseDispersion



### constrOptim2

LL2 <- calculateLL(dgeDMList[["constrOptim2"]], piXY, gamma0, mode = "constrOptim2")

plotPath <- paste0(out.dir.s, "LLpiH_", mode = "constrOptim2", "_", modeDisp,"_keep", keep, plotName, ".pdf")

plotLL(LL2, plotPath)


### constrOptim2G

LL2G <- calculateLL(dgeDMList[["constrOptim2G"]], piXY, gamma0, mode = "constrOptim2G")

plotPath <- paste0(out.dir.s, "LLpiH_", mode = "constrOptim2G", "_", modeDisp,"_keep", keep, plotName, ".pdf")

plotLL(LL2G, plotPath)



## plot the difference between LL
LLdiff <- LL2G
LLdiff[["df"]][, "LL"] <- LL2G[["df"]][, "LL"] - LL2[["df"]][, "LL"]
LLdiff[["df.est"]][, "LL"] <- NA


plotPath <- paste0(out.dir.s, "LLdifference", modeDisp,"_keep", keep, plotName, ".pdf")

plotLL(LLdiff, plotPath, range = 0)


##################################################################################
################################# plot score of piH for genes with 3 transcripts
##################################################################################
## !!!!!!!!  there was a mistake in dmScoreFunGkm1 because repRow was not used


library(plotly)


calculateSC <- function(dgeDM, piXY, gamma0, mode){
  
  g = 1
  y <- dgeDM$counts[[g]]
  transcripts <- rownames(y)
  genes <- names(y)
  ngenes <- length(y)
  
  group <- dgeDM$samples$group
  group <- as.factor(group)
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  
  igroups <- list()
  for(gr in 1:ngroups){
    # gr=2
    igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
  }
  
  SCx <- matrix(0, nrow(piXY), ngroups)
  colnames(SC1) <- lgroups
  SCy <- matrix(0, nrow(piXY), ngroups)
  colnames(SC2) <- lgroups
  SCnullx <- rep(0, nrow(piXY))
  SCnully <- rep(0, nrow(piXY))
  
  
  for(i in 1:nrow(piXY)){
    # i = 29900    
    
    if(sum(piXY[i, ]) >= 1 - sqrt(.Machine$double.eps)){
      SCx[i, ] <- NA
      SCy[i, ] <- NA
      SCnullx[i] <- NA
      SCnully[i] <- NA
    }else{
      switch(mode, 
             
             constrOptim2 = {
               for(gr in 1:ngroups){                 
                 sc <- dmScoreFunkm1(as.numeric(piXY[i, ]), gamma0, y[, igroups[[gr]], drop = FALSE])    
                 SCx[i, gr] <- sc[1]
                 SCy[i, gr] <- sc[2]
               }             
               sc <- dmScoreFunkm1(as.numeric(piXY[i, ]), gamma0, y)        
               SCnullx[i] <- sc[1]
               SCnully[i] <- sc[2]
             },
             
             constrOptim2G = {
               for(gr in 1:ngroups){
                 # gr = 1
                 sc <- dmScoreFunGkm1(as.numeric(piXY[i, ]), gamma0, y[, igroups[[gr]], drop = FALSE])    
                 SCx[i, gr] <- sc[1]
                 SCy[i, gr] <- sc[2]
               }             
               sc <- dmScoreFunGkm1(as.numeric(piXY[i, ]), gamma0, y)        
               SCnullx[i] <- sc[1]
               SCnully[i] <- sc[2]           
             })
    }
    
  }
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  colour <- c(gg_color_hue(ngroups), "orange")
  
  
  df.est <- data.frame(x = c(dgeDM$fit[[g]]$piH[1, ], dgeDM$fit.null[[g]]$piH[1, ]), y = c(dgeDM$fit[[g]]$piH[2, ], dgeDM$fit.null[[g]]$piH[2, ]), Condition = factor(c(paste0("X", lgroups), "null"), levels = c(paste0("X", lgroups), "null") ), SC = NA)
  
  df1 <- data.frame(piXY, SCx, null = SCnullx)  
  df1 <- melt(df1, id.vars = c("x", "y"), variable.name = "Condition", value.name = "SC")
  df2 <- data.frame(piXY, SCy, null = SCnully)  
  df2 <- melt(df2, id.vars = c("x", "y"), variable.name = "Condition", value.name = "SC")
  
  return(list(df1 = df1, df2 = df2, df.est = df.est, colour = colour, gamma0 = gamma0))
  
}


plotSC <- function(SC, plotPath, range = 1, limit = NULL){
  
  df1 <- SC$df1
  df2 <- SC$df2
  df.est <- SC$df.est
  colour <- SC$colour
  gamma0 <- SC$gamma0

  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  
  
  pdf(plotPath, width = 8, height = 6)
  
  for(i in 1:nrow(df.est)){
    # i = 1
    Condition <- df.est$Condition[i]
    
    df.est.tmp <- df.est[df.est$Condition == Condition,]
    df1.tmp <- df1[df1$Condition == Condition,]    
    if(is.null(limit))
      limit <- max(abs(df1.tmp$SC), na.rm= TRUE)*range
    
    sc <- ggplot(df1.tmp, aes(x=x, y=y, fill = SC)) + 
      geom_tile() +
      xlab(paste0(transcripts[1])) +
      ylab(paste0(transcripts[2])) +
      ggtitle(paste0("Condition ", Condition, " pi1 \n gamma0 = ", gamma0)) +
      geom_point(data = df.est.tmp, aes(x = x, y = y), fill = colour[i], size = 4, alpha = 1, shape = 23) +
      # scale_fill_gradientn(colours = myPalette(200))
      # scale_fill_gradientn(colours = myPalette(200), limits = c(- max(df1$SC, na.rm= TRUE), max(df1$SC, na.rm= TRUE)))
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", limits = c(-limit , limit))
    # scale_fill_gradient2(low = "red", mid = "white", high = "blue")
    
    
    print(sc)

    df2.tmp <- df2[df2$Condition == Condition,] 
    if(is.null(limit))
    limit <- max(abs(df2.tmp$SC), na.rm= TRUE)*range
    
    
    sc <- ggplot(df2.tmp, aes(x=x, y=y, fill = SC)) + 
      geom_tile() +
      xlab(paste0(transcripts[1])) +
      ylab(paste0(transcripts[2])) +
      ggtitle(paste0("Condition ", Condition, " pi2 \n gamma0 = ", gamma0)) +
      geom_point(data = df.est.tmp, aes(x = x, y = y), fill = colour[i], size = 4, alpha = 1, shape = 23) +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", limits = c(- limit, limit))
      # scale_fill_gradient2(low = "red", mid = "white", high = "blue")
    
    print(sc)
    
    
  }
  
  dev.off() 
  
  
  
}




piX <- seq(0, 1, by = 0.1)
piX <- piX[c(-1, -length(piX))]
piY <- seq(0, 1, by = 0.1)
piY <- piY[c(-1, -length(piY))]

piXY <- expand.grid(x=piX, y=piY)  
# piXY <- piXY[rowSums(piXY) < 1, ]


gamma0 <- dgeDM$tagwiseDispersion

plotPath <- paste0(out.dir.s, "SCOREpiH_", mode, "_", modeDisp,"_keep", keep, plotName, ".pdf")



SC <- calculateSC(dgeDM, piXY, gamma0, mode)

plotSC(SC, plotPath)




#### plot SCORE for "constrOptim2"  and "constrOptim2G"



piX <- seq(0.1, 0.3, by = 0.001)
piX <- piX[c(-1, -length(piX))]
piY <- seq(0.6, 0.8, by = 0.001)
piY <- piY[c(-1, -length(piY))]
piXY <- expand.grid(x=piX, y=piY)  

gamma0 <- dgeDMList[["constrOptim2"]]$tagwiseDispersion



### constrOptim2

SC2 <- calculateSC(dgeDMList[["constrOptim2"]], piXY, gamma0, mode = "constrOptim2")

plotPath <- paste0(out.dir.s, "SCOREpiH_", mode = "constrOptim2", "_", modeDisp,"_keep", keep, plotName, "2.pdf")

plotSC(SC2, plotPath, range = 0.003, limit = 20)


### constrOptim2G

SC2G <- calculateSC(dgeDMList[["constrOptim2G"]], piXY, gamma0, mode = "constrOptim2G")

plotPath <- paste0(out.dir.s, "SCOREpiH_", mode = "constrOptim2G", "_", modeDisp,"_keep", keep, plotName, "2.pdf")

plotSC(SC2G, plotPath, range = 0.003, limit = 20)



## plot the difference between SCORES
SCdiff <- SC2G
SCdiff[["df1"]][, "SC"] <- SC2G[["df1"]][, "SC"] - SC2[["df1"]][, "SC"]
SCdiff[["df2"]][, "SC"] <- SC2G[["df2"]][, "SC"] - SC2[["df2"]][, "SC"]
SCdiff[["df.est"]][, "SC"] <- NA


plotPath <- paste0(out.dir.s, "SCOREdifference", modeDisp,"_keep", keep, plotName, "2.pdf")

plotSC(SCdiff, plotPath, range = 0.05)









sc2x <- SC2[["df1"]]
sc2gx <- SC2G[["df1"]]


sc2x <- SC2[["df1"]][complete.cases(SC2[["df1"]]), ]
sc2gx <- SC2G[["df1"]][complete.cases(SC2G[["df1"]]), ]



sc2x <- sc2x[sc2x$Condition == "X0", ]
sc2gx <- sc2gx[sc2gx$Condition == "X0", ]


min(abs(sc2x$SC))
min(abs(sc2gx$SC))

sc2x[which.min(abs(sc2x$SC)), ]

sc2gx[which.min(abs(sc2gx$SC)), ]









##################################################################################
################################# plot proportions
##################################################################################
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


plotProportions <- function(dgeDM, genes2plot, plotPath){
  
  
  pdf(plotPath, width = 12, height = 7)
  
  
  for(g in 1:length(genes2plot)){
    # g = 1
    
    gene <- genes2plot[g]
    # print(gene)
    Condition <- dgeDM$samples$group 
    expr <- dgeDM$counts[[gene]]
    # colnames(expr) <- metadata$SampleName
    rownames(expr) <- subset(dgeDM$genes, gene_id==gene)$ete_id    
    
    tot <- colSums(expr)
    labels <- rownames(expr)
    prop.smp <- data.frame( ete_id =  labels, t(apply(expr, 1, function(t){ t / tot })))  
    n <- nrow(expr)  
    prop.est <- data.frame(ete_id = labels, dgeDM$fit[[gene]]$piH)
    colnames(prop.est) <- c("ete_id", colnames(dgeDM$fit[[gene]]$piH))
    prop.est.null <- data.frame(ete_id = labels, dgeDM$fit.null[[gene]]$piH)
    
    
    prop.smp.m <- melt(prop.smp, id.vars = "ete_id", variable.name = "Samples", value.name = "Proportions")  
    prop.smp.m$ete_id <- factor(prop.smp.m$ete_id, levels = unique(prop.smp.m$ete_id))
    prop.smp.m$Samples <- factor(prop.smp.m$Samples)
    prop.smp.m$Condition <- rep(Condition, each = nrow(prop.smp))
    
    prop.est.m <- melt(prop.est, id.vars = "ete_id", variable.name = "Condition", value.name = "Proportions")
    prop.est.m$ete_id <- factor(prop.est.m$ete_id, levels = unique(prop.est.m$ete_id))
    prop.est.m$Condition <- factor(prop.est.m$Condition)
    
    colnames(prop.est.null) <- c("ete_id", "Proportions")
    
    
    ### box plots with points2
    ggb <- ggplot(prop.smp.m, aes(x = ete_id, y = Proportions)) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), legend.position="none", plot.title = element_text(size=10)) +
      ggtitle(paste0(gene, "\n TagwiseDispersion = ", dgeDM$tagwiseDispersion[gene], "\n LR = ", dgeDM$table[dgeDM$table$GeneID == gene, "LR"], "\n PValue = ", dgeDM$table[dgeDM$table$GeneID == gene, "PValue"])) +     
      geom_jitter(aes(fill = Condition, colour = factor(Condition)), position = position_jitterdodge(dodge.width = 0.75), alpha = 0.5) +
      geom_boxplot(aes(colour = Condition), fill = "white", outlier.size = NA, alpha = 0) + 
#       geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, fill = Samples), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 2, shape = 18, colour = "black") +
#       geom_point(data = prop.est.null, aes(x = ete_id, y = Proportions), size = 3, shape = 18, colour = "orange") +
      geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, fill = Condition), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 3, shape = 23) +
      geom_point(data = prop.est.null, aes(x = ete_id, y = Proportions), size = 3, shape = 23, fill = "orange") +
      coord_cartesian(ylim = c(-0.1, 1.1)) 
    
    
    print(ggb)
    
    
  }
  
  
  dev.off()
  
  
  
  
}





genes2plot <- gene

plotPath <- paste0(out.dir.s, "Proportions_",mode, "_", modeDisp,"_keep", keep, plotName, ".pdf")

plotProportions(dgeDM, genes2plot, plotPath)

  


 










########################################################################################################################

# Simulate data from two group null distribution with common dispersion

########################################################################################################################

### Scenario
simPar <- list(s = paste0("GEUVADIS_", snp, "_",mode, "_", modeDisp,"_keep", keep), sample.size = nrow(dgeDM$samples), pi.org = dgeDM$fit.null[[gene]]$piH , g0.org = dgeDM$fit.null[[gene]]$gamma0, nr.genes = 1e+03, nM = dgeDM$samples$lib.size, tot = "fix")



mcCores <- 20

out.dir.s <- paste0(out.dir, "/", simPar$s, "/")
dir.create(out.dir.s, showWarnings=F, recursive=T)


sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nM, out.dir = out.dir.s , mc.cores=mcCores, save = FALSE)

sim$dge$samples$group <- factor(dgeDM$samples$group)
sim$dge$commonDispersion <- dgeDM$fit.null[[gene]]$gamma0

# save(sim, file = paste0(out.dir.s, "/sim.RData"))



##################################################################################
#### Run DM pipeline
##################################################################################


### load simulation data
# load()


################################# run DM 
dgeDMList <- list()

modeList = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")
mode <- modeList[2]


modeDispList = c("optimize", "optim", "constrOptim", "grid")
modeDisp <- modeDispList[4]


dge <- sim$dge



dgeDM <- dmEstimateTagwiseDisp(dge, group = NULL, adjust = TRUE, mode = mode, epsilon = 1e-05, maxIte = 1000, modeDisp = modeDisp, interval = c(0, 1e+5), tol = 1e-08,  initDisp = 2, initWeirMoM = TRUE, gridLength = 10, gridRange = c(-6, 3), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 10, span = 0.3, mcCores = mcCores, verbose = FALSE, plot = FALSE)

# dgeDM$tagwiseDispersion



pdf(paste0(out.dir.s, "gamma0_", mode, "_", modeDisp,"_keep", keep, ".pdf"), width = 7, height = 7)

df <- data.frame(tagwiseDispersion = dgeDM$tagwiseDispersion, meanExpr = as.character(round(dgeDM$meanExpr)))
gg <- ggplot(data = df, aes(meanExpr, tagwiseDispersion)) + geom_point(position = "jitter", colour = "mediumslateblue") + geom_hline(yintercept = dgeDM$commonDispersion)
print(gg)

dev.off()





################################# run DM


# dgeDM$tagwiseDispersion[1] <- 1000

dgeDM <- dmFit(dgeDM, group=NULL, dispersion = "tagwiseDispersion", mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

head(dgeDM$table)


dgeDMList[[mode]] <- dgeDM



pdf(paste0(out.dir.s, "HistPValuesLR_", mode, "_", modeDisp,"_keep", keep, ".pdf"))

hist(dgeDM$table$PValue, breaks = 100, col = "#1E90FF")
hist(dgeDM$table$LR, breaks = 100, col = "#1E90FF")

dev.off()






### Compare constrOptim2 with constrOptim2G
names(dgeDMList)

table(dgeDMList[[1]]$table$LR == dgeDMList[[2]]$table$LR)

tab <- merge(dgeDMList[[1]]$table, dgeDMList[[2]]$table, by = "GeneID", suffixes = c("_co2g", "_co2"))


pdf(paste0(out.dir.s, "HistPValuesLRdiff_", modeDisp,"_keep", keep, ".pdf"))
hist(tab$LR_co2g - tab$LR_co2, breaks = 100, col = "#1E90FF")
dev.off()





################################# plot proportions

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


tab <- dgeDM$table

## negative LR
tab <- tab[order(tab$LR, decreasing = FALSE), ]
genes2plot <- as.character(tab$GeneID[1:4])


plotPath <- paste0(out.dir.s, "Proportions_",mode, "_", modeDisp,"_keep", keep, "NegLR.pdf")

plotProportions(dgeDM, genes2plot, plotPath)


## FP
tab <- tab[order(tab$PValue, decreasing = FALSE), ]
genes2plot <- as.character(tab$GeneID[1:4])

plotPath <- paste0(out.dir.s, "Proportions_",mode, "_", modeDisp,"_keep", keep, "FP.pdf")

plotProportions(dgeDM, genes2plot, plotPath)




################################# plot PLL and MSE of piH


modeList = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")
mode <- modeList[2]



y <- dge$counts
genes <- names(y)
ngenes <- length(y)



group <- dge$samples$group
group <- as.factor(group)
ngroups <- nlevels(group)
lgroups <- levels(group)
nlibs <- length(group)

igroups <- list()
for(gr in 1:ngroups){
  # gr=2
  igroups[[lgroups[gr]]] <- which(group == lgroups[gr])
}

g = "g470"


gamma0 <- seq(from = 1, to = 50, by = 0.2)[-1]

PLL <- rep(0, length(gamma0))
PLLnull <- rep(0, length(gamma0))
PLLadj <- rep(0, length(gamma0))
MSE <- matrix(0, nrow(y[[g]]), length(gamma0))
MSSE <- matrix(0, nrow(y[[g]]), length(gamma0))


for(i in 1:length(gamma0)){
  # i = 10
  
  PLL[i] <- dmAdjustedProfileLikTG(gamma0 = gamma0[i], y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = FALSE, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose = FALSE)
  
  igroups.null <- list()
  igroups.null[[lgroups[1]]] <- sort(unlist(igroups))
  
  PLLnull[i] <- dmAdjustedProfileLikTG(gamma0 = gamma0[i], y = y[[g]], ngroups=1, lgroups=lgroups[1], igroups=igroups.null, adjust = FALSE, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose = FALSE)
  
  PLLadj[i] <- dmAdjustedProfileLikTG(gamma0 = gamma0[i], y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = TRUE, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose = FALSE)
  
  f <- dmOneGeneManyGroups(y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, gamma0=gamma0[i], mode = mode, epsilon = 1e-05, maxIte = 1000, verbose = FALSE)
  
  expr <- y[[g]]  
  prop.smp <- prop.table(expr, 2)
  SE <- matrix(0, nrow(prop.smp), ncol(prop.smp))
  SSE <- matrix(0, nrow(prop.smp), ncol(prop.smp))
  
  for(gr in 1:ngroups){
    # gr = 1
    ps = prop.smp[, igroups[[gr]], drop = FALSE]
    
    ps.est <- f$piH[, lgroups[gr]]
    
    ### only difference
    SE[, igroups[[gr]]] <- abs(sweep(ps, 1, ps.est, FUN = "-"))
    
    ### standardized
    SSE[, igroups[[gr]]] <- sweep(abs(sweep(ps, 1, ps.est, FUN = "-")), 1, ps.est, FUN = "/")    
    
  }
  
  MSE[, i] <- rowMeans(SE)
  rownames(MSE) <- rownames(ps)
  colnames(MSE) <- gamma0
  MSSE[, i] <- rowMeans(SSE)
  rownames(MSSE) <- rownames(ps)
  colnames(MSSE) <- gamma0
  
}




pdf(paste0(out.dir.s, "PLL_MSE_", mode, "_", modeDisp,"_keep", keep, ".pdf"), width = 8, height = 5)

df <- data.frame(gamma0=gamma0, PLL = PLL, PLLadj = PLLadj, PLLnull = PLLnull)
df <- melt(df, id.vars = "gamma0", variable.name = "Method", value.name = "PLL")

pll <- ggplot(df, aes(gamma0, PLL, group = Method, colour = Method)) + geom_line() + geom_vline(xintercept = dgeDM$tagwiseDispersion[g])
print(pll)


pll <- ggplot(df[df$gamma0 < 10 & df$gamma0 > 5, ], aes(gamma0, PLL, group = Method, colour = Method)) + geom_line() + geom_vline(xintercept = dgeDM$tagwiseDispersion[g])
print(pll)

df <- melt(MSE, varnames = c("ete_id", "gamma0"), value.name = "MSE")
mse <- ggplot(df, aes(gamma0, MSE, colour = ete_id)) + geom_line() + geom_vline(xintercept = dgeDM$tagwiseDispersion[g]) 
print(mse)


df <- melt(MSSE, varnames = c("ete_id", "gamma0"), value.name = "MSSE")
mse <- ggplot(df, aes(gamma0, MSSE, colour = ete_id)) + geom_line() + geom_vline(xintercept = dgeDM$tagwiseDispersion[g]) + stat_summary(fun.y = mean, colour="black", linetype = 2, geom="line", size = 1)
print(mse)


dev.off()































