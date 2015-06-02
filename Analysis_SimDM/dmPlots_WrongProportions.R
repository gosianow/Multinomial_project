##############################################################################

# BioC 3.0
# Created 1 June 2015

# Check the estimation of proportions with DM_0_1_4 for GEUVADIS data DM_0_1_2_Data_clean
# Use the SNPs for a gene with the highest number of significant SNPs ENSG00000196735.6 


##############################################################################

setwd("/home/gosia/Multinomial_project/Simulations_DM/")

out.dir <- "WrongProportions/"
dir.create(out.dir, showWarnings=F, recursive=T)


library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)


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

### SNP that is significant for a gene with highest number of significant SNPs
load("/home/Shared/data/seq/GEUVADIS/DM_0_1_4_sQTL_analysis/Results_Data_DM_0_1_2_clean_TagwiseDisp_gridNone_tol12_constrOptim2G/dgeSQTL_chr6.RData")

gene <- "ENSG00000196735.6"
snp <- "snp_6_32601361"


out.dir.s <- paste0(out.dir,"/GEUVADIS_", snp, "_least_significant/")
dir.create(out.dir.s, showWarnings=F, recursive=T)
out.dir.s <- paste0(out.dir.s, snp, "_")



##################################################################################
# Prepare dge object
##################################################################################




y <- dgeSQTL$counts[[gene]]

### filter out the transcripts with low proportions
## method 1
# keep <- nrow(y)
# prop <- sort(rowSums(y)/sum(y), decreasing = TRUE)
# counts <- y[rownames(y) %in% names(prop[1:keep]), ]


## method 2
min_samps <- 5
min_transcript_prop <- 0.1 # in cpm

prop <- prop.table(y, 2)
trans2keep <- rowSums(prop > min_transcript_prop) > min_samps
trans2keep

keep <- sum(trans2keep)

counts <- y[trans2keep, ]


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

plotName <- "_original"
# plotName <- "_tol10"
# plotName <- "_tolEps"

dgeDMList <- list()



modeList = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2G", "optim2NM", "FisherScoring")
mode <- modeList[3]


modeDispList = c("optimize", "optim", "constrOptim", "grid")
modeDisp <- modeDispList[4]

adjust = TRUE

dgeDM <- dmEstimateTagwiseDisp(dge, group = NULL, adjust = adjust, mode = mode, epsilon = 1e-05, maxIte = 1000, modeDisp = modeDisp, interval = c(0, 1e+5), tol = 1e-08,  initDisp = 2, initWeirMoM = TRUE, gridLength = 10, gridRange = c(-6, 6), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 10, span = 0.3, mcCores = 1, verbose = TRUE, plot = FALSE)
  

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
splineDisp


gamma0 <- sort(c(seq(from = 3, to = 5, by = 0.02)[-1], dgeDM$tagwiseDispersion), decreasing = FALSE)
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

rownames(dgeDM$table) <- dgeDM$table$GeneID
dgeDM$table


dgeDMList[[mode]] <- dgeDM

names(dgeDMList)





##################################################################################
################################# plot proportions
##################################################################################
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)




plotProportions <- function(dgeDM, genes2plot, plotPath){
  
  
  pdf(plotPath, width = 11, height = 6)
  
  
  for(g in 1:length(genes2plot)){
    # g = 1
    
    gene <- genes2plot[g]
    # print(gene)
    Condition <- dgeDM$samples$group 
    expr <- dgeDM$counts[[gene]]
    
    labels <- rownames(expr)
    
    prop.smp <- data.frame( ete_id =  labels, prop.table(expr, 2))  
    prop.est <- data.frame(ete_id = labels, dgeDM$fit[[gene]]$piH)
    colnames(prop.est) <- c("ete_id", colnames(dgeDM$fit[[gene]]$piH))
    
    prop.est.null <- data.frame(ete_id = labels, dgeDM$fit.null[[gene]]$piH)
    
    #### order transcipts by decreasing proportions 
    order.tr <- labels[order(apply(aggregate(t(prop.smp[, -1]), by = list(Condition = Condition), median)[, -1], 2, max), decreasing = TRUE)]   
    
    prop.smp.m <- melt(prop.smp, id.vars = "ete_id", variable.name = "Samples", value.name = "Proportions")  
    prop.smp.m$ete_id <- factor(prop.smp.m$ete_id, levels = order.tr)
    prop.smp.m$Samples <- factor(prop.smp.m$Samples)
    prop.smp.m$Condition <- rep(Condition, each = nrow(prop.smp))
    
    prop.est.m <- melt(prop.est, id.vars = "ete_id", variable.name = "Condition", value.name = "Proportions")
    prop.est.m$ete_id <- factor(prop.est.m$ete_id, levels = order.tr)
    prop.est.m$Condition <- factor(prop.est.m$Condition)
    prop.est.m$Samples <- factor(prop.est.m$Condition)
    
    colnames(prop.est.null) <- c("ete_id", "Proportions")
    prop.est.null$ete_id <- factor(prop.est.null$ete_id, levels = order.tr)
    prop.est.null$Samples <- "Null"
    
    main <- paste0(gene, "\n Mean Expression = ", round(dgeDM$meanExpr[gene]), " / Dispersion = ", round(dgeDM$fit.null[[gene]]$gamma0, 2), "\n LR = ", round(dgeDM$table[gene, "LR"], 4) , " / P-value = ", sprintf("%.02e", dgeDM$table[gene, "PValue"]), " / FDR = ", sprintf("%.02e", dgeDM$table[gene, "FDR"]))
    
    
    ### box plots with points2
    ggb <- ggplot(prop.smp.m, aes(x = ete_id, y = Proportions)) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), legend.position="none", plot.title = element_text(size=10)) +
      ggtitle(main) +     
      geom_jitter(aes(fill = Condition, colour = factor(Condition)), position = position_jitterdodge(dodge.width = 0.75), alpha = 0.5) +
      geom_boxplot(aes(colour = Condition), fill = "white", outlier.size = NA, alpha = 0) + 
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

  


 
































