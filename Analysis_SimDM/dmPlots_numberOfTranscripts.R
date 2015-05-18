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



### the SNP that have the most negative LR
load("/home/Shared/data/seq/GEUVADIS/DM_0_1_3_sQTL_analysis/Results_Data_DM_0_1_3_TagwiseDisp_gridNone_tolEps/dgeSQTL_chr5.RData")

gene <- "ENSG00000145730.15"
snp <- "snp_5_102094652"





##################################################################################
# Prepare dge object
##################################################################################

out.dir.s <- paste0(out.dir,"/GEUVADIS_", snp, "_mostNegative/")
dir.create(out.dir.s, showWarnings=F, recursive=T)

out.dir.s <- paste0(out.dir.s, snp, "_")



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

dgeDM <- dmEstimateTagwiseDisp(dge, group = NULL, adjust = adjust, mode = mode, epsilon = 1e-03, maxIte = 1000, modeDisp = modeDisp, interval = c(0, 1e+5), tol = 1e-10,  initDisp = 2, initWeirMoM = TRUE, gridLength = 10, gridRange = c(-6, 6), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 10, span = 0.3, mcCores = 1, verbose = TRUE, plot = FALSE)
  

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




gamma0 <- sort(c(seq(from = 6, to = 8, by = 0.02)[-1], dgeDM$tagwiseDispersion), decreasing = FALSE)
plotPath <- paste0(out.dir.s, "PLL_MSE_", mode, "_", modeDisp,"_keep", keep, plotName, ".pdf")


PLL_MSE <- calculatePLL_MSE(dgeDM, gamma0, mode)
plotPLL_MSE(PLL_MSE, plotPath)




##################################################################################
################################# run DM
##################################################################################
# dgeDM$tagwiseDispersion[1] <- 1000



dgeDM <- dmFit(dgeDM, group=NULL, dispersion = "tagwiseDispersion", mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 1)

dgeDM <- dmTest(dgeDM, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 1)

dgeDM$table



dgeDMList[[mode]] <- dgeDM



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
  
  
  df.est <- data.frame(pi = c(dgeDM$fit[[g]]$piH[1, ], dgeDM$fit.null[[g]]$piH[1, ]), LL = c(dgeDM$fit[[g]]$logLik, dgeDM$fit.null[[g]]$logLik), Condition = factor(c(paste0("X", lgroups), "LLnull")), levels = c(c(paste0("X", lgroups), "LLnull")))
  
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































