##############################################################################

# BioC 3.0
# Created 7 May 2015

# Check how changes the estimates of gamma+ and pi when reducing number of low expressed transcripts
# Take the real data expression from GEUVADIS analysis
# Produce PLL and MSE of piH plots

# Update 11 May 2015
# Simulate data with parameters from GEUVADIS real example and check the PLL and MSE of piH plots

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
load("/home/Shared/data/seq/GEUVADIS/DM_0.1.2_sQTL_analysis/Results_TagwiseDisp_gridCommonDispersion/dgeSQTL_chr5.RData")



gene <- "ENSG00000169045.13"
snp <- "snp_5_179056159"

out.dir.s <- paste0(out.dir, "GEUVADIS_", snp)

### keep two genes because "grid" did not work for one gene
# counts <- dgeSQTL$counts[c(gene, "ENSG00000113048.11")]
# counts <- do.call(rbind, counts)
# genes <- dgeSQTL$genes[dgeSQTL$genes$gene_id %in%  c(gene, "ENSG00000113048.11"), ]


counts <- dgeSQTL$counts[c(gene)]
# counts <- do.call(rbind, counts)


keep <- 48 ### max 48


counts_filtered <- lapply(counts, function(y){
  # y = counts[[1]]
  
  prop <- sort(rowSums(y)/sum(y), decreasing = TRUE)
  
  y <- y[rownames(y) %in% names(prop[1:keep]), ]
  
  ### add 1 count when 0s
  # y[ y==0 ] <- 1
  
  return(y)
})

names(counts_filtered) <- names(counts)



counts <- do.call(rbind, counts_filtered)

genes <- dgeSQTL$genes[dgeSQTL$genes$gene_id %in%  c(gene), ]
genes <- genes[genes$ete_id %in% rownames(counts), ]

genes


genotype <- dgeSQTL$genotypes[dgeSQTL$SNPs$SNP_id == snp & dgeSQTL$SNPs$gene_id == gene, ]
table(genotype)



### create the DGE object
dge <- DGEList( counts=counts , group = genotype , genes =  genes)
dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels=unique(dge$genes$gene_id)))
dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR

dge$commonDispersion <- dgeSQTL$commonDispersion #3.996774

# dgeSQTL$tagwiseDispersion[paste(snp, gene, sep="-")] # 0.1092517 - this is the values estimated with DM_0.1.2.







################################# run DM 
dgeDMList <- list()

modeList = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")
mode <- modeList[3]


modeDispList = c("optimize", "optim", "constrOptim", "grid")
modeDisp <- modeDispList[4]

adjust = TRUE

dgeDM <- dmEstimateTagwiseDisp(dge, group = NULL, adjust = adjust, mode = mode, epsilon = 1e-05, maxIte = 1000, modeDisp = modeDisp, interval = c(0, 1e+5), tol = 1e-08,  initDisp = 2, initWeirMoM = TRUE, gridLength = 10, gridRange = c(-6, 3), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 10, span = 0.3, mcCores = 1, verbose = TRUE, plot = FALSE)
  
dgeDM$tagwiseDispersion


# modeDisp <- "Orig_g01"
# dgeDM$tagwiseDispersion[1] <- dgeSQTL$tagwiseDispersion[paste(snp, gene, sep="-")]



################################# plot PLL and MSE of piH

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

g = 1


gamma0 <- seq(from = 0, to = 5, by = 0.05)[-1]

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


pll <- ggplot(df[df$gamma0 < (dgeDM$tagwiseDispersion[g] + 2) & df$gamma0 > (dgeDM$tagwiseDispersion[g] - 2), ], aes(gamma0, PLL, group = Method, colour = Method)) + geom_line() + geom_vline(xintercept = dgeDM$tagwiseDispersion[g])
print(pll)

df <- melt(MSE, varnames = c("ete_id", "gamma0"), value.name = "MSE")
mse <- ggplot(df, aes(gamma0, MSE, colour = ete_id)) + geom_line() + geom_vline(xintercept = dgeDM$tagwiseDispersion[g]) 
print(mse)


df <- melt(MSSE, varnames = c("ete_id", "gamma0"), value.name = "MSSE")
mse <- ggplot(df, aes(gamma0, MSSE, colour = ete_id)) + geom_line() + geom_vline(xintercept = dgeDM$tagwiseDispersion[g]) + stat_summary(fun.y = mean, colour="black", linetype = 2, geom="line", size = 1)
print(mse)


dev.off()



################################# run DM


# dgeDM$tagwiseDispersion[1] <- 1000

dgeDM <- dmFit(dgeDM, group=NULL, dispersion = "tagwiseDispersion", mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 1)

dgeDM <- dmTest(dgeDM, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = 1)

dgeDM$table



dgeDMList[[mode]] <- dgeDM





################################# plot proportions

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)



genes2plot <- gene




plotPath <- paste0(out.dir.s, "Proportions_",mode, "_", modeDisp,"_keep", keep, ".pdf")

plotProportions(dgeDM, genes2plot, plotPath)

  


  
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
    prop.est.null <- data.frame(ete_id = labels, dgeDM$fit.null[[gene]]$piH)
    
    
    prop.smp.m <- melt(prop.smp, id.vars = "ete_id", variable.name = "Samples", value.name = "Proportions")  
    prop.smp.m$ete_id <- factor(prop.smp.m$ete_id, levels = unique(prop.smp.m$ete_id))
    prop.smp.m$Samples <- factor(prop.smp.m$Samples)
    prop.smp.m$Condition <- rep(Condition, each = nrow(prop.smp))
    
    prop.est.m <- melt(prop.est, id.vars = "ete_id", variable.name = "Samples", value.name = "Proportions")
    prop.est.m$ete_id <- factor(prop.est.m$ete_id, levels = unique(prop.est.m$ete_id))
    prop.est.m$Samples <- factor(prop.est.m$Samples)
    
    colnames(prop.est.null) <- c("ete_id", "Proportions")
    
    
    ### box plots with points2
    ggb <- ggplot(prop.smp.m, aes(x = ete_id, y = Proportions)) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), legend.position="none", plot.title = element_text(size=10)) +
      ggtitle(paste0(gene, "\n TagwiseDispersion = ", dgeDM$tagwiseDispersion[gene], "\n LR = ", dgeDM$table[dgeDM$table$GeneID == gene, "LR"], "\n PValue = ", dgeDM$table[dgeDM$table$GeneID == gene, "PValue"])) +     
      geom_jitter(aes(fill = Condition, colour = factor(Condition)), position = position_jitterdodge(dodge.width = 0.75), alpha = 0.5) +
      geom_boxplot(aes(colour = Condition), fill = "white", outlier.size = NA, alpha = 0) + 
      geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, fill = Samples), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 2, shape = 19, colour = "black") +
      geom_point(data = prop.est.null, aes(x = ete_id, y = Proportions), size = 3, shape = 18, colour = "orange") +
      # scale_colour_manual(values=c("C1"="firebrick", "C2"="dodgerblue4", "C1b"="firebrick1", "C2b" = "dodgerblue"))  +
      coord_cartesian(ylim = c(-0.1, 1.1)) 
    
    
    print(ggb)
    
    
  }
  
  
  dev.off()
  
  
  
  
}

  



##################################################################################
# Simulate data from two group null distribution with common dispersion
##################################################################################

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































