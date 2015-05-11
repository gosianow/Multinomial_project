##############################################################################

# BioC 3.0
# Created 7 May 2015

# Check how changes the estimates of gamma+ and pi when reducing number of low expressed transcripts
# Take the real data expression from GEUVADIS analysis


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


keep <- 10


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

# dgeSQTL$tagwiseDispersion[paste(snp, gene, sep="-")] # 0.1092517







################################# run DM 
dgeDMList <- list()

modeList = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")
mode <- modeList[3]


modeDispList = c("optimize", "optim", "constrOptim", "grid")
modeDisp <- modeDispList[4]


dgeDM <- dmEstimateTagwiseDisp(dge, group = NULL, adjust = TRUE, mode = mode, epsilon = 1e-05, maxIte = 1000, modeDisp = modeDisp, interval = c(0, 1e+5), tol = 1e-08,  initDisp = 2, initWeirMoM = TRUE, gridLength = 10, gridRange = c(-6, 3), trend = c("none", "commonDispersion", "trendedDispersion")[1], priorDf = 10, span = 0.3, mcCores = 1, verbose = TRUE, plot = FALSE)
  
dgeDM$tagwiseDispersion




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


gamma0 <- seq(from = 0, to = 100, by = 4)[-1]

PLL <- rep(0, length(gamma0))
MSE <- matrix(0, nrow(y[[g]]), length(gamma0))
MSSE <- matrix(0, nrow(y[[g]]), length(gamma0))


for(i in 1:length(gamma0)){
  # i = 10
  
  PLL[i] <- dmAdjustedProfileLikTG(gamma0 = gamma0[i], y = y[[g]], ngroups=ngroups, lgroups=lgroups, igroups=igroups, adjust = TRUE, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose = FALSE)
  
  
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


# pdf(paste0(out.dir.s, "PLL_",mode, "_", modeDisp,"_keep", keep, ".pdf"), width = 12, height = 7)
# plot(gamma0, PLL, type = "l", col = "deeppink", lwd = 3)
# abline(v = dgeDM$tagwiseDispersion)
# dev.off()



pdf(paste0(out.dir.s, "PLL_MSE_", mode, "_", modeDisp,"_keep", keep, ".pdf"), width = 8, height = 5)

df <- data.frame(gamma0=gamma0, PLL = PLL)
pll <- ggplot(df, aes(gamma0, PLL)) + geom_line(colour = "deeppink") + geom_vline(xintercept = dgeDM$tagwiseDispersion)
print(pll)

df <- melt(MSE, varnames = c("ete_id", "gamma0"), value.name = "MSE")
mse <- ggplot(df, aes(gamma0, MSE, colour = ete_id)) + geom_line() + geom_vline(xintercept = dgeDM$tagwiseDispersion) 
print(mse)


df <- melt(MSSE, varnames = c("ete_id", "gamma0"), value.name = "MSSE")
mse <- ggplot(df, aes(gamma0, MSSE, colour = ete_id)) + geom_line() + geom_vline(xintercept = dgeDM$tagwiseDispersion) + stat_summary(fun.y = mean, colour="black", linetype = 2, geom="line", size = 1)
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


### Scenario parameters


### Check scenario
nBins <- 3
simPar <- list(s = "check", sample.size = 20, pi.org = rep(1, nBins)/nBins , g0.org = 100, nr.genes = 1e+04, nM = 150, tot = "uni")



### Like in GEUVADIS for snp_19_17269822-ENSG00000099331.7 - one with negative LR and two transcripts
simPar <- list(s = "GEUVADIS_snp_19_17269822", sample.size = 100, pi.org = c(0.4, 0.6) , g0.org = 7, nr.genes = 1e+04, nM = 2000, tot = "uni")



### Like in GEUVADIS for the SNP that have the most negative LR - snp_5_179056159-ENSG00000169045.13
piH <- c(0.4922335018, 0.1577883831, 0.0529355261, 0.0293984041, 0.0290847238, 0.0264224072, 0.0243318117, 0.0207018908, 0.0188256676, 0.0184684653, 0.0156071528, 0.0118260446, 0.0115589007, 0.0070278078, 0.0057829094, 0.0055859233, 0.0051052026, 0.0047757111, 0.0046757892, 0.0046102791, 0.0045847635, 0.0042218958, 0.0041134432, 0.0040331904, 0.0039526359, 0.0034897977, 0.0034433870, 0.0032155074, 0.0027587187, 0.0027284563, 0.0024924688, 0.0023142979, 0.0019838135, 0.0019837859, 0.0014639012, 0.0012791484, 0.0010303418, 0.0007425876, 0.0007168352, 0.0006505579, 0.0005237112, 0.0004658400, 0.0004015046, 0.0003565343, 0.0002020756, 0.0001042983)[1:10]

pi.org <- piH/sum(piH)

simPar <- list(s = "GEUVADIS_snp_5_179056159_10tr_df_g01", sample.size = 100, pi.org = pi.org, g0.org = 0.1, nr.genes = 1e+03, nM = 20000, tot = "uni")





mcCores <- 20




out.dir.s <- paste0(out.dir, "/", simPar$s, "/")
dir.create(out.dir.s, showWarnings=F, recursive=T)

sim <- simulate_from_DM(sample.size = simPar$sample.size, s = simPar$s, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nM, out.dir = out.dir.s , mc.cores=mcCores, save = FALSE)

sim$dge$samples$group <- as.factor(c(rep("C1", simPar$sample.size/2), rep("C2", simPar$sample.size/2)))

save(sim, file = paste0(out.dir.s, "/sim.RData"))



##################################################################################
#### Run DM pipeline, but do not estimate dispersion. Use true value as common dispersion
##################################################################################


### load simulation data
# load()



### run DM 
dgeDMList <- list()

modeList = c("constrOptim", "constrOptim2", "constrOptim2G", "optim2", "optim2NM", "FisherScoring")


### when using optim2: 
# Error in optim(par = piInit[-k], fn = dmLogLikkm1, gr = dmScoreFunkm1,  :
#                  L-BFGS-B needs finite values of 'fn'
#                In addition: Warning messages:
#                  1: In log(pi[i] * gamma0 + 1:y[i, j] - 1) : NaNs produced
#                2: In log(pi[i] * gamma0 + 1:y[i, j] - 1) : NaNs produced
### when using FisherScoring:
# * Negative piH: 0.5724435 0.1574793 0.0500953 0.04474325 0.0479419 0.01765365 0.06932008
# 0.04035913 0.0002681423 -0.0003042877
# piInit: 0.3635897 0.02963637 0.02267271 0.2320393 0.1595044 0.1328343 0.05497171 0.004189322
# 5.63362e-09 0.0005621153


mode <- modeList[3]

dge <- sim$dge

dgeDM <- dmFit(dge, group=NULL, dispersion = simPar$g0.org, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

dgeDM <- dmTest(dgeDM, mode = mode, epsilon = 1e-05, maxIte = 1000, verbose=FALSE, mcCores = mcCores)

cat("LR < 0 \n")
print(table(dgeDM$table$LR < 0))

dgeDMList[[mode]] <- dgeDM


pdf(paste0(out.dir.s, "/hist_", mode,".pdf"))

hist(dgeDM$table$PValue, breaks = 100, col = "#1E90FF")
hist(dgeDM$table$LR, breaks = 100, col = "#1E90FF")

dev.off()





save(dgeDMList, file = paste0(out.dir.s, "/dgeDMList.RData"))





### Compare constrOptim2 with constrOptim2G
names(dgeDMList)

table(dgeDMList[[1]]$table$LR == dgeDMList[[2]]$table$LR)

tab <- merge(dgeDMList[[1]]$table, dgeDMList[[2]]$table, by = "GeneID", suffixes = c("_co2g", "_co2"))

pdf(paste0(out.dir.s, "/hist_diffLR.pdf"))

hist(tab$LR_co2g - tab$LR_co2, breaks = 100, col = "#1E90FF")

dev.off()

tab[tab$LR_co2g < 0 | tab$LR_co2 < 0, ]


### LR < 0 for co2g
table(tab$LR_co2g < 0)
table(tab$LR_co2 < 0)

### more FP for co2
table(tab$FDR_co2g < 0.05)
table(tab$FDR_co2 < 0.05)







### Check the piH estimates for the genes where LR_co2g < 0

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)



### plot negative LR
tab <- tab[order(tab$LR_co2g, decreasing = FALSE), ]

### plot FP
tab <- tab[order(tab$PValue_co2g, decreasing = FALSE), ]



### only for one method
tab <- dgeDMList[[1]]$table
## negative LR
tab <- tab[order(tab$LR, decreasing = FALSE), ]
## FP
tab <- tab[order(tab$PValue, decreasing = FALSE), ]



genes2plot <- as.character(tab$GeneID[1:4])



mode <- c("constrOptim2G", "constrOptim2")[1]

dgeDM <- dgeDMList[[mode]]





pdf(paste0(out.dir.s, "/Proportions_", mode,"_FP.pdf"), width = 10, height = 5)


for(g in 1:length(genes2plot)){
  # g = 1
  gene <- genes2plot[g]
  # print(gene)
  Condition <- dgeDM$samples$group 
  expr <- dgeDM$counts[[gene]]
  # colnames(expr) <- metadata$SampleName
  rownames(expr) <- subset(dgeDM$genes, gene_id==gene)$ete_id    
  tot <- colSums(expr)
  labels <- strsplit2(rownames(expr), ":")[,2]
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
    theme(axis.text.x = element_text(angle = 20, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), legend.position="none") +
    ggtitle(paste0(gene)) +     
    geom_jitter(aes(fill = Condition, colour = factor(Condition, labels=c("C1b", "C2b"))), position = position_jitterdodge(dodge.width = 0.75), alpha = 0.5) +
    geom_boxplot(aes(colour = Condition), fill = "white", outlier.size = NA, alpha = 0) + 
    geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, fill = Samples), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 2, shape = 19, colour = "black") +
    geom_point(data = prop.est.null, aes(x = ete_id, y = Proportions), size = 3, shape = 18, colour = "orange") +
    scale_colour_manual(values=c("C1"="firebrick", "C2"="dodgerblue4", "C1b"="firebrick1", "C2b" = "dodgerblue"))  +
    coord_cartesian(ylim = c(-0.1, 1.1)) 
  
  
  print(ggb)
  
  
}


dev.off()





































