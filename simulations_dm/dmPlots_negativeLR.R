##############################################################################

# BioC 3.0
# Created 23 Apr 2015

# Investigate why the LR is sometimes negative; Compare different mode
# Update 16 May 2015
# Simulate case with 2 transcripts with params that will lead to negative LR




##############################################################################

setwd("/home/gosia/Multinomial_project/Simulations_DM/")

out.dir <- "NegativeLR/"
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
# Simulate data from two group null distribution with common dispersion

##################################################################################


### Scenario parameters


######### Check scenario
nBins <- 3
simPar <- list(s = "check", sample.size = 20, pi.org = rep(1, nBins)/nBins , g0.org = 100, nr.genes = 1e+04, nM = 150, tot = "uni")



######### Like in GEUVADIS for snp_19_17269822-ENSG00000099331.7 - one with negative LR and two transcripts
simPar <- list(s = "GEUVADIS_snp_19_17269822", sample.size = 100, pi.org = c(0.4, 0.6) , g0.org = 7, nr.genes = 1e+04, nM = 2000, tot = "uni")



######### Like in GEUVADIS for the SNP that have the most negative LR - snp_5_179056159-ENSG00000169045.13
piH <- c(0.4922335018, 0.1577883831, 0.0529355261, 0.0293984041, 0.0290847238, 0.0264224072, 0.0243318117, 0.0207018908, 0.0188256676, 0.0184684653, 0.0156071528, 0.0118260446, 0.0115589007, 0.0070278078, 0.0057829094, 0.0055859233, 0.0051052026, 0.0047757111, 0.0046757892, 0.0046102791, 0.0045847635, 0.0042218958, 0.0041134432, 0.0040331904, 0.0039526359, 0.0034897977, 0.0034433870, 0.0032155074, 0.0027587187, 0.0027284563, 0.0024924688, 0.0023142979, 0.0019838135, 0.0019837859, 0.0014639012, 0.0012791484, 0.0010303418, 0.0007425876, 0.0007168352, 0.0006505579, 0.0005237112, 0.0004658400, 0.0004015046, 0.0003565343, 0.0002020756, 0.0001042983)[1:10]

pi.org <- piH/sum(piH)

simPar <- list(s = "GEUVADIS_snp_5_179056159_10tr_df_g01", sample.size = 100, pi.org = pi.org, g0.org = 0.1092517, nr.genes = 1e+03, nM = 20000, tot = "uni")




######### Scenario with 2 transcripts and negative LR

simPar <- list(s = "negativeLR_3trans", sample.size = 100, group = factor(c(rep("C1", 80), rep("C2", 20))),  pi.org = c(0.7, 0.2, 0.1) , g0.org = 5, nr.genes = 1e+03, nM = 1000, tot = "uni")




######### simulate...


mcCores <- 20


out.dir.s <- paste0(out.dir, "/", simPar$s, "/")
dir.create(out.dir.s, showWarnings=F, recursive=T)



sim <- simulate_from_DM(s = simPar$s, sample.size = simPar$sample.size, group = simPar$group, pi.org = simPar$pi.org, g0.org = simPar$g0.org, nr.genes = simPar$nr.genes, nM = simPar$nM, tot = simPar$tot, nD = simPar$nM, out.dir = out.dir.s , mc.cores=mcCores, save = FALSE)


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

dgeDMList[[mode]] <- dgeDM

save(dgeDMList, file = paste0(out.dir.s, "/dgeDMList.RData"))




cat("LR < 0 \n")
print(table(dgeDM$table$LR < 0))
head(dgeDM$table[dgeDM$table$LR < 0, ])




pdf(paste0(out.dir.s, "/hist_", mode,".pdf"))

hist(dgeDM$table$PValue, breaks = 100, col = "#1E90FF")
hist(dgeDM$table$LR, breaks = 100, col = "#1E90FF")

dev.off()










############ Compare constrOptim2 with constrOptim2G
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







############ Check the piH estimates for the genes where LR < 0

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)



plotProportions <- function(dgeDM, genes2plot, plotPath){
  
  pdf(plotPath, width = 10, height = 5)
  
  
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
    
    
    ### box plots with points - 2 groups only
    ggb <- ggplot(prop.smp.m, aes(x = ete_id, y = Proportions)) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 20, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), legend.position="none", plot.title = element_text(size=10)) +
      ggtitle(paste0(gene, "\n TagwiseDispersion = ", dgeDM$fit[[gene]]$gamma0, "\n LR = ", dgeDM$table[dgeDM$table$GeneID == gene, "LR"], " / PValue = ", dgeDM$table[dgeDM$table$GeneID == gene, "PValue"])) +    
      geom_jitter(aes(fill = Condition, colour = factor(Condition, labels=c("C1b", "C2b"))), position = position_jitterdodge(dodge.width = 0.75), alpha = 0.5) +
      geom_boxplot(aes(colour = Condition), fill = "white", outlier.size = NA, alpha = 0) + 
      geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, fill = Samples), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 2, shape = 19, colour = "black") +
      geom_point(data = prop.est.null, aes(x = ete_id, y = Proportions), size = 3, shape = 18, colour = "orange") +
      scale_colour_manual(values=c("C1"="firebrick", "C2"="dodgerblue4", "C1b"="firebrick1", "C2b" = "dodgerblue"))  +
      coord_cartesian(ylim = c(-0.1, 1.1)) 
    
    
    print(ggb)
    
    
  }
  
  
  dev.off()
  
  
  
}






genes2plot <- head(dgeDM$table[order(dgeDM$table$LR, decreasing = FALSE), "GeneID"])

plotPath <- paste0(out.dir.s, "/Proportions_", mode,"_negativeLR.pdf")

plotProportions(dgeDM, genes2plot, plotPath)







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

plotPath <- paste0(out.dir.s, "/Proportions_", mode,"_FP.pdf")

plotProportions(dgeDM, genes2plot, plotPath)



































