#######################################################
# for Committee meeting
# Created 16 July 2014 / 
# Last updated 16 July 2014 (02 Sep 2014)

#######################################################
# BioC 2.14


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2/")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata



#######################################################
# ->>>>>> load results
#######################################################



table <- read.table("PLOTS/Table_all_results.xls", header = T, stringsAsFactors = F)

colors <- c("hotpink", "magenta", "dodgerblue3",  "orange", "brown", "brown1")
names(colors) <- c("htseq_dexseq", "htseq_dexseq_simons", "fc_voomex",  "fc_DM", "htseq_DM", "miso_DM")                                  

name <- ""

dgec <- read.table("DM/fc/dge_counts_fc.xls", header = T, stringsAsFactors = F)

dgecnf <- read.table("DM/fc/dge_counts_fc_NOT_FILTERED.xls", header = T, stringsAsFactors = F)



#######################################################
# check FP in fc_DM
#######################################################


fpg <- as.character(na.omit(table[table$adjPValue_fc_DM < 0.05 & table$status == 0, "Gene"]))
length(fpg)

col1="blue"
col2="red"



pdf("PLOTS/ExploreDM/FalsePositives.pdf", width = 8, height = 6)


for(g in fpg){
  #g <- fpg[1]
  
  dge.g <- dgec[dgec$gene.id==g,,drop=FALSE]
  
  n <- length(dge.g$ete.id)

  par(mar=c(10, 4, 4, 2) +  0.1)
  plot(rep(1, n), type="n", ylim=c(0, 15), xaxt="n", ylab="log2(counts + 1)", xlab="")
  
  axis(side=1, at=1:n, labels=dge.g$ete.id, las=2)
  
  for(i in 1:3)
    lines(log2(dge.g[,2+i]+1), col=col1, type="l", lwd=2)
  for(i in 4:6)
    lines(log2(dge.g[,2+i]+1), col=col2, type="l", lwd=2)
  
  
}


dev.off()



#######################################################
# check uniq TP in fc_DM
#######################################################

library(limma)

# info about exons
simu_info.g <- read.table("Simu_info/true_genes_simulation.txt", header=T, stringsAsFactors = F)



# info about exons
simu_info.e <- read.table("Simu_info/true_exons_simulation.txt", header=T, stringsAsFactors = F)
simu_info.e <- simu_info.e[!is.na(simu_info.e$status_exon), ]

# differentially spliced exons
ds.ex <- simu_info.e[simu_info.e$status_exon==1, c("Gene", "exon")]


# ###

tp <- as.character(na.omit(table[table$adjPValue_fc_DM < 0.05 & 
                                    table$status == 1, "Gene"]))


tp1 <- as.character(na.omit(table[table$adjPValue_htseq_DM < 0.05 & 
                                    table$status == 1 , "Gene"]))

tp2 <- as.character(na.omit(table[table$adjPValue_htseq_dexseq_simons < 0.05 & 
                                    table$status == 1 , "Gene"]))


tp3 <- as.character(na.omit(table[table$adjPValue_fc_voomex < 0.05 & 
                                    table$status == 1 , "Gene"]))



uniq.tp <- setdiff(tp, c(tp1, tp2, tp3))
length(uniq.tp) # 196



col1="royalblue3"
col2="red3"

col3="orangered"


# "FBgn0000064"


pdf("PLOTS/ExploreDM/UniqTruePositives.pdf", width = 8, height = 6)
for(g in uniq.tp[1:20]){
  
  dge.g <- dgec[dgec$gene.id==g,,drop=FALSE]
  ds.ex.g <- ds.ex[ds.ex$Gene==g, , drop=FALSE]$exon
  
  col.axis <- ifelse(dge.g$ete.id %in% ds.ex.g, col3, "grey")
  
  n <- length(dge.g$ete.id)
  
  par(mar=c(10, 4, 4, 2) +  0.1)
  plot(rep(1, n), type="n", ylim=c(0, 15), xaxt="n", ylab="log2(counts + 1)", xlab="")
  
  for(j in 1:length(dge.g$ete.id))
    axis(side=1, at=j, labels=dge.g$ete.id[j], las=2, col.axis = col.axis[j] )
  
  for(i in 1:3)
    lines(log2(dge.g[,2+i]+1), col=col1, type="l", lwd=2)
  for(i in 4:6)
    lines(log2(dge.g[,2+i]+1), col=col2, type="l", lwd=2)
  
}
dev.off()



pdf("PLOTS/ExploreDM/UniqTruePositives2.pdf", width = 14, height = 7)

for(g in uniq.tp[1]){
  
  dge.g.keep <- dgec[dgec$gene.id==g,,drop=FALSE]
  dge.g <- dgecnf[dgecnf$gene.id==g,,drop=FALSE]
  ds.ex.g <- ds.ex[ds.ex$Gene==g, , drop=FALSE]$exon
  
  col.axis <- ifelse(dge.g$ete.id %in% ds.ex.g, col3, "grey50")
  
  
  labels <- strsplit2(dge.g$ete.id, ":")[,2]
  genes.keep <- dge.g$ete.id %in% dge.g.keep$ete.id
  labels[genes.keep] <- paste0("*" ,labels[genes.keep], "*")
  
  n <- length(dge.g$ete.id)
  
  par(mar=c(8, 5, 4, 2) +  0.1)
  plot(rep(1, n), type="n", ylim=c(0, 10), xaxt="n", ylab="",xlab="", main=g,  cex.main=2, cex.axis=2)
  title(ylab = "log2(counts + 1)", mgp = c(3, 1, 0), cex.lab=1.5)
  title(xlab = "Flattened exon ID", mgp = c(6, 1, 0), cex.lab=1.5)
  
  for(j in 1:length(dge.g$ete.id))
    axis(side=1, at=j, labels=labels[j], las=2, col.axis = col.axis[j], cex.axis=2)
  
  for(i in 1:3)
    lines(log2(dge.g[,2+i]+1), col=col1, type="l", lwd=3)
  for(i in 4:6)
    lines(log2(dge.g[,2+i]+1), col=col2, type="l", lwd=3)
   
  legend("topright", c("Condition1", "Condition2"), lty=1, lwd=3, col=c(col1, col2), cex=2)
  
}
dev.off()



#######################################################
# check FN in fc_DM
#######################################################

# info about exons
simu_info.e <- read.table("Simu_info/true_exons_simulation.txt", header=T, stringsAsFactors = F)
simu_info.e <- simu_info.e[!is.na(simu_info.e$status_exon), ]

ds.ex <- simu_info.e[simu_info.e$status_exon==1, c("Gene", "exon")]

table(table(ds.ex$Gene))




###

fn_fc_DM <- as.character(na.omit(table[table$adjPValue_fc_DM > 0.05 & 
                                   table$status == 1, "Gene"]))

tp1 <- as.character(na.omit(table[table$adjPValue_htseq_DM < 0.05 & 
                                    table$status == 1 , "Gene"]))


tp2 <- as.character(na.omit(table[table$adjPValue_htseq_dexseq_simons < 0.05 & 
                                    table$status == 1 , "Gene"]))



fn <- intersect(fn_fc_DM, tp2)
fn <- setdiff(fn, tp1)



pdf("PLOTS/ExploreDM/FN.pdf", width = 8, height = 6)


for(g in fn){
  # g=fn[1]
  
  dge.g <- dgec[dgec$gene.id==g,,drop=FALSE]
  ds.ex.g <- ds.ex[ds.ex$Gene==g, , drop=FALSE]$exon
  
  col.axis <- ifelse(dge.g$ete.id %in% ds.ex.g, "magenta", "grey")
  
  n <- length(dge.g$ete.id)
  
  par(mar=c(10, 4, 4, 2) +  0.1)
  plot(rep(1, n), type="n", ylim=c(0, 15), xaxt="n", ylab="log2(counts + 1)", xlab="")
  
  for(j in 1:length(dge.g$ete.id))
  axis(side=1, at=j, labels=dge.g$ete.id[j], las=2, col.axis = col.axis[j] )
  
  for(i in 1:3)
    lines(log2(dge.g[,2+i]+1), col=col1, type="l", lwd=2)
  for(i in 4:6)
    lines(log2(dge.g[,2+i]+1), col=col2, type="l", lwd=2)
  
  
}


dev.off()



###


fn_fc_DM <- as.character(na.omit(table[table$adjPValue_fc_DM >= 0.05 & 
                                         table$status == 1, "Gene"]))

length(fn_fc_DM)

tp_fc_DM <- as.character(na.omit(table[table$adjPValue_fc_DM < 0.05 & 
                                         table$status == 1, "Gene"]))

length(tp_fc_DM)



DSexons <- rep(0, length(fn_fc_DM))


for(j in 1:length(fn_fc_DM)){
  
  g <- fn_fc_DM[j]
  
  dge.g <- dgec[dgec$gene.id==g,,drop=FALSE]
  ds.ex.g <- ds.ex[ds.ex$Gene==g, , drop=FALSE]$exon
  
  DSexons[j] <- sum(dge.g$ete.id %in% ds.ex.g)
  
  
}

table(DSexons)




#######################################################
# check TP in fc_DM after filtering
#######################################################



tp <- as.character(na.omit(table[table$adjPValue_fc_DM < 0.05 & 
                                   table$status == 1, "Gene"]))


length(tp)



DSexons <- rep(0, length(tp))
TESTexons <- rep(0, length(tp))


for(j in 1:length(tp)){
  # j=1
  g <- tp[j]
  
  dge.g <- dgec[dgec$gene.id==g,,drop=FALSE]
  ds.ex.g <- ds.ex[ds.ex$Gene==g, , drop=FALSE]$exon

  DSexons[j] <- sum(dge.g$ete.id %in% ds.ex.g)
  TESTexons[j] <- dim(dge.g)[1]
  
}

table(DSexons)

table(DSexons, TESTexons)





##############################################################################################################
# trend for DM: overall gene expression vs dispersion
##############################################################################################################

source("/home/gosia/R/R_Multinomial_project/Analysis_SimDroV2/Analysis_for_Committee_meeting1/dirmult.R")
source("/home/gosia/R/R_Multinomial_project/Analysis_SimDroV2/Analysis_for_Committee_meeting1/glmFitDM.R")



fpg <- as.character(na.omit(table[table$adjPValue_fc_DM < 0.05 & table$status == 0, "Gene"]))
length(fpg)


name1 <- "fc"
load(paste0("DM/fc/DM_",name1,"_results.RData"))

name1 <- "htseq"
load(paste0("DM/htseq/DM_",name1,"_results.RData"))


glmlrt <- g2lrt

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
table.gh0 <- gh0



pdf(paste0("PLOTS_back2_PhD_CM/ExploreDM/", name1, "_TREND_meanY01_vs_gh0.pdf"), width = 21, height = 7)
# par(mfrow=c(1,3), cex=1.5)
par(mfrow=c(1,3), cex.main=1.2, cex.axis = 1.5, cex.lab=1.5, cex=1.5)

smoothScatter(log(table.gh0$meanY0), log(table.gh0$gh0.0f), col=2, main="Y1", xlab="log mean gene expression", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white", blues9)) )
smoothScatter(log(table.gh0$meanY1) , log(table.gh0$gh0.1f), col=3, main="Y2", xlab="log mean gene expression", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white", blues9)))
smoothScatter(log(table.gh0$meanY) , log(table.gh0$gh0.0n), col=4, main="Y12", xlab="log mean gene expression", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white",blues9)))

dev.off()



table.gh0 <- merge(gh0, table, by=1, all.y=TRUE)

table.gh0 <- table.gh0[!is.na(table.gh0$adjPValue_fc_DM), ]

table.gh0 <- table.gh0[table.gh0$adjPValue_fc_DM < 0.05 & table.gh0$status ==1, ]


pdf(paste0("PLOTS/ExploreDM/", name, "TREND_meanY01_vs_gh0_TP.pdf"), width = 21, height = 7)
par(mfrow=c(1,3), cex.main=1.2, cex.axis = 1.5, cex.lab=1.5, cex=1.5)

smoothScatter(log(table.gh0$meanY0), log(table.gh0$gh0.0f), col=2, main="Y1", xlab="log mean gene expression", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white", blues9)) )
smoothScatter(log(table.gh0$meanY1) , log(table.gh0$gh0.1f), col=3, main="Y2", xlab="log mean gene expression", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white", blues9)))
smoothScatter(log(table.gh0$meanY) , log(table.gh0$gh0.0n), col=4, main="Y12", xlab="log mean gene expression", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white",blues9)))

dev.off()



#######################################################
# trend for DM: dispersion vs. entropy
#######################################################

source("/home/gosia/R/R_Multinomial_project/Analysis_SimV2/Analysis_for_Committee_meeting1/dirmult.R")
source("/home/gosia/R/R_Multinomial_project/Analysis_SimV2/Analysis_for_Committee_meeting1/glmFitDM.R")



name1 <- "fc"
load(paste0("DM/fc/DM_",name1,"_results.RData"))
glmlrt <- g2lrt




dDmultinom <- function (x, gamma, log = FALSE){  
  
  size <- sum(x)
  gammaP <- sum(gamma)
  
  r <- lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(gammaP +1) - lgamma(size + gammaP +1) + sum(lgamma(x + gamma +1) - lgamma(gamma)) 
  
  if (log)
    r
  else exp(r)
}


entropy <- function(Y, gamma){  ##### something must be wrong with this calculations because Entropy H > 0 
  
  lp <- rep(0, nrow(Y))
  
  for(i in 1:nrow(Y)){
    
    lp[i] <- dDmultinom(x=Y[i,], gamma=gamma, log = TRUE)
    
  }
  
  H <- - sum(exp(lp) * lp)
  
  return(H)
  
}



gh0.tmp <-  mclapply(names(glmlrt$fit.full), function(g){
  # g="FBgn0000064"
  
  if(!is.null(glmlrt$fit.full[[g]])){           
    
    
    meanY <- entropy(glmlrt$fit.null[[g]]$Y, glmlrt$fit.null[[g]]$gh[1,])
    meanY0 <- entropy(glmlrt$fit.full[[g]]$Y[1:3,], glmlrt$fit.full[[g]]$gh[1,])
    meanY1 <- entropy(glmlrt$fit.full[[g]]$Y[4:6,], glmlrt$fit.full[[g]]$gh[4,])
    
    
    gh0.0f <- glmlrt$fit.full[[g]]$gh0[1]
    gh0.1f <- glmlrt$fit.full[[g]]$gh0[4]
    gh0.0n <- glmlrt$fit.null[[g]]$gh0[1]
    
    return(data.frame(gh0.0f=gh0.0f, gh0.1f=gh0.1f, gh0.0n=gh0.0n, meanY=meanY, meanY0=meanY0, meanY1=meanY1))
  }    
  return(data.frame(gh0.0f=NA, gh0.1f=NA, gh0.0n=NA, meanY=NA, meanY0=NA, meanY1=NA))
  
}, mc.cores=10)


gh0 <- do.call(rbind, gh0.tmp)
gh0 <- data.frame(GeneID=names(glmlrt$fit.full), gh0)

head(gh0)
table.gh0 <- gh0



pdf(paste0("PLOTS/ExploreDM/", name, "TREND_hg0_vs_entropy.pdf"), width = 21, height = 7)
# par(mfrow=c(1,3), cex=1.5)
par(mfrow=c(1,3), cex.main=1.2, cex.axis = 1.5, cex.lab=1.5, cex=1.5)

smoothScatter(log(-table.gh0$meanY0), log(table.gh0$gh0.0f), col=2, main="Y1", xlab="log (-Entropy)", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), xlim=c(0, 200), colramp=colorRampPalette(c("white", blues9)) )
smoothScatter(log(-table.gh0$meanY1) , log(table.gh0$gh0.1f), col=3, main="Y2", xlab="log (-Entropy)", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), xlim=c(0, 200),colramp=colorRampPalette(c("white", blues9)))
smoothScatter(log(-table.gh0$meanY) , log(table.gh0$gh0.0n), col=4, main="Y12", xlab="log (-Entropy)", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), xlim=c(0, 200), colramp=colorRampPalette(c("white",blues9)))

dev.off()





pdf(paste0("PLOTS/ExploreDM/", name, "TREND_hg0_vs_entropy2.pdf"), width = 21, height = 7)
# par(mfrow=c(1,3), cex=1.5)
par(mfrow=c(1,3), cex.main=1.2, cex.axis = 1.5, cex.lab=1.5, cex=1.5)

smoothScatter(-table.gh0$meanY0, log(table.gh0$gh0.0f), col=2, main="Y1", xlab="Entropy", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), xlim=c(0,100) ,colramp=colorRampPalette(c("white", blues9)) )
smoothScatter(-table.gh0$meanY1 , log(table.gh0$gh0.1f), col=3, main="Y2", xlab="Entropy", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white", blues9)))
smoothScatter(-table.gh0$meanY , log(table.gh0$gh0.0n), col=4, main="Y12", xlab="Entropy", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white",blues9)))

dev.off()








table.gh0 <- merge(gh0, table, by=1, all.y=TRUE)

table.gh0 <- table.gh0[!is.na(table.gh0$adjPValue_fc_DM), ]

table.gh0 <- table.gh0[table.gh0$adjPValue_fc_DM < 0.05 & table.gh0$status ==1, ]


pdf(paste0("PLOTS/ExploreDM/", name, "TREND_meanY01_vs_gh0_TP.pdf"), width = 21, height = 7)
par(mfrow=c(1,3), cex.main=1.2, cex.axis = 1.5, cex.lab=1.5, cex=1.5)

smoothScatter(log(table.gh0$meanY0), log(table.gh0$gh0.0f), col=2, main="Y1", xlab="log mean gene expression", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white", blues9)) )
smoothScatter(log(table.gh0$meanY1) , log(table.gh0$gh0.1f), col=3, main="Y2", xlab="log mean gene expression", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white", blues9)))
smoothScatter(log(table.gh0$meanY) , log(table.gh0$gh0.0n), col=4, main="Y12", xlab="log mean gene expression", ylab="log gamma +", nrpoints = Inf, ylim=c(0, 25), colramp=colorRampPalette(c("white",blues9)))

dev.off()



##############################################################################################################
# check DM estimates 
##############################################################################################################


source("/home/gosia/R/R_Multinomial_project/Analysis_SimV2/Analysis_for_Committee_meeting1/dirmult.R")
source("/home/gosia/R/R_Multinomial_project/Analysis_SimV2/Analysis_for_Committee_meeting1/glmFitDM.R")

name1 <- "fc"
load(paste0("DM/fc/DM_",name1,"_results.RData"))

fitf <- g2lrt$fit.full
fitn <- g2lrt$fit.null
tbl <- g2lrt$table

# for FP

fpg <- as.character(na.omit(table[table$adjPValue_fc_DM < 0.05 & table$status == 0, "Gene"]))
length(fpg)

tbl[tbl$GeneID %in% fpg, ]

fitf.fp <- fitf[fpg]
fitn.fp <- fitn[fpg]

fitf.fp[["FBgn0039101"]]
fitn.fp[["FBgn0039101"]]





# for FN

fitf[["FBgn0086657"]]
fitn[["FBgn0086657"]]

table.gh0[table.gh0$GeneID=="FBgn0086657" , ]













