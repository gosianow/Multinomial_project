#######################################################
# Created 22 Jul 2014 
# run on simulated data: Simulations 12, 14 from Katarina
# ROC plots
# histograms
# run new edgeR and voomex 
# compare with Katarina's results voomex, edgeR, DEXSeq
#######################################################
# BioC 2.14

setwd("/home/gosia/Multinomial_project/Simulations")

# load information about simulation 
simu_info <- read.table("Simu_info/simu14_reads_final.txt", header=T)
head(simu_info)


# create metadata file
metadata <- data.frame(SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "R",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)



# ### BioC 2.13
# library("DEXSeq")
# # run read.HTSeqCounts to create ExonCountSet object
# 
# ecs <- read.HTSeqCounts(countfiles = paste0("MISO_counts/", metadata$SampleName, "_miso_sim14.txt"), design = metadata[,c("condition")])
# 
# # ecs <- read.HTSeqCounts(countfiles = paste0("DEXSeq_counts/sim14", metadata$SampleName, "rno_snfb.txt"), design = metadata[,c("condition")])
# 
# sampleNames(ecs) = metadata$SampleName
# 
# counts <- counts(ecs)
# gene.id <- fData(ecs)$geneID
# ete.id <- fData(ecs)$exonID


library(DEXSeq)

ecs <- DEXSeqDataSetFromHTSeq(countfiles = paste0("MISO_counts/", metadata$SampleName, "_miso_sim14.txt"), sampleData= data.frame(row.names = metadata$SampleName, condition = metadata$condition), design = ~sample + exon + condition:exon)

counts <- featureCounts(ecs)
gene.id <- unlist(lapply(seq(rownames(counts)), function(i){e <- rownames(counts)[i];strsplit(e, ":")[[1]][1] }))
ete.id <- unlist(lapply(seq(rownames(counts)), function(i){e <- rownames(counts)[i];strsplit(e, ":")[[1]][2] }))



library("edgeR")
# create DGEList object
dge <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))
# filter counts by log-cpm (for 4 vs 3 comparison)
keep <- rowSums(cpm(dge)>1) > 3
dge <- dge[keep,]
# normalisation
dge <- calcNormFactors(dge)


# # for MISO: filter events with only 0s: 0,0,0,0,0
# ete.id.tmp <- gsub("E", "", dge$genes$ete.id)
# ete.id.tmp <- strsplit(ete.id.tmp, ",")
# keep2 <- unlist( lapply(ete.id.tmp, function(ev){ sum(as.numeric(ev)) > 0 }) )
# dge <- dge[keep2,]


# # for MISO: filter events with only 1s: 1,1,1,1
# ete.id.tmp <- gsub("E", "", dge$genes$ete.id)
# ete.id.tmp <- strsplit(ete.id.tmp, ",")
# keep3 <- unlist( lapply(ete.id.tmp, function(ev){ !sum(as.numeric(ev)) == length(ev) }) )
# dge <- dge[keep3,]


## add 1 to counts
#dge$counts <- dge$counts + 1


## save dge
dir.create("PLOTS", showWarnings=F, recursive=T)
#name1 <- "MISO_p1"
#name1 <- "MISO_filt0s"
name1 <- "MISO"

write.table(dge$counts, paste0("PLOTS/dge_counts_",name1,".xls"), quote=F, sep="\t", row.names=T, col.names=T)

design <- model.matrix(~condition, data=metadata)
design



source("/home/gosia/R/R_Multinomial_project/dirmult.R")
source("/home/gosia/R/R_Multinomial_project/glmFitDM.R")



#######################################################
# run GLM 
#######################################################

# name <- "optim_p1_"
# 
# 
# glmfit <- glmFitDM(dge, design, optimization="optim")
# 
# glmlrt <- glmLRTDM(glmfit, coef=2)
# 
# 
# head(glmlrt$table)
# 
# # save(glmlrt, file=paste0("PLOTS/", name, "glmlrt.RData"))


#######################################################
# run G2
#######################################################


name <- "MISO_G2_EST_"


g2fit <- g2FitDM(dge)

g2lrt <- g2LRTDM(g2fit)



head(g2lrt$table)

glmlrt <- g2lrt

# save(g2lrt, file=paste0("PLOTS/", name, "g2lrt.RData"))


#######################################################
# run voomex
#######################################################
library(limma)

out.dir <- "Results_voomex"
dir.create(out.dir)

pdf(paste0(out.dir, "/MDS.pdf"))
plotMDS(dge)
dev.off()


pdf(paste0(out.dir, "/voom.pdf"))
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
plotSA(fit)
dev.off()

vex <- diffSplice(fit, geneid=fit$genes$gene.id)

# top table
top.exon <- topSplice(vex, coef=2, level="exon", n=Inf)
top.gene <- topSplice(vex, coef=2, level="gene", n=Inf) 


write.table(top.exon, paste0(out.dir, "/MISO_voom_sim14_exon.xls"), sep="\t", row.names=FALSE, quote=F) 
write.table(top.gene, paste0(out.dir, "/MISO_voom_sim14_gene.xls"), sep="\t", row.names=FALSE, quote=F) 


#######################################################
# run edgeR
#######################################################

out.dir <- "Results_edgeR"
dir.create(out.dir)

sv <- spliceVariants(dge, geneID=dge$genes$gene.id)

edgeR_res <- topTags(sv, n=nrow(sv$table), adjust.method="BH", sort.by="PValue")

write.table(edgeR_res, paste0(out.dir, "/MISO_edgeR_sim14.xls"), sep="\t", row.names=FALSE, quote=F) 



#######################################################
# merge results with simu_info
#######################################################


table.tmp <- glmlrt$table
dim(table.tmp)

simu_info.tmp <- simu_info[!duplicated(simu_info$Gene),c("Gene", "status", "n", "rpkm")]
dim(simu_info.tmp)
head(simu_info.tmp)


table <- merge(table.tmp, simu_info.tmp, by=1, all.x=TRUE, sort=FALSE)

## remove genes with NA
table[!complete.cases(table), ]

table <- table[complete.cases(table), ]
dim(table)


# write.table(table, paste0("PLOTS/",name,"table.xls"), quote=F, sep="\t", row.names=F, col.names=T)


glmlrt$counts[glmlrt$genes$gene.id == "FBgn0000052", ]


# ### add other results produced by K.
# 
# ## MISO
# out_edgeR.org <- read.table(paste0("Results/MISO_edgeR_sim14.xls"), header=T)
# out_voom.org <- read.table("Results/MISO_voom_sim14_gene.xls", header=T)
# 
# # ## DEXSeq counts
# # out_edgeR.org <- read.table(paste0("Results/DEXSeq_edgeR_sim14.xls"), header=T)
# # out_voom.org <- read.table("Results/DEXSeq_voom_sim14_gene.xls", header=T)
# 
# 
# out_edgeR <- out_edgeR.org[,c("GeneID", "PValue", "FDR")]
# colnames(out_edgeR) <- c("GeneID", "PValue_edgeR", "adjPValue_edgeR")
# 
# out_voom <- out_voom.org[,c("GeneID", "P.Value", "FDR")]
# colnames(out_voom) <-  c("GeneID", "PValue_voom", "adjPValue_voom")
# 
# table.all <- merge(table, out_edgeR, by=1, all.x=TRUE, sort=FALSE)
# table.all <- merge(table.all, out_voom, by=1, all.x=TRUE, sort=FALSE)
# 
# table.all <- table.all[complete.cases(table.all), ]
# 
# write.table(table.all, paste0("PLOTS/",name,"table_all.xls"), row.names=F, col.names=T, quote=F, sep="\t")



### Results voomex & edgeR ESTimated

out_voom.org <- read.table("Results_voomex/MISO_voom_sim14_gene.xls", header=T)
out_edgeR.org <- read.table(paste0("Results_edgeR/MISO_edgeR_sim14.xls"), header=T)

out_edgeR <- out_edgeR.org[,c("GeneID", "PValue", "FDR")]
colnames(out_edgeR) <- c("GeneID", "PValue_edgeR", "adjPValue_edgeR")

out_voom <- out_voom.org[,c("GeneID", "P.Value", "FDR")]
colnames(out_voom) <-  c("GeneID", "PValue_voom", "adjPValue_voom")

table.all <- merge(table, out_edgeR, by=1, all.x=TRUE, sort=FALSE)
table.all <- merge(table.all, out_voom, by=1, all.x=TRUE, sort=FALSE)

table.all <- table.all[complete.cases(table.all), ]

write.table(table.all, paste0("PLOTS/",name,"table_all.xls"), row.names=F, col.names=T, quote=F, sep="\t")


### add DEXSeq results from Katarina

out_dexseq.org <- read.table("Results/DEXSeq_MISO_sim14_genelevel.xls", header=T)

out_dexseq <- out_dexseq.org[,c("Gene", "padjust")]
colnames(out_dexseq) <- c("GeneID",  "adjPValue_dexseq")

table.all <- merge(table.all, out_dexseq, by=1, all.x=TRUE, sort=FALSE)

table.all <- table.all[complete.cases(table.all), ]

write.table(table.all, paste0("PLOTS/",name,"table_all.xls"), row.names=F, col.names=T, quote=F, sep="\t")




#######################################################
# histograms of p-values
#######################################################


pdf(paste0("PLOTS/", name, "hist_pvalues.pdf"), w=15, h=5)
par(mfrow=c(1,3))
hist(table.all$PValue, col="orange", breaks=50, main="DM", cex.main=3, cex.lab=1.5, cex.axis = 1.5)
hist(table.all$PValue_edgeR, col="green", breaks=50, main="edgeR", cex.main=3, cex.lab=1.5, cex.axis = 1.5)
hist(table.all$PValue_voom, col="grey", breaks=50, main="voomex", cex.main=3, cex.lab=1.5, cex.axis = 1.5)
dev.off()


#######################################################
# ROC plot
#######################################################


# install.packages("/home/gosia/R/packages/benchTools_0.0-2.tar.gz", dependencies=TRUE, lib="/home/gosia/R/libraries/3.1.0")

library(benchTools)


table$score <- 1-table$PValue

score <- table$score
label <- table$status
score.X <- 1-table$FDR


# pdf(paste0("PLOTS/",name,"roc.pdf"))
# r <- rocX(score, label)
# plot(r, col="orange")
# r <- rocX(score, label, score.X = score.X, threshold.X = 0.95)
# plot(r, col="orange")
# dev.off()


#### ROC with other results


score <- as.matrix(1-table.all[,c("PValue","PValue_edgeR", "PValue_voom")])
label <- table.all$status

r <- rocX(score, label)

pdf(paste0("PLOTS/",name,"roc_all.pdf"))
plot(r, col=c("orange", "green", "blue"), cex.main=3, cex.lab=1.5, cex.axis = 2, lwd=3)
legend("bottomright", c("DM", "edgeR", "voom"), col=c("orange", "green", "blue"), lty=1, lwd=2.5)
dev.off()



#######################################################
# ROC plot - stratify by df
#######################################################

table(table$df)

p.df <- c( 1, 2, 3, 4, 5 )


pdf(paste0("PLOTS/",name,"rocSTR.pdf"))
r <- rocX(score, label)
plot(r, col="orange", lwd=2, lty=2)
leg <- "all"
lty <- 2
col <- "orange"

indx <- table$df <= p.df[1]
r <- rocX(score[indx], label[indx])
lines(r[[1]]@x.values[[1]], r[[1]]@y.values[[1]], col=1, lty=1, lwd=1)
leg <- c( leg, paste0("df <= ", p.df[1]) )
lty <- c(lty, 1)
col <- c(col, 1)

for( i in 1:(length(p.df)-1) ){  
  
  indx <- table$df > p.df[i] & table$df <= p.df[i+1]
  r <- rocX(score[indx], label[indx])
  lines(r[[1]]@x.values[[1]], r[[1]]@y.values[[1]], col=i+1, lty=1, lwd=1) 
  leg <- c(leg, paste0("df > ", p.df[i], " & df <= ", p.df[i+1]))
  lty <- c(lty, 1)
  col <- c(col, i+1)
  
}

indx <- table$df > p.df[length(p.df)]
r <- rocX(score[indx], label[indx])
lines(r[[1]]@x.values[[1]], r[[1]]@y.values[[1]], col=length(p.df)+1, lty=1, lwd=1)
leg <- c(leg, paste0("df > ", p.df[length(p.df)]))
lty <- c(lty, 1)
col <- c(col, length(p.df)+1)

legend("bottomright", leg, lty=lty, col=col, lwd=2)

hist(table$df, col="darkcyan", breaks=range(table$df)[2]+1)

dev.off()



###

score <- 1-table.all[,c("PValue_edgeR")]
label <- table.all$status


p.df <- c( 1, 2, 3, 4, 5 )


pdf(paste0("PLOTS/",name,"rocSTR_edgeR.pdf"))
r <- rocX(score, label)
plot(r, col="green", lwd=2, lty=2, main="edgeR")
leg <- "all"
lty <- 2
col <- "green"

indx <- table$df <= p.df[1]
r <- rocX(score[indx], label[indx])
lines(r[[1]]@x.values[[1]], r[[1]]@y.values[[1]], col=1, lty=1, lwd=1)
leg <- c( leg, paste0("df <= ", p.df[1]) )
lty <- c(lty, 1)
col <- c(col, 1)

for( i in 1:(length(p.df)-1) ){  
  
  indx <- table$df > p.df[i] & table$df <= p.df[i+1]
  r <- rocX(score[indx], label[indx])
  lines(r[[1]]@x.values[[1]], r[[1]]@y.values[[1]], col=i+1, lty=1, lwd=1) 
  leg <- c(leg, paste0("df > ", p.df[i], " & df <= ", p.df[i+1]))
  lty <- c(lty, 1)
  col <- c(col, i+1)
  
}

indx <- table$df > p.df[length(p.df)]
r <- rocX(score[indx], label[indx])
lines(r[[1]]@x.values[[1]], r[[1]]@y.values[[1]], col=length(p.df)+1, lty=1, lwd=1)
leg <- c(leg, paste0("df > ", p.df[length(p.df)]))
lty <- c(lty, 1)
col <- c(col, length(p.df)+1)

legend("bottomright", leg, lty=lty, col=col, lwd=2)

dev.off()







score <- 1-table.all[,c("PValue_voom")]
label <- table.all$status


p.df <- c( 1, 2, 3, 4, 5 )


pdf(paste0("PLOTS/",name,"rocSTR_voom.pdf"))
r <- rocX(score, label)
plot(r, col="grey", lwd=2, lty=2, main="voom")
leg <- "all"
lty <- 2
col <- "grey"

indx <- table$df <= p.df[1]
r <- rocX(score[indx], label[indx])
lines(r[[1]]@x.values[[1]], r[[1]]@y.values[[1]], col=1, lty=1, lwd=1)
leg <- c( leg, paste0("df <= ", p.df[1]) )
lty <- c(lty, 1)
col <- c(col, 1)

for( i in 1:(length(p.df)-1) ){  
  
  indx <- table$df > p.df[i] & table$df <= p.df[i+1]
  r <- rocX(score[indx], label[indx])
  lines(r[[1]]@x.values[[1]], r[[1]]@y.values[[1]], col=i+1, lty=1, lwd=1) 
  leg <- c(leg, paste0("df > ", p.df[i], " & df <= ", p.df[i+1]))
  lty <- c(lty, 1)
  col <- c(col, i+1)
  
}

indx <- table$df > p.df[length(p.df)]
r <- rocX(score[indx], label[indx])
lines(r[[1]]@x.values[[1]], r[[1]]@y.values[[1]], col=length(p.df)+1, lty=1, lwd=1)
leg <- c(leg, paste0("df > ", p.df[length(p.df)]))
lty <- c(lty, 1)
col <- c(col, length(p.df)+1)

legend("bottomright", leg, lty=lty, col=col, lwd=2)

dev.off()




#######################################################
# TPR vs FPR
#######################################################


label <- table$status
adj.p.value <- table$FDR


TPRvsFDR <- function(label, adj.p.value, FDR.cut.off=seq(from=0.01, to=1, by=0.01)){
  
  n <- length(label)
  q <- length(FDR.cut.off)
  TPR <- rep(0, q)
  FDR <- rep(0, q)
  
  for(i in 1:q){
    # i=1
    label.est <- as.numeric(adj.p.value <= FDR.cut.off[i])
    
    TP <- sum(label==1 & label.est==1)
    FP <- sum(label==0 & label.est==1)
    FN <- sum(label==1 & label.est==0)
    
    TPR[i] <- TP/(TP+FN)
    FDR[i] <- FP/(FP+TP)
      
  }  
  return(cbind( FDR , TPR))  
}





pdf(paste0("PLOTS/", name, "TPRvsachievedFDR.pdf"))

FDR.cut.off=c(0.01, 0.05, 0.1)

tf <- TPRvsFDR(label, adj.p.value, FDR.cut.off=FDR.cut.off)
plot(tf, xlim=c(0,0.3), ylim=c(0,1), pch=19, col="orange", cex.main=1.5, cex.lab=1.5, cex.axis = 1.5, cex=2 )
abline(v=FDR.cut.off, lty=3)

tf <- TPRvsFDR(table.all$status, table.all$adjPValue_edgeR, FDR.cut.off=FDR.cut.off)
points(tf, col="green", cex=2, pch=19)

tf <- TPRvsFDR(table.all$status, table.all$adjPValue_voom, FDR.cut.off=FDR.cut.off)
points(tf, col="blue", cex=2, pch=19)

tf <- TPRvsFDR(table.all$status, table.all$adjPValue_dexseq, FDR.cut.off=FDR.cut.off)
points(tf, col="magenta", cex=2, pch=19)


legend("bottomright", c("DM", "edgeR", "voom", "DEXSeq"), col=c("orange", "green", "blue", "magenta"), pch=19)


dev.off()




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





pdf(paste0("PLOTS/", name, "TREND_meanY01_vs_gh0.pdf"), w=30, h=10)
par(mfrow=c(1,3), cex.main=1.5, cex.axis = 1.3, cex=2 )

plot(table.gh0$meanY0, table.gh0$gh0.0f, col=2, log="xy", main="G1", xlab="", ylab="")
plot(table.gh0$meanY1 , table.gh0$gh0.1f, col=3, log="xy", main="G2", xlab="", ylab="")
plot(table.gh0$meanY , table.gh0$gh0.0n, col=4, log="xy", main="null: G1 + G2", xlab="", ylab="")

dev.off()


pdf(paste0("PLOTS/", name, "TRENDg1_meanY01_vs_gh0.pdf"))

plot(table.gh0$meanY0, table.gh0$gh0.0f, col=2, log="xy", xlab="mean gene expression", ylab=expression(paste(gamma, "+")), cex.axis = 1.5, cex.lab=1.5, cex=2)

dev.off()



#######################################################
# venn diagrams
#######################################################

table.all <- read.table("PLOTS/MISO_G2_EST_table_all.xls", header=TRUE)

FDR.cutoff <- 0.05


library(VennDiagram)



colors.results <- c("orange", "green", "blue", "magenta", "grey")
names(colors.results) <- c("DM", "edgeR", "voomex", "DEXSeq", "True")


plot.venn <- function(venne.genes, colors.results, venn.methods, name=""){
  
  colors.col <- colors.results[venn.methods]
  
  if(length(venn.methods)==4)
    colors.col <- colors.col[c(1, 3, 4 ,2)]
  
  venn.d <- venn.diagram(venne.genes[venn.methods], filename=NULL, fill = colors.results[venn.methods], col=colors.col, cex=1, cat.cex=1.4, lwd=2, lty=1, alpha=0.5, margin=0.05)
  
  pdf(paste0("PLOTS/", name, "Venn_Diagr.pdf"))
  grid.draw(venn.d)
  dev.off()
  
  
}

## all discoveries

venne.genes <- list()

venne.genes$DM <- table.all[table.all$FDR < FDR.cutoff, "GeneID"]
venne.genes$edgeR <- table.all[table.all$adjPValue_edgeR < FDR.cutoff, "GeneID"]
venne.genes$voomex <- table.all[table.all$adjPValue_voom < FDR.cutoff, "GeneID"]


plot.venn(venne.genes, colors.results, venn.methods = c("DM", "edgeR", "voomex"), name=name)


## TP discoveries

venne.genes <- list()

venne.genes$DM <- table.all[table.all$status == 1 & table.all$FDR < FDR.cutoff, "GeneID"]
venne.genes$edgeR <- table.all[table.all$status == 1 &table.all$adjPValue_edgeR < FDR.cutoff, "GeneID"]
venne.genes$voomex <- table.all[table.all$status == 1 &table.all$adjPValue_voom < FDR.cutoff, "GeneID"]
venne.genes$DEXSeq <- table.all[table.all$status == 1 &table.all$adjPValue_dexseq < FDR.cutoff, "GeneID"]

plot.venn(venne.genes, colors.results, venn.methods = c("DM", "edgeR", "voomex", "DEXSeq"), name=paste0(name, "TP_"))



## FP discoveries

venne.genes <- list()

venne.genes$DM <- table.all[table.all$status == 0 & table.all$FDR < FDR.cutoff, "GeneID"]
venne.genes$edgeR <- table.all[table.all$status == 0 &table.all$adjPValue_edgeR < FDR.cutoff, "GeneID"]
venne.genes$voomex <- table.all[table.all$status == 0 &table.all$adjPValue_voom < FDR.cutoff, "GeneID"]


plot.venn(venne.genes, colors.results, venn.methods = c("DM", "edgeR", "voomex"), name=paste0(name, "FP_"))



## 

venne.genes <- list()

venne.genes$DM <- table.all[ table.all$FDR < FDR.cutoff, "GeneID"]
venne.genes$edgeR <- table.all[table.all$adjPValue_edgeR < FDR.cutoff, "GeneID"]
venne.genes$voomex <- table.all[table.all$adjPValue_voom < FDR.cutoff, "GeneID"]
venne.genes$DEXSeq <- table.all[table.all$adjPValue_dexseq < FDR.cutoff, "GeneID"]
venne.genes$True <- table.all[table.all$status == 1, "GeneID"]


plot.venn(venne.genes, colors.results, venn.methods = c("DM", "edgeR", "voomex", "DEXSeq", "True"), name=paste0(name, "All_"))















































