#######################################################
# run on REAL data
#######################################################


setwd("/home/gosia/Multinomial_project/kaessmann_hsa_tissue")


# create metadata file
metadata <- read.table("metadata.xls", header=T, sep="\t")

metadata <- metadata[c(2, 4, 6, 7:9),]
metadata

library("edgeR")
library("DEXSeq")

# run read.HTSeqCounts to create ExonCountSet object
metadata$condition <- as.factor(as.character(metadata$TissueType))


ecs <- read.HTSeqCounts(countfiles = paste0("MISO_counts/MISO_counts_", metadata$SampleName, ".txt"), design = metadata[,c("condition")])


sampleNames(ecs) = metadata$SampleName


counts <- counts(ecs)
gene.id <- fData(ecs)$geneID
ete.id <- fData(ecs)$exonID

# create DGEList object
dge <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))
# normalisation
dge <- calcNormFactors(dge)
# filter counts by log-cpm (for 4 vs 3 comparison)
keep <- rowSums(cpm(dge)>1)>=5
dge <- dge[keep,]

# for MISO: filter events with only 0s: 0,0,0,0,0
ete.id.tmp <- gsub("E", "", dge$genes$ete.id)
ete.id.tmp <- strsplit(ete.id.tmp, ",")
keep2 <- unlist( lapply(ete.id.tmp, function(ev){ sum(as.numeric(ev)) > 0 }) )
dge <- dge[keep2,]


# # for MISO: filter events with only 1s: 1,1,1,1
# ete.id.tmp <- gsub("E", "", dge$genes$ete.id)
# ete.id.tmp <- strsplit(ete.id.tmp, ",")
# keep3 <- unlist( lapply(ete.id.tmp, function(ev){ !sum(as.numeric(ev)) == length(ev) }) )
# dge <- dge[keep3,]



dge$counts <- dge$counts + 1


dir.create("PLOTS", showWarnings=F, recursive=T)
name1 <- "MISO_p1"


write.table(dge$counts, paste0("PLOTS/dge_counts_",name1,".xls"), quote=F, sep="\t", row.names=T, col.names=T)


design <- model.matrix(~condition, data=metadata)
design


save(dge, design, file=paste0("PLOTS/", "dge_design",name1,".RData"))


source("/home/gosia/R/R_Multinomial_project/glmFitDM.R")



#######################################################
# run G2
#######################################################


name <- "MISO_G2_p1_no0ev_"


g2fit <- g2FitDM(dge)

g2lrt <- g2LRTDM(g2fit)



head(g2lrt$table)


save(g2lrt, file=paste0("PLOTS/", name, "g2lrt.RData"))

#######################################################
# run GLM 
#######################################################

name <- "MISO_optim_p1_"


glmfit <- glmFitDM(dge, design, optimization="optim")

glmlrt <- glmLRTDM(glmfit, coef=2)


head(glmlrt$table)


save(glmlrt, file=paste0("PLOTS/", name, "glmlrt.RData"))



#######################################################
# histograms of p-values
#######################################################

glmlrt <- g2lrt


pdf(paste0("PLOTS/", name, "hist_pvalues.pdf"))
hist(glmlrt$table$PValue, col="orange", breaks=50)
dev.off()


#######################################################
# trend: overall gene expression vs dispersion
#######################################################

gh0.tmp <-  mclapply(names(glmlrt$fit.full), function(g){
  
  if(!is.null(glmlrt$fit.full[[g]])){           
    
    meanY <- mean(rowSums(glmlrt$fit.full[[g]]$Y))
    gh0.0f <- mean(glmlrt$fit.full[[g]]$gh0[1])
    gh0.1f <- mean(glmlrt$fit.full[[g]]$gh0[4])
    gh0.0n <- mean(glmlrt$fit.null[[g]]$gh0[1])
    max.th.0n <- max(glmlrt$fit.null[[g]]$th)
    
    
    return(data.frame(gh0.0f=gh0.0f, gh0.1f=gh0.1f, gh0.0n=gh0.0n, meanY=meanY, max.th.0n=max.th.0n))
  }    
  return(data.frame(gh0.0f=NA, gh0.1f=NA, gh0.0n=NA, meanY=NA, max.th.0n=NA))
  
}, mc.cores=10)


gh0 <- do.call(rbind, gh0.tmp)

gh0 <- data.frame(GeneID=names(glmlrt$fit.full), gh0)

head(gh0)



table.gh0 <- gh0

table.gh0[!complete.cases(table.gh0), ]

table.gh0 <- table.gh0[complete.cases(table.gh0), ]
head(table.gh0)





pdf(paste0("PLOTS/", name, "TREND_meanY_vs_gh0.pdf"))

plot(table.gh0$meanY, table.gh0$gh0.0f, col=2, log="xy", main="log(meanY) vs log(gh0)")
plot(table.gh0$meanY , table.gh0$gh0.1f, col=3, log="xy")
plot(table.gh0$meanY , table.gh0$gh0.0n, col=4, log="xy")


hist(log10(table.gh0$gh0.0f),breaks=50)
hist(log10(table.gh0$gh0.1f),breaks=50)
hist(log10(table.gh0$gh0.0n),breaks=50)

dev.off()



pdf(paste0("PLOTS/", name, "TREND_gh0_vs_max_th.pdf"))

plot(table.gh0$max.th.0n, table.gh0$gh0.0n, col=4, log="xy")

dev.off()





#### plot likelihood from glm and g2

glm.lh <-  mclapply(names(glmlrt$fit.full), function(g){
  
  if(!is.null(glmlrt$fit.full[[g]])){           
    
    loglikh.f <- glmlrt$fit.full[[g]]$loglikh
    loglikh.n <- glmlrt$fit.null[[g]]$loglikh
    
    
    return(data.frame(loglikh.f=loglikh.f, loglikh.n=loglikh.n))
  }    
  return(data.frame(loglikh.f=NA, loglikh.n=NA))
  
}, mc.cores=10)



g2.lh <-  mclapply(names(g2lrt$fit.full), function(g){
  
  if(!is.null(g2lrt$fit.full[[g]])){           
    
    loglikh.f <- g2lrt$fit.full[[g]]$loglikh
    loglikh.n <- g2lrt$fit.null[[g]]$loglikh
    
    
    return(data.frame(loglikh.f=loglikh.f, loglikh.n=loglikh.n))
  }    
  return(data.frame(loglikh.f=NA, loglikh.n=NA))
  
}, mc.cores=10)



glm.lh <- do.call(rbind, glm.lh)

glm.lh <- data.frame(GeneID=names(glmlrt$fit.full), glm.lh)

head(glm.lh)



g2.lh <- do.call(rbind, g2.lh)

g2.lh <- data.frame(GeneID=names(g2lrt$fit.full), g2.lh)

head(g2.lh)


lh <- merge(glm.lh, g2.lh, by=1)
lh <- lh[complete.cases(lh), ]

sum(lh[,2] == lh[, 4])


pdf(paste0("PLOTS/", name, "LOGLH_glm_vs_g2.pdf"))

plot(lh[,2], lh[, 4], main="full")
plot(lh[,3], lh[, 5], main="null")

plot(-lh[,2], -lh[, 4], main="full", log="xy")
plot(-lh[,3], -lh[, 5], main="null", log="xy")


dev.off()


#### plot pvalues glm vs g2

pvs <- merge(glmlrt$table, g2lrt$table, by=1)
head(pvs)


write.table(pvs, paste0("PLOTS/", name, "TABLE_glm_vs_g2.xls"), sep="\t", row.names=F)


pdf(paste0("PLOTS/", name, "PVALUES_glm_vs_g2.pdf"))

plot()

dev.off()




#######################################################
# check 
#######################################################

genes.top <-  table.gh0[order(table.gh0$gh0.0n, decreasing=T), 1]


gene <- genes.top[16]
gene


glmlrt$fit.null[[gene]]

g2lrt$fit.null[[gene]]


glmlrt$fit.full[[gene]]
g2lrt$fit.full[[gene]]







































