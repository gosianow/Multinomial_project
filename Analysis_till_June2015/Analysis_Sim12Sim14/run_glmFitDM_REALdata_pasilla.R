#######################################################
# run on REAL data
#######################################################


setwd("/home/gosia/Multinomial_project/brooks_pasilla")


# create metadata file
metadata <- read.table("metadata.xls", header=T, sep="\t")[1:6,]
head(metadata)

library("edgeR")
library("DEXSeq")

# run read.HTSeqCounts to create ExonCountSet object
metadata$condition <- as.factor(metadata$condition)


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
keep <- rowSums(cpm(dge)>1)>=3
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


# dir.create("PLOTS", showWarnings=F, recursive=T)


dge$counts <- dge$counts + 1


write.table(dge$counts, "PLOTS/dge_counts_MISOp1_no0ev.xls", quote=F, sep="\t", row.names=T, col.names=T)

design <- model.matrix(~condition, data=metadata)
design


save(dge, design, file=paste0("PLOTS/", "dge_designp1_no0ev.RData"))


source("/home/gosia/R/R_Multinomial_project/dirmult.R")
source("/home/gosia/R/R_Multinomial_project/glmFitDM.R")



#######################################################
# run G2
#######################################################


name <- "G2_p1_no0ev_"


g2fit <- g2FitDM(dge)

g2lrt <- g2LRTDM(g2fit)



head(g2lrt$table)


save(g2lrt, file=paste0("PLOTS/", name, "g2lrt.RData"))



#######################################################
# histograms of p-values
#######################################################

# glmlrt <- g2lrt


pdf(paste0("PLOTS/", name, "hist_pvalues_MISO.pdf"))
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


hist(log10(table.gh0$gh0.0f))
hist(log10(table.gh0$gh0.1f))
hist(log10(table.gh0$gh0.0n))

dev.off()



pdf(paste0("PLOTS/", name, "TREND_gh0_vs_max_th.pdf"))

plot(table.gh0$max.th.0n, table.gh0$gh0.0n, col=4, log="xy")

dev.off()



#######################################################
# check 
#######################################################

genes.top <-  table.gh0[table.gh0$gh0.0n > 1e+6, 1]

gene <- genes.top[6]


#glmlrt$fit.full[[gene]]
glmlrt$fit.null[[gene]]

























