#######################################################
# 
# Created 30 July 2014 / Last updated 31 July 2014
# 
# Update 31 July 2014:
# + comparison for stratified number of features/exons
# + count number of exons and number of DE exons

#######################################################
# BioC 2.14


setwd("/home/gosia/Multinomial_project/Simulations_hsapiens_V2")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


#######################################################
# run DM on featureCounts data with different filtering
#######################################################

source("/home/gosia/R/R_Multinomial_project/DM/g2FitDM_dirmult_v2.R")

out.dir <- "DM_filtering/fc"
dir.create(out.dir, showWarnings=F, recursive=T)

library(limma)

fc <- read.table("featureCounts/featureCounts.txt")

counts <- fc[,2:7]
colnames(counts) <- metadata$SampleName1
gene.id <- strsplit2(fc[,1], ":")[,1]
ete.id <- fc[,1]

library("edgeR")
dge.org <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))

name1 <- "fc"
dge <- dge.org
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)



######### standard filtering
name1 <- "fc"
dge <- dge.org
keep <- rowSums(cpm(dge)>1) >= 3
dge <- dge[keep,]
dge$counts[ dge$counts == 0 ] <- 1

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))



######### standard filtering and keep 0s
name1 <- "fc_keep0s"
dge <- dge.org
keep <- rowSums(cpm(dge)>1) >= 3
dge <- dge[keep,]

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit, mc.cores=30)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))




######### filtering 
name1 <- "fc_g0_s4"
dge <- dge.org
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]
dge$counts[ dge$counts == 0 ] <- 1

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit, mc.cores=30)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))




######### filtering / keep 0s
name1 <- "fc_g0_s4_keep0s"
dge <- dge.org
keep <- rowSums(cpm(dge) > 0) >= 4
dge <- dge[keep,]

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit, mc.cores=30)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))




######### NO filtering
name1 <- "fc_NOfiltering"
dge <- dge.org

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit, mc.cores=30)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))




#######################################################
# run DM on MISO data
#######################################################

source("/home/gosia/R/R_Multinomial_project/DM/g2FitDM_dirmult_v2.R")

out.dir <- "DM_filtering/miso"
dir.create(out.dir, showWarnings=F, recursive=T)

library(limma)

miso <- read.table("MISOCounts/MISO_counts_all_samps.txt", header = T)

counts <- miso[,2:7]
gene.id <- strsplit2(miso[,1], ":")[,1]
ete.id <- miso[,1]


library("edgeR")
dge.org <- DGEList(counts=counts, group = metadata$condition, genes=data.frame(gene.id=gene.id, ete.id=ete.id))

name1 <- "miso"
dge <- dge.org
write.table(data.frame(dge$genes,dge$counts), paste0(out.dir, "/dge_counts_",name1,"_NOT_FILTERED.xls"), quote=F, sep="\t", row.names=F, col.names=T)



##################### DM analysis 

### standard filtering / keep all events
name1 <- "miso"
dge <- dge.org
keep <- rowSums(cpm(dge)>1) >= 3
dge <- dge[keep,]
dge$counts[ dge$counts == 0 ] <- 1

# run DM G2
g2fit <- g2FitDM(dge)
g2lrt <- g2LRTDM(g2fit)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))




### standard filtering / keep all events / keep 0s
name1 <- "miso_keep0s"
dge <- dge.org
keep <- rowSums(cpm(dge)>1) >= 3
dge <- dge[keep,]

# run DM G2
g2fit <- g2FitDM(dge)
g2lrt <- g2LRTDM(g2fit)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))



###
name1 <- "miso_g0_s4"
dge <- dge.org
keep <- rowSums(cpm(dge)>0) >= 4
dge <- dge[keep,]
dge$counts[ dge$counts == 0 ] <- 1

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit, mc.cores=30)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))



### 
name1 <- "miso_g0_s4_keep0s"
dge <- dge.org
keep <- rowSums(cpm(dge)>0) >= 4
dge <- dge[keep,]

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit, mc.cores=30)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))







### 
name1 <- "miso_g0_s4_filter0ev"
dge <- dge.org
keep <- rowSums(cpm(dge)>0) >= 4
dge <- dge[keep,]
# for MISO: filter events with only 0s: 0,0,0,0,0
keep2 <- unlist( lapply(strsplit(strsplit2(dge$genes$ete.id, ":")[,2], ","), function(ev){ sum(as.numeric(ev)) > 0 }) )
dge <- dge[keep2,]

dge$counts[ dge$counts == 0 ] <- 1

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit, mc.cores=30)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))




### 
name1 <- "miso_g0_s4_filter0ev_keep0s"
dge <- dge.org
keep <- rowSums(cpm(dge)>0) >= 4
dge <- dge[keep,]
# for MISO: filter events with only 0s: 0,0,0,0,0
keep2 <- unlist( lapply(strsplit(strsplit2(dge$genes$ete.id, ":")[,2], ","), function(ev){ sum(as.numeric(ev)) > 0 }) )
dge <- dge[keep2,]

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit, mc.cores=30)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))








### 
name1 <- "miso_g0_s4_filter0ev1ev"
dge <- dge.org
keep <- rowSums(cpm(dge)>0) >= 4
dge <- dge[keep,]
# for MISO: filter events with only 0s: 0,0,0,0,0
keep2 <- unlist( lapply(strsplit(strsplit2(dge$genes$ete.id, ":")[,2], ","), function(ev){ sum(as.numeric(ev)) > 0 }) )
dge <- dge[keep2,]
# for MISO: filter events with only 1s: 1,1,1,1
keep3 <- unlist( lapply(strsplit(strsplit2(dge$genes$ete.id, ":")[,2], ","), function(ev){ !sum(as.numeric(ev)) == length(ev) }) )
dge <- dge[keep3,]

dge$counts[ dge$counts == 0 ] <- 1

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit, mc.cores=30)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))




### 
name1 <- "miso_g0_s4_filter0ev1ev_keep0s"
dge <- dge.org
keep <- rowSums(cpm(dge)>0) >= 4
dge <- dge[keep,]
# for MISO: filter events with only 0s: 0,0,0,0,0
keep2 <- unlist( lapply(strsplit(strsplit2(dge$genes$ete.id, ":")[,2], ","), function(ev){ sum(as.numeric(ev)) > 0 }) )
dge <- dge[keep2,]
# for MISO: filter events with only 1s: 1,1,1,1
keep3 <- unlist( lapply(strsplit(strsplit2(dge$genes$ete.id, ":")[,2], ","), function(ev){ !sum(as.numeric(ev)) == length(ev) }) )
dge <- dge[keep3,]

# run DM G2
g2fit <- g2FitDM(dge, mc.cores=30)
g2lrt <- g2LRTDM(g2fit, mc.cores=30)

write.table(g2lrt$table, paste0(out.dir, "/DM_",name1,"_results.xls"), quote=F, sep="\t", row.names=F, col.names=T)
save(dge, g2fit, g2lrt, file=paste0(out.dir, "/DM_",name1,"_results.RData"))





#######################################################
# load information about simulation 
#######################################################

simu_info.e <- read.table("Simu_info/true_exons_simulation.txt", header=T)
simu_info.g <- read.table("Simu_info/true_genes_simulation.txt", header=T)

head(simu_info.g)
head(simu_info.e)

simu_info1 <- simu_info.g[!duplicated(simu_info.g$Gene),c("Gene", "status", "num", "rpk.mean")]

###################### count number of defferentially spliced exons

simu_info.e.spl <- split(simu_info.e, simu_info.e$Gene)

num <- sapply(simu_info.e.spl, function(g){
  # g = simu_info.e.spl[[2]]
  
  n <- sum(g$status_exon==1)
  
  if(is.na(n))
    n <- 0
  
  return(n)
})

table(num, useNA = "always")

simu_info2 <- data.frame(Gene = names(num), num.diff.ex=num)


simu_info <- merge(simu_info1, simu_info2, by = "Gene", all=TRUE)
simu_info$num.diff.ex[is.na(simu_info$num.diff.ex)] <- 0

write.table(simu_info, "Simu_info/simu_info.xls", quote = F, sep = "\t", row.names = F, col.names = T)

###################### add info about the number of exons for each gene

fc <- read.table("featureCounts/featureCounts.txt")

library(limma)
gene.id <- strsplit2(fc[,1], ":")[,1]

num.ex <- as.data.frame(table(gene.id))
colnames(num.ex) <- c("Gene", "num.ex")


simu_info <- merge(simu_info, num.ex, by="Gene", all.x = TRUE)

table(simu_info$num.ex, useNA = "always") # NAs when overlapping genes

simu_info$num.ex[is.na(simu_info$num.ex)] <- 0

write.table(simu_info, "Simu_info/simu_info.xls", quote = F, sep = "\t", row.names = F, col.names = T)



###################### some stats

all( simu_info$num.ex >= simu_info$num.diff.ex )

## genes with no exons at all
sum( simu_info$status == 1 & simu_info$num.ex == 0 )

sum( simu_info$num.ex == 0 )


## genes with no DE exons
sum( simu_info$status == 1 & simu_info$num.diff.ex == 0 )





#######################################################
# calculate gene-level p-values for DEXSeq with Simes' method
#######################################################

dexseq.ex <- read.table("Results_from_Katarina/dexseq_exon_results.txt", header = T, stringsAsFactors = FALSE)

dexseq.ex <- dexseq.ex[!is.na(dexseq.ex$pvalue) & (duplicated(dexseq.ex$geneID, fromLast = T) | duplicated(dexseq.ex$geneID, fromLast = F) ) ,c("geneID", "pvalue")]

dexseq.ex.spl <- split(dexseq.ex, dexseq.ex$geneID)

dexseq.g <- data.frame(geneID=names(dexseq.ex.spl), pvalue=0, padjust=0)

dexseq.g$pvalue <- sapply(dexseq.ex.spl, function(g){
  
  n <- nrow(g)
  pv <- min(sort(g$pvalue,decreasing = FALSE)*n/(1:n))
  
})


dexseq.g$padjust <- p.adjust(dexseq.g$pvalue, method = "BH")


write.table(dexseq.g, "Results_from_Katarina/dexseq_gene_Simes_results.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



# comapre with DEXSeq results

dexseq.g.org <- read.table("Results_from_Katarina/dexseq_gene_results.txt", header = T, stringsAsFactors = FALSE)


d <- merge(dexseq.g.org, dexseq.g, by = 1, all=TRUE)

sum(is.na(d$padjust.x))
sum(is.na(d$padjust.y))
sum(is.na(d$padjust.y) && is.na(d$padjust.x))

d$padjust.x[is.na(d$padjust.x)] <- 1.1
d$padjust.y[is.na(d$padjust.y)] <- 1.1


pdf("Results_from_Katarina/dexseq_gene_results.pdf")
smoothScatter(d$padjust.x, d$padjust.y, nrpoints = Inf, pch=19, xlab="DEXSeq", ylab="Simes")
abline(a=0, b=1, col="red")
smoothScatter(log(d$padjust.x), log(d$padjust.y), nrpoints = Inf, xlab="DEXSeq", ylab="Simes")
abline(a=0, b=1, col="red")
dev.off()



#######################################################
# merge Diff Expr results with simu_info
#######################################################

results <- list()

####################### results produced by Katarina

rt <- read.table("Results_from_Katarina/dexseq_gene_results.txt", header = T, stringsAsFactors = F)
head(rt)
rt <- rt[,c("Gene", "padjust")]
colnames(rt) <- c("Gene",  "adjPValue_htseq_dexseq")

results[["htseq_dexseq"]] <- rt



rt <- read.table("Results_from_Katarina/dexseq_gene_Simes_results.txt", header = T, stringsAsFactors = F)
head(rt)
rt <- rt[,c("geneID","pvalue" ,"padjust")]
colnames(rt) <- c("Gene", "PValue_htseq_dexseq_simes", "adjPValue_htseq_dexseq_simes")

results[["htseq_dexseq_simes"]] <- rt


rt <- read.table("Results_from_Katarina/voom_gene_results.txt", header = T, stringsAsFactors = F)
head(rt)
rt <- rt[,c("gene.id","P.Value" ,"FDR")]
colnames(rt) <- c("Gene", "PValue_fc_voomex", "adjPValue_fc_voomex")

results[["fc_voomex"]] <- rt


## problems when adding this results (same genes have multiple p-values)
# rt <- read.table("Results_from_Katarina/cuffdiff_results.txt", header = T, stringsAsFactors = F)
# head(rt)
# rt <- rt[,c("ensembl_gene_id","p_value" ,"q_value")]
# colnames(rt) <- c("Gene", "PValue_cuffdiff", "adjPValue_cuffdiff")
# 
# results[["cuffdiff"]] <- rt


####################### DM results / different filtering for fc

filtering <- filtering.fc <- paste0("fc_", c("", "keep0s_", "g0_s4_", "g0_s4_keep0s_", "NOfiltering_"))

for(i in seq_len(length(filtering))){
  
  rt <- read.table(paste0("DM_filtering/fc/DM_",filtering[i],"results.xls"), header = T, stringsAsFactors = F)
  head(rt)
  rt <- rt[,c("GeneID","PValue" ,"FDR")]
  colnames(rt) <- c("Gene", paste0("PValue_",filtering[i],"DM"), paste0("adjPValue_",filtering[i],"DM"))
  
  results[[paste0(filtering[i],"DM")]] <- rt
  
  
}


####################### DM results / different filtering for miso



# filtering <- filtering.miso <- paste0("miso_", c("", "keep0s_", "g0_s4_", "g0_s4_keep0s_", "g0_s4_filter0ev_", "g0_s4_filter0ev_keep0s_", "g0_s4_filter0ev1ev_", "g0_s4_filter0ev1ev_keep0s_"))
# 
# for(i in seq_len(length(filtering))){
#   # i=2
#   cat(filtering[i], "\n")
#   rt <- read.table(paste0("DM_filtering/miso/DM_",filtering[i],"results.xls"), header = T, stringsAsFactors = F)
#   head(rt)
#   rt <- rt[,c("GeneID","PValue" ,"FDR")]
#   colnames(rt) <- c("Gene", paste0("PValue_",filtering[i],"DM"), paste0("adjPValue_",filtering[i],"DM"))
#   
#   results[[paste0(filtering[i],"DM")]] <- rt
#   
#   
# }


####################### DM results /  for htseq


# rt <- read.table("DM/htseq/DM_htseq_results.xls", header = T, stringsAsFactors = F)
# head(rt)
# rt <- rt[,c("GeneID","PValue" ,"FDR")]
# colnames(rt) <- c("Gene", "PValue_htseq_DM", "adjPValue_htseq_DM")
# 
# results[["htseq_DM"]] <- rt


######################### merge results into one table

out.dir <- "PLOTS_filtering"
dir.create(out.dir, showWarnings=F, recursive=T)

table <- simu_info
for(i in 1:length(results)){
  table <- merge(table, results[[i]], by="Gene", all=TRUE)  
}

table <- unique(table)
table(table$status, useNA="always")

table <- table[!is.na(table$status), ]


write.table(table, paste0(out.dir,"/Table_all_results.xls"),  quote = F, sep = "\t", row.names = F, col.names = T)



## check table 

dim(table[!complete.cases(table), ])
colnames(table)


# # double check when adding cuffdiff results 
# table <- table[!is.na(table$Gene), ]
# table.dupl <- table[duplicated(table$Gene, fromLast = T) | duplicated(table$Gene, fromLast = F), ]
# write.table(table.dupl, "PLOTS/Table_all_results_dupl.xls",  quote = F, sep = "\t", row.names = F, col.names = T)
# table <- unique(table)
# table.dupl <- table[duplicated(table$Gene, fromLast = T) | duplicated(table$Gene, fromLast = F), ]
# write.table(table.dupl, "PLOTS/Table_all_results_dupl.xls",  quote = F, sep = "\t", row.names = F, col.names = T)



sum(!is.na(table$PValue_fc_voomex))
sum(!is.na(table$PValue_fc_DM))


          



##############################################################################################################
# ->>>>>> load results
##############################################################################################################

out.dir <- "PLOTS_filtering"

table <- read.table(paste0(out.dir,"/Table_all_results.xls"), header = T, stringsAsFactors = F)

filtering.fc <- paste0("fc_", c("", "keep0s_", "g0_s4_", "g0_s4_keep0s_", "NOfiltering_"))
filtering.miso <- paste0("miso_", c("", "keep0s_", "g0_s4_", "g0_s4_keep0s_", "g0_s4_filter0ev_", "g0_s4_filter0ev_keep0s_", "g0_s4_filter0ev1ev_", "g0_s4_filter0ev1ev_keep0s_"))


library(RColorBrewer)

colors <- c("pink", "magenta", "dodgerblue3", "darkorchid")
n.colors <- c("htseq_dexseq", "htseq_dexseq_simes", "fc_voomex", "htseq_DM")                                  
names(colors) <- n.colors


###

n.colors.miso <- paste0(filtering.miso, "DM")   
colors.miso <- as.character( colorRampPalette(c("brown1","brown4"))(length(n.colors.miso)) )                   
names(colors.miso) <- n.colors.miso

n.colors.fc <- paste0(filtering.fc, "DM")   
colors.fc <- as.character( colorRampPalette(c("orange1","orange4"))(length(n.colors.fc)) )                   
names(colors.fc) <- n.colors.fc

###

# n.colors.miso <- paste0(filtering.miso, "DM")   
# colors.miso <- as.character( rainbow(length(n.colors.miso)) )                   
# names(colors.miso) <- n.colors.miso
# 
# n.colors.fc <- paste0(filtering.fc, "DM")   
# colors.fc <- as.character( rainbow (length(n.colors.fc)) )                   
# names(colors.fc) <- n.colors.fc


ramp <- colorRampPalette(brewer.pal(12,"Paired"))(12)

# ramp <- colorRampPalette(brewer.pal(12,"Set3"))(12)
# plot(1:12, pch=19, cex=4, col=ramp)
# dev.off()


n.colors.miso <- paste0(filtering.miso, "DM")   
colors.miso <- ramp[1:length(n.colors.miso) ]               
names(colors.miso) <- n.colors.miso

n.colors.fc <- paste0(filtering.fc, "DM")   
colors.fc <- ramp[1:length(n.colors.fc) ]                   
names(colors.fc) <- n.colors.fc


###

colors.org <- c(colors, colors.miso, colors.fc)
n.colors.org <- names(colors.org)
name <- ""


colors.org <- colors.org[names(results)]
n.colors.org <- names(colors.org)


colors <- colors.org
n.colors <- n.colors.org


#######################################################
# histograms of p-values
#######################################################


colors <- colors.org[-1]
n.colors <- n.colors.org[-1]

pdf(paste0(out.dir, "/", name, "hist_pvalues.pdf"))


for(i in seq_len(length(colors))){
  # i=4
  cat(i, "\n")
  hist(table[, paste0("PValue_", n.colors[i]) ], col=colors[i], breaks=50, main=n.colors[i], cex.main=2, cex.lab=1.5, cex.axis = 1.5)
  
}

dev.off()



#######################################################
# ROCx plot  --- my function
#######################################################




ROCx <- function(pvs, apvs, status, add=FALSE, ...){
  # status==1, means true DS, status!=1, means not DS
  
  status.org <- status
  P.org <- sum(status.org==1)
  N.org <- sum(status.org==0)
  
  NAs <- is.na(pvs)
  pvs <- pvs[!NAs]
  status <- status[!NAs]
  
  P <- sum(status==1)
  N <- sum(status==0)
  
  ord <- order(pvs, decreasing = FALSE, na.last = TRUE)
  status <- status[ord]
  pvs <- pvs[ord]
  
  TPRv <- cumsum(status) / P
  FPRv <- cumsum(!status) / N
  
  
  if(add==FALSE){
    plot(FPRv, TPRv, xlab="False positive rate", ylab="True positive rate", type="l", ...)
  }else{   
    lines(FPRv, TPRv, type="l", ...)
  }
  
  TPR <- sum(apvs[status.org==1] < 0.05 , na.rm = T) / P.org
  points((approx(TPRv, FPRv, xout=TPR)$y), TPR, ...)  
  
}



ROCx_notNorm <- function(pvs, apvs, status, add=FALSE, ...){
  # status==1, means true DS, status!=1, means not DS

  status.org <- status
  P <- sum(status.org==1)
  N <- sum(status.org==0)
  
  NAs <- is.na(pvs)
  pvs <- pvs[!NAs]
  status <- status[!NAs]
  
  ord <- order(pvs, decreasing = FALSE, na.last = TRUE)
  status <- status[ord]
  pvs <- pvs[ord]
  
  TPRv <- cumsum(status) / P
  FPRv <- cumsum(!status) / N

  
  if(add==FALSE){
    plot(FPRv, TPRv, xlab="False positive rate", ylab="True positive rate", type="l", ...)
  }else{   
    lines(FPRv, TPRv, type="l", ...)
  }
  
  TPR <- sum(apvs[status.org==1] < 0.05 , na.rm = T) / P
  points((approx(TPRv, FPRv, xout=TPR)$y), TPR, ...)  
  
}




#######################################################
# generate FINAL ROCx plots 
#######################################################


colors <- colors.org
n.colors <- n.colors.org


############################ repeat previous plots

pdf(paste0(out.dir, "/",name,"rocX_all_DM_notNorm.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col=colors["fc_DM"], pch=4, xlim=c(0,1), ylim=c(0,1), lwd=3, cex=3, cex.lab=1.5, cex.axis=1.5)

pvs <- table[,c("PValue_htseq_DM")]
apvs <- table[,c("adjPValue_htseq_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["htseq_DM"], pch=4,  lwd=3, cex=3)

pvs <- table[,c("PValue_miso_DM")]
apvs <- table[,c("adjPValue_miso_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["miso_DM"], pch=4, lwd=3, cex=3)

legend("bottomright", c("fc_DM", "htseq_DM", "miso_DM"), col=colors[c("fc_DM", "htseq_DM", "miso_DM")], lty=1, lwd=3, cex=1.5)
dev.off()




#methods <- c("fc_DM", "htseq_DM", "fc_voomex", "htseq_dexseq_simes")
methods <- c("fc_DM", "fc_voomex", "htseq_dexseq_simes")

pdf(paste0(out.dir, "/",name,"rocX_all_notNorm.pdf"))

pvs <- table[,paste0("PValue_", methods[1])]
apvs <- table[,paste0("adjPValue_", methods[1])]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col=colors[methods[1]], pch=4, xlim=c(0,1), ylim=c(0,1), lwd=3, cex=3, cex.lab=1.5, cex.axis=1.5)

for(i in methods[-1]){
  pvs <- table[,paste0("PValue_", i)]
  apvs <- table[,paste0("adjPValue_", i)]
  status <- table$status
  ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors[i], pch=4, xlim=c(0,1), lwd=4, cex=3)
}

legend("bottomright", methods , col=colors[methods], lty=1, lwd=4, cex=1, bty="n")
dev.off()






############################ filtered data


pdf(paste0(out.dir, "/",name,"rocX_fc_filt_notNorm.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col=colors["fc_DM"], pch=4, ylim=c(0,1), xlim=c(0,0.7), lwd=4, cex=3, cex.lab=1.5, cex.axis=1.5)

pvs <- table[,c("PValue_fc_keep0s_DM")]
apvs <- table[,c("adjPValue_fc_keep0s_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["fc_keep0s_DM"], pch=4, xlim=c(0,1), lwd=4, cex=3, lty=2)

pvs <- table[,c("PValue_fc_g0_s4_DM")]
apvs <- table[,c("adjPValue_fc_g0_s4_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["fc_g0_s4_DM"], pch=4, xlim=c(0,1), lwd=4, cex=3)

pvs <- table[,c("PValue_fc_g0_s4_keep0s_DM")]
apvs <- table[,c("adjPValue_fc_g0_s4_keep0s_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["fc_g0_s4_keep0s_DM"], pch=4, xlim=c(0,1), lwd=4, cex=3, lty=2)

legend("bottomright", c("fc_DM", "fc_keep0s_DM", "fc_g0_s4_DM", "fc_g0_s4_keep0s_DM"), col=colors[c("fc_DM", "fc_keep0s_DM", "fc_g0_s4_DM", "fc_g0_s4_keep0s_DM")], lwd=4, cex=1.5, lty=c(1, 2, 1, 2))
dev.off()




pdf(paste0(out.dir, "/",name,"rocX_fc_filt_notNorm2.pdf"))

pvs <- table[,c("PValue_fc_DM")]
apvs <- table[,c("adjPValue_fc_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col=colors["fc_DM"], pch=4, ylim=c(0.4,0.9), xlim=c(0,0.7), lwd=4, cex=3, cex.lab=1.5, cex.axis=1.5)

pvs <- table[,c("PValue_fc_keep0s_DM")]
apvs <- table[,c("adjPValue_fc_keep0s_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["fc_keep0s_DM"], pch=4, xlim=c(0,1), lwd=4, cex=3, lty=2)

pvs <- table[,c("PValue_fc_g0_s4_DM")]
apvs <- table[,c("adjPValue_fc_g0_s4_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["fc_g0_s4_DM"], pch=4, xlim=c(0,1), lwd=4, cex=3)

pvs <- table[,c("PValue_fc_g0_s4_keep0s_DM")]
apvs <- table[,c("adjPValue_fc_g0_s4_keep0s_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["fc_g0_s4_keep0s_DM"], pch=4, xlim=c(0,1), lwd=4, cex=3, lty=2)


pvs <- table[,c("PValue_fc_NOfiltering_DM")]
apvs <- table[,c("adjPValue_fc_NOfiltering_DM")]
status <- table$status
ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors["fc_NOfiltering_DM"], pch=4, xlim=c(0,1), lwd=4, cex=3, lty=2)


legend("bottomright", c("fc_DM", "fc_keep0s_DM", "fc_g0_s4_DM", "fc_g0_s4_keep0s_DM", "fc_NOfiltering_DM"), col=colors[c("fc_DM", "fc_keep0s_DM", "fc_g0_s4_DM", "fc_g0_s4_keep0s_DM", "fc_NOfiltering_DM")], lwd=4, cex=1.5, lty=c(1, 2, 1, 2))
dev.off()




methods <- c("miso_DM", "miso_keep0s_DM", "miso_g0_s4_DM", "miso_g0_s4_keep0s_DM", "miso_g0_s4_filter0ev_DM", "miso_g0_s4_filter0ev_keep0s_DM", "miso_g0_s4_filter0ev1ev_DM", "miso_g0_s4_filter0ev1ev_keep0s_DM")

pdf(paste0(out.dir, "/",name,"rocX_miso_filt_notNorm.pdf"))

pvs <- table[,paste0("PValue_", methods[1])]
apvs <- table[,paste0("adjPValue_", methods[1])]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col=colors[methods[1]], pch=4, xlim=c(0,0.15), ylim=c(0.5,0.8), lwd=3, cex=3, cex.lab=1.5, cex.axis=1.5)

for(i in methods[-1]){
  pvs <- table[,paste0("PValue_", i)]
  apvs <- table[,paste0("adjPValue_", i)]
  status <- table$status
  ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors[i], pch=4, xlim=c(0,1), lwd=4, cex=3)
}

legend("bottomright", methods , col=colors[methods], lty=1, lwd=4, cex=1, bty="n")
dev.off()



############################ filtered data with previous methods


methods <- c("htseq_dexseq_simes", "fc_voomex", "fc_DM", "fc_g0_s4_keep0s_DM")

pdf(paste0(out.dir, "/",name,"rocX_all_filt_notNorm.pdf"))

pvs <- table[,paste0("PValue_", methods[1])]
apvs <- table[,paste0("adjPValue_", methods[1])]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col=colors[methods[1]], pch=4, ylim=c(0,1), xlim=c(0,1), lwd=4, cex=3, cex.lab=1.5, cex.axis=1.5)

for(i in methods[-1]){
  pvs <- table[,paste0("PValue_", i)]
  apvs <- table[,paste0("adjPValue_", i)]
  status <- table$status
  ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors[i], pch=4, xlim=c(0,1), lwd=4, cex=3)
}

legend("bottomright", methods , col=colors[methods], lty=1, lwd=4, cex=1.5)
dev.off()



methods <- c("htseq_dexseq_simes", "fc_voomex", "fc_DM", "fc_g0_s4_DM")

pdf(paste0(out.dir, "/",name,"rocX_all_filt_notNorm2.pdf"))

pvs <- table[,paste0("PValue_", methods[1])]
apvs <- table[,paste0("adjPValue_", methods[1])]
status <- table$status
ROCx_notNorm(pvs, apvs, status, col=colors[methods[1]], pch=4, ylim=c(0,1), xlim=c(0,1), lwd=4, cex=3, cex.lab=1.5, cex.axis=1.5)

for(i in methods[-1]){
  pvs <- table[,paste0("PValue_", i)]
  apvs <- table[,paste0("adjPValue_", i)]
  status <- table$status
  ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors[i], pch=4, xlim=c(0,1), lwd=4, cex=3)
}

legend("bottomright", methods , col=colors[methods], lty=1, lwd=4, cex=1.5)
dev.off()



#######################################################
# compare the results between add1s and keep0s 
#######################################################


filtering.fc <- paste0("fc_", c("", "keep0s_", "g0_s4_", "g0_s4_keep0s_", "NOfiltering_"))
filtering.miso <- paste0("miso_", c("", "keep0s_", "g0_s4_", "g0_s4_keep0s_", "g0_s4_filter0ev_", "g0_s4_filter0ev_keep0s_", "g0_s4_filter0ev1ev_", "g0_s4_filter0ev1ev_keep0s_"))


##### fc
sum(! table[, paste0("PValue_", filtering.fc[1], "DM")] == table[, paste0("PValue_", filtering.fc[2], "DM")], na.rm = TRUE)

sum(! table[, paste0("PValue_", filtering.fc[3], "DM")] == table[, paste0("PValue_", filtering.fc[4], "DM")], na.rm = TRUE)

sum(is.na(table[, paste0("PValue_", filtering.fc[3], "DM")]))
sum(is.na(table[, paste0("PValue_", filtering.fc[4], "DM")]))


##### miso

sum(! table[, paste0("PValue_", filtering.miso[1], "DM")] == table[, paste0("PValue_", filtering.miso[2], "DM")], na.rm = TRUE)

sum(! table[, paste0("PValue_", filtering.miso[3], "DM")] == table[, paste0("PValue_", filtering.miso[4], "DM")], na.rm = TRUE)

sum(! table[, paste0("PValue_", filtering.miso[5], "DM")] == table[, paste0("PValue_", filtering.miso[6], "DM")], na.rm = TRUE)

sum(! table[, paste0("PValue_", filtering.miso[7], "DM")] == table[, paste0("PValue_", filtering.miso[8], "DM")], na.rm = TRUE)



#######################################################
# TPR vs acieved FDR --- my functions
#######################################################



TPRvsFDR <- function(status, apvs, FDR.cut.off=c(0.01, 0.05, 0.1), pch=c( 22, 23, 24), col="red", cex.axis=1.5, cex=2.5, add=FALSE, ...){
  
  apvs[is.na(apvs)] <- 1
  
  n <- length(status)
  q <- length(FDR.cut.off)
  TPR <- rep(0, q)
  FDR <- rep(0, q)
  
  for(i in 1:q){
    # i=1
    status.est <- as.numeric(apvs < FDR.cut.off[i])
    
    TP <- sum(status==1 & status.est==1)
    FP <- sum(status==0 & status.est==1)
    FN <- sum(status==1 & status.est==0)
    
    TPR[i] <- TP/(TP+FN)
    FDR[i] <- FP/(FP+TP)
    
  }  
  tf <- cbind( FDR , TPR)
  
  bg <- rep(col,q)
  bg[FDR > FDR.cut.off] <- "white"
  
  if(add==FALSE){
   
    plot(1:2, type="n", xlab="Achieved FDR", ylab="TPR", xaxt="n", col=col, cex.axis=cex.axis, ...)
    axis(side=1, at=(2:10)/10, labels=(2:10)/10, las=1, col.ticks="grey", col.axis="grey", cex.axis=cex.axis)
    axis(side=1, at=FDR.cut.off, labels=FDR.cut.off, las=1, cex.axis=cex.axis)

    for(i in 1:q)
    lines(rep(FDR.cut.off[i], 50), seq(-0.1,1.1,length.out=50), type="b", pch=pch[i], cex=0.5, bg=1)
    
    lines(tf, type="l", col=col, ...)
    points(tf, pch=pch, bg=bg, col=col, cex=cex, ...)
    
  }else{    
    lines(tf, type="l", col=col, ...)
    points(tf, pch=pch, bg=bg, col=col, cex=cex,...)
  }
  
  
  
}



TPRvsFDRextend <- function(status, apvs, FDR.cut.off=c(0.01, 0.05, 0.1), pch=c( 22, 23, 24), col="red", cex.axis=1.5, cex=2.5, add=FALSE, ...){
  
  apvs[is.na(apvs)] <- 1
  
  FDR.cut.off.org <- FDR.cut.off
  FDR.cut.off <- c(FDR.cut.off.org, seq(max(FDR.cut.off.org), 1, 0.01))

  n <- length(status)
  q.org <- length(FDR.cut.off.org)
  q <- length(FDR.cut.off)
  TPR <- rep(0, q)
  FDR <- rep(0, q)
  
  for(i in 1:q){
    # i=1
    status.est <- as.numeric(apvs < FDR.cut.off[i])
    
    TP <- sum(status==1 & status.est==1)
    FP <- sum(status==0 & status.est==1)
    FN <- sum(status==1 & status.est==0)
    
    TPR[i] <- TP/(TP+FN)
    FDR[i] <- FP/(FP+TP)
    
  }  
  tf <- cbind( FDR , TPR)
  
  bg <- rep(col,q.org)
  bg[FDR[1:q.org] > FDR.cut.off.org] <- "white"
  
  if(add==FALSE){
    
    plot(1:2, type="n", xlab="Achieved FDR", ylab="TPR", xaxt="n", col=col, cex.axis=cex.axis, ...)
    axis(side=1, at=(2:10)/10, labels=(2:10)/10, las=1, col.ticks="grey", col.axis="grey", cex.axis=cex.axis)
    axis(side=1, at=FDR.cut.off.org, labels=FDR.cut.off.org, las=1, cex.axis=cex.axis)
    
    for(i in 1:q.org)
      lines(rep(FDR.cut.off[i], 50), seq(-0.1,1.1,length.out=50), type="b", pch=pch[i], cex=0.5, bg=1)
    
    lines(tf, type="l", col=col, ...)
    points(tf[1:q.org,], pch=pch, bg=bg, col=col, cex=cex, ...)
    
  }else{    
    lines(tf, type="l", col=col, ...)
    points(tf[1:q.org,], pch=pch, bg=bg, col=col, cex=cex,...)
  }
  
  
  
}


#######################################################
# TPR vs acieved FDR plots
#######################################################


# methods <- c("htseq_dexseq_simes", "fc_voomex", "htseq_DM", "fc_DM")
methods <- c("htseq_dexseq_simes", "fc_voomex", "fc_DM")

pdf(paste0(out.dir, "/", name, "TPRachievedFDR_all.pdf"))

apvs <- table[,paste0("adjPValue_", methods[1])]
status <- table$status
TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.5), ylim=c(0,1), lwd=4, cex=2.5, cex.lab=1.5, cex.axis=1.5)

for( i in methods[-1]){
  apvs <- table[,paste0("adjPValue_", i)]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
}

legend("bottomright", methods, col=colors[methods], lty=1, lwd=4, cex=1.5, bty="n")
dev.off()




methods <- c("htseq_dexseq","htseq_dexseq_simes", "fc_voomex", "fc_DM", "fc_g0_s4_keep0s_DM")

pdf(paste0(out.dir, "/", name, "TPRachievedFDR_all_filt.pdf"))

apvs <- table[,paste0("adjPValue_", methods[1])]
status <- table$status
TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.7), ylim=c(0,1), lwd=4, cex=2.5, cex.lab=1.5, cex.axis=1.5)

for( i in methods[-1]){
  apvs <- table[,paste0("adjPValue_", i)]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
}

legend("bottomright", methods, col=colors[methods], lty=1, lwd=4, cex=1.5, bg="white" )
dev.off()



methods <- c("htseq_dexseq","htseq_dexseq_simes", "fc_voomex", "fc_DM", "fc_g0_s4_DM")

pdf(paste0(out.dir, "/", name, "TPRachievedFDR_all_filt2.pdf"))

apvs <- table[,paste0("adjPValue_", methods[1])]
status <- table$status
TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.7), ylim=c(0,1), lwd=4, cex=2.5, cex.lab=1.5, cex.axis=1.5)

for( i in methods[-1]){
  apvs <- table[,paste0("adjPValue_", i)]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
}

legend("bottomright", methods, col=colors[methods], lty=1, lwd=4, cex=1.5, bg="white" )
dev.off()



################################### TPRvsFDRextend



methods <- c("htseq_dexseq","htseq_dexseq_simes", "fc_voomex", "fc_DM", "fc_g0_s4_keep0s_DM")

pdf(paste0(out.dir, "/", name, "TPRachievedFDRextend_all_filt.pdf"))

apvs <- table[,paste0("adjPValue_", methods[1])]
status <- table$status
TPRvsFDRextend(status, apvs, col=colors[methods[1]], xlim=c(0,1), ylim=c(0,1), lwd=4, cex=1, cex.lab=1.5, cex.axis=1.5)

for( i in methods[-1]){
  apvs <- table[,paste0("adjPValue_", i)]
  status <- table$status
  TPRvsFDRextend(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2)
}

legend("bottomright", methods, col=colors[methods], lty=1, lwd=4, cex=1.5, bg="white" )
dev.off()



methods <- c("htseq_dexseq","htseq_dexseq_simes", "fc_voomex", "fc_DM", "fc_g0_s4_DM")

pdf(paste0(out.dir, "/", name, "TPRachievedFDRextend_all_filt2.pdf"))

apvs <- table[,paste0("adjPValue_", methods[1])]
status <- table$status
TPRvsFDRextend(status, apvs, col=colors[methods[1]], xlim=c(0,1), ylim=c(0,1), lwd=4, cex=1, cex.lab=1.5, cex.axis=1.5)

for( i in methods[-1]){
  apvs <- table[,paste0("adjPValue_", i)]
  status <- table$status
  TPRvsFDRextend(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2)
}

legend("bottomright", methods, col=colors[methods], lty=1, lwd=4, cex=1.5, bg="white" )
dev.off()



#######################################################
# venn diagrams --- my function
#######################################################
library(VennDiagram)

FDR.cutoff <- 0.05


plot.venn <- function(venne.genes, colors, venn.methods, name="", cex=1, cat.cex=1.4, lwd=2, lty=1, alpha=0.5, margin=0.05){
  
  colors.col <- colors[venn.methods]
  
  if(length(venn.methods)==4)
    colors.col <- colors.col[c(1, 3, 4 ,2)]
  
  venn.d <- venn.diagram(venne.genes[venn.methods], filename=NULL, fill = colors[venn.methods], col=colors.col, cex=cex, cat.cex=cat.cex, lwd=lwd, lty=lty, alpha=alpha, margin=margin)
  
  pdf(paste0(out.dir, "/", name, "Venn_Diagr.pdf"))
  grid.draw(venn.d)
  dev.off()
  
  
}


#######################################################
# venn diagrams
#######################################################

## TP as a separate circle

venne.genes <- list()

for(i in n.colors){ 
  venne.genes[[i]] <- na.omit(table[ table[,paste0("adjPValue_", i)] < FDR.cutoff, "Gene"]) 
}

venne.genes$True <- table[table$status == 1, "Gene"]




plot.venn(venne.genes, colors=c(colors[c("fc_DM", "fc_voomex", "htseq_dexseq_simes")], True="grey"), venn.methods = c("fc_DM", "fc_voomex", "htseq_dexseq_simes", "True"), name=paste0(name, "All_"), margin=0.1, cat.cex=1.4, cex=1.7)



plot.venn(venne.genes, colors=c(colors[c("fc_g0_s4_keep0s_DM", "fc_voomex", "htseq_dexseq_simes")], True="grey"), venn.methods = c("fc_g0_s4_keep0s_DM", "fc_voomex", "htseq_dexseq_simes", "True"), name=paste0(name, "All_filt_"), margin=0.1, cat.cex=1.4, cex=1.7)


plot.venn(venne.genes, colors=c(colors[c("fc_g0_s4_DM", "fc_voomex", "htseq_dexseq_simes")], True="grey"), venn.methods = c("fc_g0_s4_DM", "fc_voomex", "htseq_dexseq_simes", "True"), name=paste0(name, "All_filt2_"), margin=0.1, cat.cex=1.4, cex=1.7)





##############################################################################################################
# stratified results 
##############################################################################################################

table.org <- table

c1 <- table(table.org[,c("num.ex", "status")])

pdf(paste0(out.dir, "/",name,"barplot_STRAT.pdf"), h=7, w=14)
barplot(t(c1)+1, main="Genes distribution by nr. of exons",xlab="Number of Exons", col=c("darkblue","red"), legend = colnames(c1), xlim=c(0, 100))
dev.off()


c2 <- table(table.org[table.org$status==1,c("num.ex", "num.diff.ex")])

# library(gplots)
# library(RColorBrewer)
# 
# pdf(paste0(out.dir, "/",name,"barplot_STRAT2.pdf"))
# heatmap.2(c2, dendrogram="none", trace="none", xlim=c(0,30), ylim=c(0,30))
# dev.off()


strat <- c(2, 6, 10, 14, 18, 22, 30, 38, 50, 400)


##################### ROCx plots


methods <- c("htseq_dexseq_simes", "fc_voomex", "fc_g0_s4_keep0s_DM")
pdf(paste0(out.dir, "/",name,"rocX_STRAT_all_filt_notNorm.pdf"))

for(j in 1:(length(strat)-1) ){
  
  table <- table.org[ table.org$num.ex >= strat[j] & table.org$num.ex < strat[j+1], ]
  
  pvs <- table[,paste0("PValue_", methods[1])]
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  ROCx_notNorm(pvs, apvs, status, col=colors[methods[1]], pch=4, ylim=c(0,1), xlim=c(0,1), lwd=4, cex=3, cex.lab=1.5, cex.axis=1.5, main=paste0( "Number of exons: ", paste0( strat[j]:( strat[j+1] -1 ), collapse = ", ") ) )
  
  for(i in methods[-1]){
    pvs <- table[,paste0("PValue_", i)]
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors[i], pch=4, xlim=c(0,1), lwd=4, cex=3)
  }
  
  legend("bottomright", methods , col=colors[methods], lty=1, lwd=4, cex=1.5)
  
}

dev.off()




methods <- c("htseq_dexseq_simes", "fc_voomex", "fc_g0_s4_DM")
pdf(paste0(out.dir, "/",name,"rocX_STRAT_all_filt_notNorm2.pdf"))

for(j in 1:(length(strat)-1) ){
  
  table <- table.org[ table.org$num.ex >= strat[j] & table.org$num.ex < strat[j+1], ]
  
  pvs <- table[,paste0("PValue_", methods[1])]
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  ROCx_notNorm(pvs, apvs, status, col=colors[methods[1]], pch=4, ylim=c(0,1), xlim=c(0,1), lwd=4, cex=3, cex.lab=1.5, cex.axis=1.5, main=paste0( "Number of exons: ", paste0( strat[j]:( strat[j+1] -1 ), collapse = ", ") ) )
  
  for(i in methods[-1]){
    pvs <- table[,paste0("PValue_", i)]
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    ROCx_notNorm(pvs, apvs, status,  add=TRUE, col=colors[i], pch=4, xlim=c(0,1), lwd=4, cex=3)
  }
  
  legend("bottomright", methods , col=colors[methods], lty=1, lwd=4, cex=1.5)
  
}

dev.off()


##################### TPRvsFDR plots

methods <- c("htseq_dexseq","htseq_dexseq_simes", "fc_voomex", "fc_g0_s4_keep0s_DM")
pdf(paste0(out.dir, "/", name, "TPRachievedFDR_STRAT_all_filt.pdf"))

for(j in 1:(length(strat)-1) ){
  
  table <- table.org[ table.org$num.ex >= strat[j] & table.org$num.ex < strat[j+1], ]
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.7), ylim=c(0,1), lwd=4, cex=2.5, cex.lab=1.5, cex.axis=1.5, main=paste0("Number of exons: ", paste0(strat[j]:( strat[j+1] -1 ), collapse = ", ")) )  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
  } 
  legend("bottomright", methods, col=colors[methods], lty=1, lwd=4, cex=1.5, bg="white" ) 
}
dev.off()



methods <- c("htseq_dexseq","htseq_dexseq_simes", "fc_voomex", "fc_g0_s4_DM")
pdf(paste0(out.dir, "/", name, "TPRachievedFDR_STRAT_all_filt2.pdf"))

for(j in 1:(length(strat)-1) ){
  
  table <- table.org[ table.org$num.ex >= strat[j] & table.org$num.ex < strat[j+1], ]
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.7), ylim=c(0,1), lwd=4, cex=2.5, cex.lab=1.5, cex.axis=1.5, main=paste0("Number of exons: ", paste0(strat[j]:( strat[j+1] -1 ), collapse = ", ")) )
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
  }  
  legend("bottomright", methods, col=colors[methods], lty=1, lwd=4, cex=1.5, bg="white" ) 
}
dev.off()
























































































