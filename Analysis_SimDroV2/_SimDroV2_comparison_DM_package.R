#######################################################
# 
# Created 23 Oct 2014 

# comprare the DM version 5 commonDispersion and adjustement performance with other methods and DM dirmult

# Update 23 Oct 2014:



#######################################################
# BioC 2.14

setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2")




#######################################################
# calculate gene-level p-values for DEXSeq with Simes' method
#######################################################


DEXSeqSimes <- function(dexseqExonFile, dexseqGeneFile, dexseqGeneFileSimes){
  
  dexseq.ex <- read.table(dexseqExonFile, header = T, stringsAsFactors = FALSE)
  head(dexseq.ex)
  
  dexseq.ex <- dexseq.ex[!is.na(dexseq.ex$pvalue), c("geneID", "pvalue")]
  
  dexseq.ex.spl <- split(dexseq.ex, dexseq.ex$geneID)
  
  dexseq.g <- data.frame(geneID=names(dexseq.ex.spl), pvalue=0, padjust=0)
  
  dexseq.g$pvalue <- sapply(dexseq.ex.spl, function(g){
    
    n <- nrow(g)
    pv <- min(sort(g$pvalue,decreasing = FALSE)*n/(1:n))
    
  })
  
  dexseq.g$padjust <- p.adjust(dexseq.g$pvalue, method = "BH")
  
  write.table(dexseq.g, dexseqGeneFileSimes, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  # comapre with DEXSeq results
  
  dexseq.g.org <- read.table(dexseqGeneFile, header = T, stringsAsFactors = FALSE)
  
  
  d <- merge(dexseq.g.org, dexseq.g, by = 1, all=TRUE)
  
  sum(is.na(d$padjust.x))
  sum(is.na(d$padjust.y))
  sum(is.na(d$padjust.y) && is.na(d$padjust.x))
  
  d$padjust.x[is.na(d$padjust.x)] <- 1.1
  d$padjust.y[is.na(d$padjust.y)] <- 1.1
  
  
  pdf(gsub(".txt", ".pdf", dexseqGeneFile, fixed = TRUE))
  plot(d$padjust.x, d$padjust.y, pch=19, xlab="DEXSeq", ylab="Simes")
  abline(a=0, b=1, col="red")
  smoothScatter(d$padjust.x, d$padjust.y, nrpoints = Inf, pch=19, xlab="DEXSeq", ylab="Simes")
  abline(a=0, b=1, col="red")
  smoothScatter(log(d$padjust.x), log(d$padjust.y), nrpoints = Inf, xlab="DEXSeq", ylab="Simes", main="log scale")
  abline(a=0, b=1, col="red")
  dev.off()
  
  
  
}


dexseqExonFile <- "Results_from_Katarina/dexseq_exon_version_1.10.8_results.txt"
dexseqGeneFile <- "Results_from_Katarina/dexseq_gene_version_1.10.8_results.txt"
dexseqGeneFileSimes <- "Results_from_Katarina/dexseq_gene_Simes_version_1.10.8_results.txt"


DEXSeqSimes(dexseqExonFile, dexseqGeneFile, dexseqGeneFileSimes)


dexseqExonFile <- "Results_from_Katarina/dexseq_exon_featureCounts_version_1.10.8_results.txt"
dexseqGeneFile <- "Results_from_Katarina/dexseq_gene_featureCounts_version_1.10.8_results.txt"
dexseqGeneFileSimes <- "Results_from_Katarina/dexseq_gene_Simes_featureCounts_version_1.10.8_results.txt"


DEXSeqSimes(dexseqExonFile, dexseqGeneFile, dexseqGeneFileSimes)


dexseqExonFile <- "Results_from_Katarina/dexseq_exon_Gordon_featureCounts_version_1.10.8_results.txt"
dexseqGeneFile <- "Results_from_Katarina/dexseq_gene_Gordon_featureCounts_version_1.10.8_results.txt"
dexseqGeneFileSimes <- "Results_from_Katarina/dexseq_gene_Simes_Gordon_featureCounts_version_1.10.8_results.txt"


DEXSeqSimes(dexseqExonFile, dexseqGeneFile, dexseqGeneFileSimes)




#######################################################
# merge Diff Expr results 
#######################################################

results <- list()

####################### results produced by Katarina

rt <- read.table("Results_from_Katarina/dexseq_gene_version_1.10.8_results.txt", header = T, stringsAsFactors = F)
head(rt)

rt <- rt[,c("Gene", "padjust")]
colnames(rt) <- c("Gene",  "adjPValue_htseq_dexseq")

results[["htseq_dexseq"]] <- rt



rt <- read.table("Results_from_Katarina/dexseq_gene_featureCounts_version_1.10.8_results.txt", header = T, stringsAsFactors = F)
head(rt)

rt <- rt[,c("Gene", "padjust")]
colnames(rt) <- c("Gene",  "adjPValue_fc_dexseq")

results[["fc_dexseq"]] <- rt



rt <- read.table("Results_from_Katarina/dexseq_gene_Gordon_featureCounts_version_1.10.8_results.txt", header = T, stringsAsFactors = F)
head(rt)

rt <- rt[,c("Gene", "padjust")]
colnames(rt) <- c("Gene",  "adjPValue_fcG_dexseq")

results[["fcG_dexseq"]] <- rt



rt <- read.table("Results_from_Katarina/voom_gene_results.txt", header = T, stringsAsFactors = F)
head(rt)

rt <- rt[,c("gene.id","P.Value" ,"FDR")]
colnames(rt) <- c("Gene", "PValue_fc_voomex", "adjPValue_fc_voomex")

results[["fc_voomex"]] <- rt



rt <- read.table("Results_from_Katarina/voom_gene_Gordon_results.txt", header = T, stringsAsFactors = F)
head(rt)

rt <- rt[,c("gene.id","P.Value" ,"FDR")]
colnames(rt) <- c("Gene", "PValue_fcG_voomex", "adjPValue_fcG_voomex")

results[["fcG_voomex"]] <- rt



# # problems when adding this results (same genes have multiple p-values)
# # take min adj p-value 
# rt <- read.table("Results_from_Katarina/cuffdiff_results.txt", header = T, stringsAsFactors = F)
# head(rt)
# 
# rt <- rt[,c("ensembl_gene_id","p_value" ,"q_value")]
# colnames(rt) <- c("Gene", "PValue_cuffdiff", "adjPValue_cuffdiff")
# 
# o <- order(rt$adjPValue_cuffdiff, decreasing = FALSE)
# 
# rt <- rt[o, ]
# rt <- rt[!duplicated(rt$Gene, fromLast = FALSE), ]
# 
# results[["cuffdiff"]] <- rt




####################### DM v5 results / on htseq, fc, fcG

counting <- c("htseq", "fcG", "rsem" ,"fc")


for(j in counting){
  
  files <- list.files(path = paste0("DM_v5/", j, "/"), pattern = "_results.xls" )
  filtering  <- gsub(pattern = "_results.xls", replacement = "", x = files)
   
  for(i in seq_len(length(files))){
    
    rt <- read.table(paste0(paste0("DM_v5/", j, "/"), files[i]), header = T, stringsAsFactors = F)
    head(rt)
    rt <- rt[,c("GeneID","PValue" ,"FDR")]
    colnames(rt) <- c("Gene", paste0("PValue_",filtering[i]), paste0("adjPValue_",filtering[i]))
    
    results[[filtering[i]]] <- rt  
    
  }
  
}




#######################################################
# load information about simulation 
#######################################################

simu_info.e <- read.table("Simu_info/true_exons_simulation.txt", header=T)
simu_info.g <- read.table("Simu_info/true_genes_simulation.txt", header=T)

head(simu_info.g)
head(simu_info.e)

simu_info1 <- simu_info.g[!duplicated(simu_info.g$gene_id),c("gene_id", "status", "num", "rpk.mean")]

###################### count number of defferentially spliced exons

simu_info.e.spl <- split(simu_info.e, simu_info.e$gene_id)

num <- sapply(simu_info.e.spl, function(g){
  # g = simu_info.e.spl[[2]]
  
  n <- sum(g$status_exon==1)
  
  if(is.na(n))
    n <- 0
  
  return(n)
})

table(num, useNA = "always")

simu_info2 <- data.frame(gene_id = names(num), num.diff.ex=num)


simu_info <- merge(simu_info1, simu_info2, by = "gene_id", all=TRUE)
simu_info$num.diff.ex[is.na(simu_info$num.diff.ex)] <- 0

write.table(simu_info, "Simu_info/simu_info.xls", quote = F, sep = "\t", row.names = F, col.names = T)


#######################################################
# merge Diff Expr results with simu_info into one table
#######################################################


out.dir <- "PLOTS_DM_v5/"
dir.create(out.dir, showWarnings=F, recursive=T)

simu_info <- read.table("Simu_info/simu_info.xls", header=T)
dim(simu_info)

table <- simu_info

for(i in 1:length(results)){
  print(dim(results[[i]]))
  table <- merge(table, results[[i]], by.x="gene_id", by.y="Gene", all = TRUE)  
  table <- unique(table)
  print(dim(table))
}


write.table(table, paste0(out.dir,"/Table_all_results.xls"),  quote = F, sep = "\t", row.names = F, col.names = T)


## check table 

dim(table)
head(table)

dim(table[!complete.cases(table), ])
colnames(table)

    


##############################################################################################################
# ->>>>>> load results
##############################################################################################################

out.dir <- "PLOTS_DM_v5/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)


table <- read.table(paste0(out.dir,"/Table_all_results.xls"), header = T, stringsAsFactors = F)
tableNames <- colnames(table)

library(RColorBrewer)


allMethods <- gsub("adjPValue_", "",tableNames[grep(pattern = "adjPValue_", tableNames)])
allMethods

# rampNr <-  (length(allMethods) %/% 9 + 1) * 9 
# ramp <- colorRampPalette(brewer.pal(9,"Set1"))(rampNr)
# 
# set.seed(42)
# ord <- sample(rampNr, length(allMethods))
# ramp <- ramp[ord]


ramp <- colorRampPalette(brewer.pal(9,"Set1"))(length(allMethods))
# set.seed(42)
# ord <- sample.int(length(allMethods))
# ramp <- ramp[ord]

colors.org <- ramp
n.colors.org <- allMethods
names(colors.org) <- n.colors.org

# colors.org[c("htseq_dexseq", "htseq_dexseq_simes", "fc_voomex",  "cuffdiff", "dirmultDM_fc_g0_s4_keep0s", "fc_g0_s4_keep0s_subsetInf_DM4", "fc_g0_s4_keep0s_subsetInf_DM4adj")] <- c("pink", "magenta", "blue",  "red", "slateblue", "chartreuse3", "orange")


colors <- colors.org
n.colors <- n.colors.org

name <- ""



# pdf("Colors.pdf", width = 10, height = 5)
# 
# sets <- c("Set1", "Set2", "Set3",  "Accent", "Dark2", "Paired", "Pastel1", "Pastel2")
# nrs <- c(9, 8, 12, 8, 8, 12, 9, 8)
# 
# nr <- nrs[1]
# set <- sets[1]
# ramp <- colorRampPalette(brewer.pal(nr,set))(nr*2)
# plot(rep(1, length(ramp)), (1:length(ramp)) *4 , pch=19, cex=4, col=ramp, xlim=c(0,9), ylim=c(0, 100), xaxt="n")
# axis(side=1, at=1:8, labels=sets, las=1)
# 
# for(i in 2:length(sets)){
#   
#   nr <- nrs[i]
#   set <- sets[i]
#   ramp <- colorRampPalette(brewer.pal(nr,set))(nr*2)
#   points(rep(i, length(ramp)), (1:length(ramp)) *4 , pch=19, cex=4, col=ramp)
#    
# }
# 
# dev.off()


#######################################################
# histograms of p-values
#######################################################

colors <- colors.org
n.colors <- n.colors.org
colors <- colors.org[!grepl("dexseq",n.colors.org)]
n.colors <- n.colors.org[!grepl("dexseq",n.colors.org)]



pdf(paste0(out.dir, "/", name, "hist_pvalues.pdf"))

for(i in seq_len(length(colors))){
  # i=4
  cat(i, "\n")
  hist(table[, paste0("PValue_", n.colors[i]) ], col=colors[i], breaks=50, main=n.colors[i], cex.main=2, cex.lab=1.5, cex.axis = 1.5, xlab="P-values")
  
}

dev.off()



# pdf(paste0(out.dir, "/", name, "hist_pvalues-DEXSeqExon.pdf"))
# rt <- read.table("Results_from_Katarina/dexseq_exon_results.txt", header = T, stringsAsFactors = F)
# head(rt)
# hist(rt[, paste0("pvalue") ], col="pink", breaks=50, main="DEXSeq exon level", cex.main=2, cex.lab=1.5, cex.axis = 1.5, xlab="P-values")
# dev.off()





#######################################################
# generate ROCx plots 
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/ROCx.R")

colors <- colors.org
n.colors <- n.colors.org

n.colors

############ plots


methodsListLeg <- methodsList <- list(
  DMvoomex = n.colors[!grepl("dexseq", n.colors)],
  DMvoomex2 = n.colors[!grepl("dexseq", n.colors) & !grepl("fcG", n.colors)],
  DM = n.colors[!grepl("dexseq", n.colors) & !grepl("fc", n.colors)])



methodsList <- list(
  IMLS = c("fc_g0_s4_keep0s_subsetInf_DM5adj", "fc_g0_s4_keep0s_subsetInf_DM5TGcoadj"))

methodsListLeg <- list(
  IMLS = c("fc_DM_common", "fc_DM_tagwise"))


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
    
    pdf(paste0(out.dir, "/",name,"rocX_notNorm_",j,".pdf"))
  
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
  
  methodsLeg <- methodsListLeg[[j]]
  
  legend("bottomright", methodsLeg , col=colors[methods], lty=1, lwd=4, cex=1)
  dev.off()
  
  
}





#######################################################
# generate TPR vs achieved FDR plots
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/TPRvsFDR.R")

colors <- colors.org
n.colors <- n.colors.org

n.colors

############################## plots

methodsListLeg <- methodsList <- list(
  htseq = n.colors[grep("htseq", n.colors)],
  fcGfc =  n.colors[grep("fc", n.colors)] ,
  fchtseq = n.colors[!grepl("fcG", n.colors)],
  rsemhtseq = n.colors[!grepl("fc", n.colors) & grepl("adj", n.colors) & !grepl("M2", n.colors)])


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  pdf(paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf"))
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.4), ylim=c(0.5,1), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
  }
  
  methodsLeg <- methodsListLeg[[j]]
  
  legend("bottomright", methodsLeg, col=colors[methods], lty=1, lwd=4, cex=1, bty="o", bg="white", box.col="white")
  dev.off()
    
}

############################## plots

methodsListLeg <- methodsList <- list(
  fc = n.colors[grepl("fc_", n.colors) & !grepl("voomex", n.colors)])


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  pdf(paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf"))
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.3), ylim=c(0.7,0.9), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
  }
  
  methodsLeg <- methodsListLeg[[j]]
  
  legend("bottomright", methodsLeg, col=colors[methods], lty=1, lwd=4, cex=1, bty="o", bg="white", box.col="white")
  dev.off()
  
}

############################## plots

methodsListLeg <- methodsList <- list(
  DM = n.colors[grepl("fc_", n.colors) & grepl("DM", n.colors) & !grepl("M2", n.colors) & !grepl("M4", n.colors) & !grepl("M8", n.colors) & grepl("adj", n.colors)])


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  pdf(paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf"))
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.3), ylim=c(0.8,0.9), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
  }
  
  methodsLeg <- methodsListLeg[[j]]
  
  legend("bottomright", methodsLeg, col=colors[methods], lty=1, lwd=4, cex=1, bty="o", bg="white", box.col="white")
  dev.off()
  
}



############################## plots


methodsList <- list(
  dexseq = c("htseq_dexseq", "htseq_g0_s4_keep0s_subsetInf_DM5adj", "rsem_g0_s4_keep0s_subsetInf_DM5adj"))
methodsListLeg <- list(
  dexseq = c("htseq_dexseq", "htseq_DM", "rsem_DM"))


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  pdf(paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf"))
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.4), ylim=c(0.6,1), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
  }
  
  methodsLeg <- methodsListLeg[[j]]
  
  legend("bottomright", methodsLeg, col=colors[methods], lty=1, lwd=4, cex=1.5, bty="o", bg="white", box.col="white")
  dev.off()
  
}


############################## plots

methodsList <- list(
  IMLS = c("htseq_dexseq", "fc_g0_s4_keep0s_subsetInf_DM5adj", "fc_g0_s4_keep0s_subsetInf_DM5TGcoadj", "rsem_g0_s4_keep0s_subsetInf_DM5adj"))

methodsListLeg <- list(
  IMLS = c("htseq_dexseq", "fc_DM_common", "fc_DM_tagwise", "rsem_DM_common"))



for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  pdf(paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf"))
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.3), ylim=c(0.7,1), lwd=7, cex=3.3, cex.lab=1.45, cex.axis=1.5)
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=7, cex=3.3)
  }
  
  methodsLeg <- methodsListLeg[[j]]
  
  legend("bottomright", methodsLeg, col=colors[methods], lty=1, lwd=4, cex=1.5, bty="o", bg="white", box.col="white")
  dev.off()
  
}



#######################################################
# venn diagrams --- my function
#######################################################
library(VennDiagram)

FDR.cutoff <- 0.05


plotVenn <- function(venne.genes, colors, venn.methods, methodsLeg, name1="", name2="", cex=1, cat.cex=1.4, lwd=2, lty=1, alpha=0.5, margin=0.05){
  
  colors.col <- colors[venn.methods]
  
  if(length(venn.methods)==4)
    colors.col <- colors.col[c(1, 3, 4 ,2)]
  
  venn.d <- venn.diagram(venne.genes[venn.methods], filename=NULL, fill = colors[venn.methods], col=colors.col, cex=cex, cat.cex=cat.cex, lwd=lwd, lty=lty, alpha=alpha, margin=margin, category.names=methodsLeg)
  
  pdf(paste0(out.dir, "/", name1, "Venn_Diagr_", name2,".pdf"))
  grid.draw(venn.d)
  dev.off()
  
  
}


#######################################################
# generate venn diagrams 
#######################################################


## TP as a separate circle
venne.genes <- list()
for(i in n.colors){ 
  venne.genes[[i]] <- na.omit(table[ table[,paste0("adjPValue_", i)] < FDR.cutoff, "gene_id"]) 
}
venne.genes$True <- table[table$status == 1, "gene_id"]



methodsListLeg <- methodsList <- list(
  DM = n.colors[grepl("fc_", n.colors) & grepl("DM", n.colors) & !grepl("M2", n.colors) & !grepl("M4", n.colors) & !grepl("M8", n.colors) & !grepl("adj", n.colors)],
  DMadj = n.colors[grepl("fc_", n.colors) & grepl("DM", n.colors) & !grepl("M2", n.colors) & !grepl("M4", n.colors) & !grepl("M8", n.colors) & grepl("adj", n.colors)])


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  methodsLeg <- methodsListLeg[[j]]
  
  plotVenn(venne.genes, colors=c(colors[methods], True="grey"), venn.methods = c(methods, "True"), methodsLeg = c(methodsLeg, "TRUE"),  name1=name, name2=j, margin=0.1, cat.cex=0.8, cex=1.7)
  
}












































































