#######################################################
# 
# Created 3 Feb 2014 

# comprare the DM_0.1.2: tagwise Dipsersion and filtering; using htseq counts and bitseq counts
# For DM: Filtering005

# Update 5 Feb 2014:

# BioC 3.0

#######################################################


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding")

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata




#######################################################
# calculate gene-level p-values for DEXSeq with Simes' method
#######################################################


DEXSeqSimes <- function(dexseqExonFile, dexseqGeneFile, dexseqGeneFileSimes){
  
  dexseq.ex <- read.table(dexseqExonFile, header = T, as.is = TRUE, sep = "\t", fill = TRUE)
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


dexseqExonFile <- "Results_from_Katarina/dexseq_1.10.8_exon_htseq_results.txt"
dexseqGeneFile <- "Results_from_Katarina/dexseq_1.10.8_gene_htseq_results.txt"
dexseqGeneFileSimes <- "Results_from_Katarina/dexseq_1.10.8_gene_Simes_htseq_results.txt"


DEXSeqSimes(dexseqExonFile, dexseqGeneFile, dexseqGeneFileSimes)




#######################################################
# merge Diff Expr results 
#######################################################

results <- list()

####################### results produced by Katarina

### DEXSeq
rt <- read.table("Results_from_Katarina/dexseq_1.10.8_gene_htseq_results.txt", header = T, as.is = TRUE)
head(rt)

rt <- rt[,c("Gene", "padjust")]
colnames(rt) <- c("Gene",  "adjPValue_htseq_dexseq")

results[["htseq_dexseq"]] <- rt

sum(results[["htseq_dexseq"]]$adjPValue_htseq_dexseq < 0.05)

### Cuffdiff

rt <- read.table("Results_from_Katarina/cuffdiff_version_2.2.1_results.txt", header = T, as.is = TRUE)
head(rt)

rt <- rt[,c("gene_id", "FDR")]
colnames(rt) <- c("Gene",  "adjPValue_cuffdiff")

results[["cuffdiff"]] <- rt

####################### Other results 

### DEXSeq on bitseq counts

rt <- read.table("DEXSeq_1.10.8/bitseq_1.10.0/dexseq_1.10.8_gene_bitseq_results.txt", header = T, as.is = TRUE)
head(rt)

rt <- rt[,c("Gene", "padjust")]
colnames(rt) <- c("Gene",  "adjPValue_bitseq_dexseq")

results[["bitseq_dexseq"]] <- rt




####################### DM_0.1.2 results on htseq counts and bitseq 

count.method <- c("htseq", "bitseq_1.10.0")

# Filtering <- "Filtering005/"
# Filtering <- "Filtering0/"
Filtering <- "Filtering_DEXSeq/"


for(j in count.method){
  # j = count.method[2]
  
  out.path <- paste0("DM_0.1.2/", j, "/", Filtering)
  
  files <- list.files(path = out.path, pattern = "_results.xls" )
  filtering  <- gsub(pattern = "_results.xls", replacement = "", x = files)
   
  for(i in seq_len(length(files))){
    
    rt <- read.table(paste0(out.path, files[i]), header = T, stringsAsFactors = F)
    head(rt)
    rt <- rt[,c("GeneID","PValue" ,"FDR")]
    colnames(rt) <- c("Gene", paste0("PValue_",filtering[i]), paste0("adjPValue_",filtering[i]))
    
    results[[filtering[i]]] <- rt  
    
  }
  
}

names(results)


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


# out.dir <- "PLOTS_DM_0.1.2_Filtering005/"
# out.dir <- "PLOTS_DM_0.1.2_Filtering0/"
out.dir <- "PLOTS_DM_0.1.2_Filtering_DEXSeq/"
dir.create(out.dir, showWarnings=F, recursive=T)



simu_info <- read.table("Simu_info/simu_info.xls", header=T)
dim(simu_info)

table <- simu_info

for(i in 1:length(results)){
  print(dim(results[[i]]))
  table <- merge(table, results[[i]], by.x="gene_id", by.y="Gene", all.x = TRUE)  
  table <- unique(table)
  print(dim(table))
}


write.table(table, paste0(out.dir,"/Table_all_results.xls"),  quote = F, sep = "\t", row.names = F, col.names = T)


## check table 

dim(table)
head(table)
colnames(table)

    


##############################################################################################################
# ->>>>>> load results
##############################################################################################################

library(ggplot2)
library(reshape2)
library(gridExtra)


# out.dir <- "PLOTS_DM_0.1.2_Filtering005/"
# out.dir <- "PLOTS_DM_0.1.2_Filtering0/"
out.dir <- "PLOTS_DM_0.1.2_Filtering_DEXSeq/"

dir.create(out.dir, showWarnings=F, recursive=T)



table <- read.table(paste0(out.dir,"/Table_all_results.xls"), header = T, stringsAsFactors = F)
tableNames <- colnames(table)


allMethods <- gsub("adjPValue_", "",tableNames[grep(pattern = "adjPValue_", tableNames)])
allMethods

### generate ggplot hue colors 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

ramp <- gg_color_hue(length(allMethods))


colors.org <- ramp
n.colors.org <- allMethods
names(colors.org) <- n.colors.org

colors <- colors.org
n.colors <- n.colors.org

name <- ""






#######################################################
# histograms of p-values
#######################################################

colors <- colors.org
n.colors <- n.colors.org
colors <- colors.org[grepl("DM",n.colors.org)]
n.colors <- n.colors.org[grepl("DM",n.colors.org)]



pdf(paste0(out.dir, "/", name, "hist_pvalues.pdf"), 5, 5)

for(i in seq_len(length(colors))){
  # i=1
  cat(i, "\n")
  
  df <- data.frame(PValues = table[, paste0("PValue_", n.colors[i]) ])
  
  ggp <- ggplot(df, aes(x = PValues)) +
    theme_bw() +
    ggtitle(n.colors[i]) +
    geom_histogram(binwidth = 0.02, colour=colors[i], fill=colors[i]) +
    theme(text = element_text(size=15))
   
  print(ggp)
  
}

dev.off()



#######################################################
# generate ROCx plots 
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/ROCx.R")

colors <- colors.org
n.colors <- n.colors.org

n.colors




methodsListLeg <- methodsList <- list(
  DM = n.colors[grepl("DM", n.colors)], 
  bitseq = n.colors[grepl("bitseq", n.colors) & !grepl("dexseq", n.colors)])


for(j in names(methodsList)){
  # j = names(methodsList)[2]
  
  methods <- methodsList[[j]]
  
  scores <- list()
  
  for(i in methods)
  scores[[i]] <- data.frame(pvs = table[,paste0("PValue_", i)], apvs = table[ ,paste0("adjPValue_", i)], status = table$status)
  
  names(scores) <- methods
  
  
  pdf(paste0(out.dir, "/",name,"rocX_notNorm_",j,".pdf"), width = 7, height = 7)
  
  ggROCx_notNorm(scores = scores, colours = colors[methods])

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
  All = n.colors)


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  pdf(paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf"))
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  
  NAs <- is.na(status)
  apvs <- apvs[!NAs]
  status <- status[!NAs]
  
  
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.5), ylim=c(0.4,1), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    NAs <- is.na(status)
    apvs <- apvs[!NAs]
    status <- status[!NAs]
    
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
  }
  
  methodsLeg <- methodsListLeg[[j]]
  
  legend("bottomright", methodsLeg, col=colors[methods], lty=1, lwd=4, cex=1, bty="o", bg="white", box.col="white")
  dev.off()
    
}


############################## plots

methodsListLeg <- methodsList <- list(
  bitseq = n.colors[grepl("bitseq", n.colors)])


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  colors <- gg_color_hue(length(methods))
  names(colors) <- methods
  
  pdf(paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf"))
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  
  NAs <- is.na(status)
  apvs <- apvs[!NAs]
  status <- status[!NAs]
  
  
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.4), ylim=c(0.8,1), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    NAs <- is.na(status)
    apvs <- apvs[!NAs]
    status <- status[!NAs]
    
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
  }
  
  methodsLeg <- methodsListLeg[[j]]
  
  legend("bottomright", methodsLeg, col=colors[methods], lty=1, lwd=4, cex=1, bty="o", bg="white", box.col="white")
  dev.off()
  
}



############################## plots

methodsListLeg <- methodsList <- list(
  htseq = n.colors[grepl("htseq", n.colors)])


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  colors <- gg_color_hue(length(methods))
  names(colors) <- methods
  
  pdf(paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf"))
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  
  NAs <- is.na(status)
  apvs <- apvs[!NAs]
  status <- status[!NAs]
  
  
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.3), ylim=c(0.8,1), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    NAs <- is.na(status)
    apvs <- apvs[!NAs]
    status <- status[!NAs]
    
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
  }
  
  methodsLeg <- methodsListLeg[[j]]
  
  legend("bottomright", methodsLeg, col=colors[methods], lty=1, lwd=4, cex=1, bty="o", bg="white", box.col="white")
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












































































