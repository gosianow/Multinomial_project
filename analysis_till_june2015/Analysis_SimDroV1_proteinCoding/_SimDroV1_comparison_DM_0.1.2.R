#######################################################

# Created 3 Feb 2014 
# BioC 3.0

# comprare the DM_0.1.2: tagwise Dipsersion and filtering; using htseq counts and bitseq counts
# For DM: Filtering005, Filtering0, Filtering_DEXSeq

# Update 5 Feb 2014:

# Update 4 Mar 2015 
# - add dispersion vs mean plot with marked FP
# - add plot of raw and DM estimated ratios for DM FP


#######################################################


setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding")

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata





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
### remove duplicated genes
rt <- rt[order(rt$FDR, decreasing = FALSE),]
rt <- rt[!duplicated(rt$gene_id), ]

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

 Filtering <- "Filtering005/"
# Filtering <- "Filtering0/"
# Filtering <- "Filtering_DEXSeq/"


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


out.dir <- "PLOTS_DM_0.1.2_Filtering005/"
# out.dir <- "PLOTS_DM_0.1.2_Filtering0/"
# out.dir <- "PLOTS_DM_0.1.2_Filtering_DEXSeq/"
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



tdp <- table[duplicated(table$gene_id, fromLast = T) | duplicated(table$gene_id, fromLast = F),]
dim(tdp)


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



### load table with all results

table <- read.table(paste0(out.dir,"/Table_all_results.xls"), header = T, stringsAsFactors = F)
tableNames <- colnames(table)


### assign colors to different methods

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
  
  
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.3), ylim=c(0.4,1), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
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
  
  
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.3), ylim=c(0.9,1), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
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
  
  
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.3), ylim=c(0.6,1), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
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

  out.name <- gsub(pattern = "\\.", replacement = "_" ,paste0( name1, "Venn_Diagr_", name2))
  
  pdf(paste0(out.dir, "/", out.name, ".pdf"))
  grid.draw(venn.d)
  dev.off()
  
  
}


#######################################################
# generate venn diagrams 
#######################################################

colors <- colors.org
n.colors <- n.colors.org

n.colors


## TP as a separate circle
venne.genes <- list()
for(i in n.colors){ 
  venne.genes[[i]] <- na.omit(table[ table[,paste0("adjPValue_", i)] < FDR.cutoff, "gene_id"]) 
}
venne.genes$True <- na.omit(table[table$status == 1, "gene_id"])

names(venne.genes)


count.method <- "htseq"
count.method <- "bitseq"

DM.methods <- n.colors[grepl(count.method, n.colors) & grepl("DM", n.colors)]

methodsListLeg <- methodsList <- list()

for(i in 1:length(DM.methods)){
  methodsListLeg[[DM.methods[i]]] <- methodsList[[DM.methods[i]]] <- c(DM.methods[i], paste0(count.method, "_dexseq"))
}




for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  methodsLeg <- methodsListLeg[[j]]
  
  plotVenn(venne.genes, colors=c(colors[methods], True="grey"), venn.methods = c(methods, "True"), methodsLeg = c(methodsLeg, "TRUE"),  name1=name, name2=j, margin=0.1, cat.cex=2, cex=1.7)
  
}



##############################################################################################################
# Plots of FP called by DM
##############################################################################################################

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)


#######################################################
# plot dispersion vs mean & mark FP
#######################################################

FilteringList <- c("Filtering_DEXSeq", "Filtering0", "Filtering005")
count.methodList <- c("bitseq", "htseq")
count.method.fileList <- c("bitseq_1.10.0", "htseq")

for(f in 1:length(FilteringList)){
  
  Filtering <- FilteringList[f]
  
  out.dir <- file.path(paste0("PLOTS_DM_0.1.2_", Filtering))
  dir.create(out.dir, showWarnings=F, recursive=T)
  
  ### load table with all results
  
  table <- read.table(paste0(out.dir,"/Table_all_results.xls"), header = T, stringsAsFactors = F)
  tableNames <- colnames(table)
   
  ### assign colors to different methods
  
  allMethods <- gsub("adjPValue_", "",tableNames[grep(pattern = "adjPValue_", tableNames)])
  
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
  
  
  ## TP as a separate circle
  venne.genes <- list()
  for(i in allMethods){ 
    venne.genes[[i]] <- na.omit(table[ table[,paste0("adjPValue_", i)] < FDR.cutoff, "gene_id"]) 
  }
  venne.genes$True <- na.omit(table[table$status == 1, "gene_id"])
  
  
  for(cm in 1:length(count.method.fileList)){
    
    count.method.file <- count.method.fileList[cm]
    count.method <- count.methodList[cm]
       
    res.path <- file.path("DM_0.1.2", count.method.file, Filtering)
      
    DM.methods <- allMethods[grepl(count.method, allMethods) & grepl("DM.*TG", allMethods)]
    print(DM.methods)
    
    for(m in 1:length(DM.methods)){
      
      method <- DM.methods[m]
      method
      method.ref <- paste0(count.method, "_dexseq")
      method.ref
      
      ### list of FP to plot
#       FP2plot <- setdiff(setdiff(venne.genes[[method]], venne.genes[["True"]]), venne.genes[[method.ref]])
      FP2plot <- setdiff(venne.genes[[method]], venne.genes[["True"]])
      length(FP2plot)
      
      ### load DM pipeline results
      dgeDM.rdata <- file.path(res.path, paste0(gsub( "\\.", "-",method), "_dgeDM.RData"))
      print(dgeDM.rdata)
      load(dgeDM.rdata)
      dgeDM
      
      commonDispersion <- dgeDM$commonDispersion
      tagwiseDispersion <- dgeDM$tagwiseDispersion
      meanExpr <- dgeDM$meanExpr
    
      
      pdf(file.path(out.dir, paste0("TREMDgammaVSmean_", gsub( "\\.", "_",method),".pdf")))
      
      # smoothScatter(log10(meanExpr), log10(tagwiseDispersion), xlab="log10 mean gene expression", ylab="log10 gamma +", nrpoints = Inf, colramp=colorRampPalette(c("white", "darkgrey")), pch = 19, cex=0.6, las = 2)
      # points(log10(meanExpr[FP2plot]), log10(tagwiseDispersion[FP2plot]), pch = 19, cex=1, col= colors[method])
      # abline(h = log10(commonDispersion), col = "grey", lwd = 3, lty = 3)
      
      df <- data.frame(meanExpr = log10(meanExpr+1), tagwiseDispersion = log10(tagwiseDispersion))
      ggp <- ggplot(df, aes(x = meanExpr, y = tagwiseDispersion)) +
        theme_bw() +
        labs(x = "log10 mean expression", y = "log10 tagwise gamma+") +
        theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold")) +
        geom_point(size = 2) +
        geom_point(data = df[names(meanExpr) %in% FP2plot, ], aes(x = meanExpr, y = tagwiseDispersion), colour = colors[method], size = 3) +
        geom_hline(aes(yintercept=log10(commonDispersion)), colour = "grey", linetype="dashed", size = 1)
      print(ggp)
      
      
      dev.off()
      
      
    }
    
  }
  
}



#######################################################
# plot expression ratios of FP called by DM
#######################################################

nr.topG <- 20

FilteringList <- c("Filtering_DEXSeq", "Filtering0", "Filtering005")
count.methodList <- c("bitseq", "htseq")
count.method.fileList <- c("bitseq_1.10.0", "htseq")

for(f in 1:length(FilteringList)){
  # f = 1
  Filtering <- FilteringList[f]
  
  out.dir <- file.path(paste0("PLOTS_DM_0.1.2_", Filtering))
  dir.create(out.dir, showWarnings=F, recursive=T)
  
  ### load table with all results
  
  table <- read.table(paste0(out.dir,"/Table_all_results.xls"), header = T, stringsAsFactors = F)
  tableNames <- colnames(table)
  rownames(table) <- table$gene_id
  allMethods <- gsub("adjPValue_", "",tableNames[grep(pattern = "adjPValue_", tableNames)])
  

  ## TP as a separate circle
  venne.genes <- list()
  for(i in allMethods){ 
    venne.genes[[i]] <- na.omit(table[ table[,paste0("adjPValue_", i)] < FDR.cutoff, "gene_id"]) 
  }
  venne.genes$True <- na.omit(table[table$status == 1, "gene_id"])
  
  
  for(cm in 1:length(count.method.fileList)){
    # cm = 1
    count.method.file <- count.method.fileList[cm]
    count.method <- count.methodList[cm]
    
    res.path <- file.path("DM_0.1.2", count.method.file, Filtering)
    
    DM.methods <- allMethods[grepl(count.method, allMethods) & grepl("DM.*TG", allMethods)]
    print(DM.methods)
    
    for(m in 1:length(DM.methods)){
      
      method <- DM.methods[m]
      method
      method.ref <- paste0(count.method, "_dexseq")
      method.ref
      
      ### list of FP to plot
        FP2plot <- setdiff(setdiff(venne.genes[[method]], venne.genes[["True"]]), venne.genes[[method.ref]])
#       FP2plot <- setdiff(venne.genes[[method]], venne.genes[["True"]])
      print(length(FP2plot))
      
      ### order FP genes by significance
      FP2plot <- FP2plot[order(table[table$gene_id %in% FP2plot, paste0("adjPValue_", method)], decreasing = FALSE)]

      ### load DM pipeline results
      dgeDM.rdata <- file.path(res.path, paste0(gsub( "\\.", "-",method), "_dgeDM.RData"))
      print(dgeDM.rdata)
      load(dgeDM.rdata)
      # dgeDM
      
      commonDispersion <- dgeDM$commonDispersion
      tagwiseDispersion <- dgeDM$tagwiseDispersion
      meanExpr <- dgeDM$meanExpr
      
      
      pdf(file.path(out.dir, paste0("ExonProportions_topFP_", gsub( "\\.", "_",method),".pdf")), width = 7, height = 3)
      
      for(g in 1:min(nr.topG, length(FP2plot))){
        # g = 1
        gene <- FP2plot[g]
        # print(gene)
        expr <- dgeDM$counts[[gene]]
        colnames(expr) <- metadata$SampleName
        rownames(expr) <- subset(dgeDM$genes, gene_id==gene)$ete_id    
        tot <- colSums(expr)
        labels <- strsplit2(rownames(expr), ":")[,2]
        prop.smp <- data.frame( ete_id =  labels, t(apply(expr, 1, function(t){ t / tot })))  
        n <- nrow(expr)  
        prop.est <- data.frame(ete_id = labels, dgeDM$fit[[gene]]$piH)
        
        prop.smp.m <- melt(prop.smp, id.vars = "ete_id", variable.name = "Samples", value.name = "Proportions")  
        Condition <- prop.smp.m$Samples 
        levels(Condition) <- substr(levels(Condition), 1,2)
        
        prop.est.m <- melt(prop.est, id.vars = "ete_id", variable.name = "Samples", value.name = "Proportions")
        
#         ### line plots
#         ggp <- ggplot() +
#           theme_bw() +
#           geom_line(data = prop.smp.m, aes(x = ete_id, y = Proportions, group = factor(Samples), colour = Condition )) +
#           geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, group = factor(Samples), colour = Samples ), size = 3) +
#           theme(axis.text.x = element_text(angle = 30, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=14, face="bold"), legend.position="none") +
#           ggtitle(paste0(gene))
#         
#         print(ggp)
        
        ### bar plots
        ggb <- ggplot(prop.smp.m, aes(x = ete_id, y = Proportions, fill = factor(Samples))) +
          theme_bw() + 
          theme(axis.text.x = element_text(angle = 20, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=14, face="bold"), legend.position="none") +
          ggtitle(paste0(gene, " FDR = ", sprintf("%.02e",table[gene, paste0("adjPValue_", method)]), "\n log10 gamma = ", round(log10(tagwiseDispersion[gene]), digits = 2), " log10 mean expr = ", round(log10(meanExpr[gene] +1), 2))) +
          geom_bar(stat = "identity", position=position_dodge()) +       
          scale_fill_manual(values=c(rep("firebrick1", 4), rep("dodgerblue", 4))) +
          geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, colour = factor(Samples), fill = factor(Samples)), position = position_jitterdodge(jitter.width = 0, jitter.height = 0 ), shape = 18, size = 3) +
          scale_colour_manual(values=c(rep("firebrick", 1), rep("dodgerblue4", 1))) +
          coord_cartesian(ylim = c(0, 1)) 

        print(ggb)
        
        
        
      }
      
      
      dev.off()
      
      
    }
    
  }
  
}











































