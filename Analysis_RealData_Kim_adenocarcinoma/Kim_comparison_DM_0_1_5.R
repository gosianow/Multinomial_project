######################################################
# BioC 3.1

# Created 22 June 2015 
# Comparison of DM_0.1.5 with DEXSeq on Kim_adenocarcinoma data using htseq and kallisto counts 

######################################################


setwd("/home/Shared/data/seq/Kim_adenocarcinoma/")

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)

DM_out <- "DM_0_1_5/"
DM_plots_out <- "DM_0_1_5_PLOTS/"



##########################################################################
# load metadata
##########################################################################

metadata <- read.table("3_metadata/Malgorzata Nowicka2014-11-04GSE37764.csv", stringsAsFactors=F, sep=",", header=T) 
metadata <- metadata[metadata$X == "RNA-seq",]

metadata$sampleName <- metadata$ids
metadata$condition <- metadata$Tissue.Type

metadata <- metadata[order(metadata$condition), ]

metadata


#######################################################
# merge Diff Expr results 
#######################################################

# model <- "diff_out_Full_tumorVSnormal"
# 
# count.methodList <- c("htseq", "kallisto")
# count.method <- count.methodList[1]
# 
# filter.methodList <- c("Filtering_DEXSeq", "Filtering_min6prop0_01min9cpm1maxNA")  
# filter.method <- filter.methodList[1]



results <- list()

####################### results produced by Charlotte


rt <- read.table(paste0("4_results/DEXSeq_1.10.8/", model,"/DEXSeq_", count.method, "_gene_results.txt"), header = T, as.is = TRUE)
head(rt)

colnames(rt) <- c("gene_id",  paste0("adjPValue.",count.method,"_dexseq"))

head(rt)
results[[paste0( count.method, "_dexseq" )]] <- rt




####################### DM_0.1.5 results on htseq counts and bitseq 



res.path <- paste0(DM_out, model, "/",count.method, "/", filter.method, "/")

files <- list.files(path = res.path, pattern = "_results.xls" )
res  <- gsub(pattern = "_results.xls", replacement = "", x = files)
res 

for(i in 1:length(files)){
  # i = 1
  
  rt <- read.table(paste0(res.path, files[i]), header = TRUE, as.is = TRUE)
  head(rt)
  
  rt <- rt[,c("geneID", "pValue" ,"FDR")]
  
  colnames(rt) <- c("gene_id", paste0("PValue.",res[i]), paste0("adjPValue.",res[i]))
  head(rt)
  
  results[[res[i]]] <- rt  
  
}


names(results)

results <- Reduce(function(...) merge(..., by = "gene_id", all=TRUE, sort = FALSE), results)




#######################################################
# merge Diff Expr results into one table
#######################################################


plots.path <- paste0(DM_plots_out, model, "/", count.method, "/", filter.method, "/")
dir.create(plots.path, showWarnings=F, recursive=T)

table <- results
write.table(table, paste0(plots.path,"/table_all_results.xls"),  quote = F, sep = "\t", row.names = F, col.names = T)


##############################################################################################################
# ->>>>>> load results
##############################################################################################################

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


# model <- "diff_out_Full_tumorVSnormal"
# 
# count.methodList <- c("htseq", "kallisto")
# count.method <- count.methodList[1]
# 
# filter.methodList <-  c("Filtering_DEXSeq", "Filtering_min6prop0_01min9cpm1maxNA")  
# filter.method <- filter.methodList[1]


plots.path <- paste0(DM_plots_out, model, "/", count.method, "/", filter.method, "/")

name <- ""




### load table with all results

table <- read.table(paste0(plots.path,"/table_all_results.xls"), header = T, stringsAsFactors = F)
tableNames <- colnames(table)

### assign colors to different methods

allMethods <- data.frame(names = gsub("adjPValue.", "",tableNames[grep(pattern = "adjPValue.", tableNames)]), stringsAsFactors = FALSE)
rownames(allMethods) <- allMethods$names
allMethods


### generate ggplot hue colors 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}


colorb <- function(n) {
  colorRampPalette(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))(n)
}


allMethods <- allMethods[order(allMethods$names), , drop = FALSE]

ramp <- colorb(nrow(allMethods))

allMethods$colors <- ramp

allMethods



#######################################################
# histograms of p-values
#######################################################


colors <- allMethods[grepl("DM",allMethods$names), ]

pdf(paste0(plots.path, "/", name, "hist_pvalues.pdf"), 5, 5)

for(i in seq_len(nrow(colors))){
  # i=1
  cat(i, "\n")
  
  PValues = table[, paste0("PValue.", colors$names[i]) ] 
  hist(PValues, breaks = 100, main = paste0(colors$names[i]), col = colors$colors[i], cex.main = 1)
  
}

dev.off()



#######################################################
# generate venn diagrams 
#######################################################


source("/home/gosia/R/R_Multinomial_project/Plot_Functions/plotVenn_0_1_4.R")


colors <- allMethods

FDR.cutoff <- 0.05

## TP as a separate circle
venn.list <- list()
for(i in colors$names){ 
  venn.list[[i]] <- na.omit(table[ table[,paste0("adjPValue.", i)] < FDR.cutoff, "gene_id"]) 
}
names(venn.list)


colors <- allMethods[grepl("DM", allMethods$names), ]


methodsListLeg <- methodsList <- list()
for(i in colors$names){
  methodsList[[i]] <- c(paste0(count.method, "_dexseq"), i)
  methodsListLeg[[i]] <- c(paste0(count.method, "_dexseq"), paste0(count.method, "_DM"))
}


colors <- allMethods

for(j in names(methodsList)){
  # j = names(methodsList)[1]
  methods <- methodsList[[j]]
  methodsLeg <- methodsListLeg[[j]]
  
  plotPath =  paste0(plots.path, "/", name, "Venn_Diagram_", gsub( "\\.", "_",j) ,".pdf")
  
  venn.colors <- colors[methods, "colors"]
  names(venn.colors) <- methods
  
  venn.colorsLeg <- venn.colors
  names(venn.colorsLeg) <- methodsLeg
  venn.listLeg <- venn.list[methods]
  names(venn.listLeg) <- methodsLeg
  
  plotVenn2(venn.listLeg, venn.colors = venn.colorsLeg, venn.subset = c(methodsLeg), margin = 0, cat.cex=1.2, cex=1.7, plotPath = plotPath)
  
}




##################################################################################
################################# plot proportions
##################################################################################
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


plotProportions <- function(dgeDM, genes2plot, plotPath){
  
  rownames(dgeDM$table) <- dgeDM$table$geneID
  
  pdf(plotPath, width = 12, height = 7)
  
  
  for(g in 1:length(genes2plot)){
    # g = 1
    
    gene <- genes2plot[g]
    # print(gene)
    Condition <- dgeDM$samples$group 
    expr <- dgeDM$counts[[gene]]
    
    labels <- strsplit2(rownames(expr), ":")[, 2]
    
    prop.smp <- data.frame( ete_id =  labels, prop.table(expr, 2))  
    prop.est <- data.frame(ete_id = labels, dgeDM$fitFull[[gene]]$piH)
    
    prop.est.null <- data.frame(ete_id = labels, dgeDM$fitNull[[gene]]$piH)
    
    #### order transcipts by decreasing proportions 
    order.tr <- labels[order(apply(aggregate(t(prop.smp[, -1]), by = list(Condition = Condition), median)[, -1], 2, max), decreasing = TRUE)]   
    
    prop.smp.m <- melt(prop.smp, id.vars = "ete_id", variable.name = "Samples", value.name = "Proportions")  
    prop.smp.m$ete_id <- factor(prop.smp.m$ete_id, levels = order.tr)
    prop.smp.m$Samples <- factor(prop.smp.m$Samples)
    prop.smp.m$Condition <- rep(Condition, each = nrow(prop.smp))
    
    prop.est.m <- melt(prop.est, id.vars = "ete_id", variable.name = "Condition", value.name = "Proportions")
    prop.est.m$ete_id <- factor(prop.est.m$ete_id, levels = order.tr)
    prop.est.m$Condition <- factor(prop.est.m$Condition)
    prop.est.m$Samples <- factor(prop.est.m$Condition)
    
    colnames(prop.est.null) <- c("ete_id", "Proportions")
    prop.est.null$ete_id <- factor(prop.est.null$ete_id, levels = order.tr)
    prop.est.null$Samples <- "Null"
    
    main <- paste0(gene, "\n Mean Expression = ", ifelse(is.null(dgeDM$meanExpr), NA, round(dgeDM$meanExpr[gene])), " / Dispersion = ", round(dgeDM$fitNull[[gene]]$gamma0, 2), "\n LR = ", round(dgeDM$table[gene, "LR"], 4) , " / P-value = ", sprintf("%.02e", dgeDM$table[gene, "pValue"]), " / FDR = ", sprintf("%.02e", dgeDM$table[gene, "FDR"]))
    
    
    #     ### box plots with points2
    #     ggb <- ggplot(prop.smp.m, aes(x = ete_id, y = Proportions)) +
    #       theme_bw() + 
    #       theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), legend.position="none", plot.title = element_text(size=10)) +
    #       ggtitle(main) +     
    #       geom_jitter(aes(fill = Condition, colour = factor(Condition)), position = position_jitterdodge(dodge.width = 0.75), alpha = 0.5) +
    #       geom_boxplot(aes(colour = Condition), fill = "white", outlier.size = NA, alpha = 0) + 
    #       geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, fill = Condition), position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 3, shape = 23) +
    #       geom_point(data = prop.est.null, aes(x = ete_id, y = Proportions), size = 3, shape = 23, fill = "orange") +
    #       coord_cartesian(ylim = c(-0.1, 1.1)) 
    #     
    #     
    #     print(ggb)
    #     
    #         
    #     ### line plots
    #     ggp <- ggplot() +
    #       theme_bw() +
    #       geom_line(data = prop.smp.m, aes(x = ete_id, y = Proportions, group = Samples, colour = Condition )) +
    #       geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, group = Condition, fill = Condition ), size = 3, shape = 23) +
    #       theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), legend.position="none", plot.title = element_text(size=10)) +
    #       ggtitle(main)+
    #       coord_cartesian(ylim = c(-0.1, 1.1))
    #     
    #     print(ggp)
    
    
    
    ### bar plots
    ggb <- ggplot(prop.smp.m, aes(x = ete_id, y = Proportions, group = Samples, fill = Condition)) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=12, face="bold"), legend.position="none", plot.title = element_text(size=10)) +
      ggtitle(main) +
      geom_bar(stat = "identity", position=position_dodge()) +       
      geom_point(data = prop.est.m, aes(x = ete_id, y = Proportions, fill = Condition), position = position_jitterdodge(jitter.width = 0, jitter.height = 0 ), size = 3, shape = 23) +
      geom_point(data = prop.est.null, aes(x = ete_id, y = Proportions), size = 3, shape = 23, fill = "orange") +
      coord_cartesian(ylim = c(-0.1, 1.1))
    
    print(ggb)
    
    
  }
  
  
  dev.off()
  
  
  
  
}




#######################################################
# plot expression ratios of FP and FN called by DM
#######################################################


out.dir.prop <- paste0(plots.path, "/Proportions/")
dir.create(out.dir.prop, showWarnings=F, recursive=T)

FDR.cutoff <- 0.05


### generate list of significant genes for each method
venne.genes <- list()
for(i in allMethods$names){ 
  venne.genes[[i]] <- na.omit(table[ table[, paste0("adjPValue.", i)] < FDR.cutoff, "gene_id"]) 
}
names(venne.genes)


res_path <- paste0(DM_out, model, "/" ,count.method, "/", filter.method, "/")

DM.methods <- allMethods[grepl("DM", allMethods$names), "names"]


for(m in 1:length(DM.methods)){
  # m = 1
  
  method <- DM.methods[m]
  method
  method.ref <- paste0(count.method, "_dexseq")
  method.ref
  
  ### list of FP
  FP <- setdiff(venne.genes[[method]], venne.genes[[method.ref]])
  # FP2plot <- setdiff(venne.genes[[method]], venne.genes[["True"]])
  cat("Number of FP: ",length(FP), "\n")
  
  ### order FP genes by significance
  table.FP <- table[table$gene_id %in% FP, c("gene_id", paste0("PValue.", method), paste0("adjPValue.", method))]
  table.FP <- table.FP[order(table.FP[, 2], decreasing = FALSE), ]
  
  genes2plot <- table.FP[1:min(20, nrow(table.FP)), "gene_id"]
  
  ### load DM pipeline results
  dgeDM.rdata <- file.path(res_path, paste0(gsub( "\\.", "-",method), "_dgeDM.RData"))
  print(dgeDM.rdata)
  load(dgeDM.rdata)
  # dgeDM
  
  
  plotPath <- file.path(out.dir.prop, paste0("ExonProportions_topFP_", gsub( "\\.", "_",method),".pdf"))
  
  plotProportions(dgeDM, genes2plot, plotPath)
  
  
  ### list of FN
  FN <- setdiff(venne.genes[[method.ref]], venne.genes[[method]])
  cat("Number of FN: ",length(FN), "\n")
  
  ### order FN genes by significance
  table.FN <- table[table$gene_id %in% FN, c("gene_id", paste0("PValue.", method), paste0("adjPValue.", method))] 

  nas <- table(is.na(table.FN[, 2]))
  print(nas) 
  if(! "FALSE" %in% names(nas))
    next
  
  table.FN <- table.FN[order(table.FN[, 2], decreasing = TRUE), ]
  
  genes2plot <- table.FN[1:min(20, nas["FALSE"]), "gene_id"]
  
  plotPath <- file.path(out.dir.prop, paste0("ExonProportions_topFN_", gsub( "\\.", "_",method),".pdf"))
  
  plotProportions(dgeDM, genes2plot, plotPath)
  
  
}





#######################################################
# Run comparison in more automatic way... 
#######################################################



# model <- "diff_out_Full_tumorVSnormal"
# 
# count.methodList <- c("htseq", "kallisto")
# count.method <- count.methodList[1]
# 
# filter.methodList <- c("Filtering_DEXSeq", "Filtering_min6prop0_01min9cpm1maxNA")  
# filter.method <- filter.methodList[1]
# 
# 
# source("/home/gosia/R/R_Multinomial_project/Analysis_RealData_Kim_adenocarcinoma/Kim_comparison_DM_0_1_5.R")
















