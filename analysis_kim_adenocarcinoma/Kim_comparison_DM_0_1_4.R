######################################################
# BioC 3.0

# Created 11 June 2015 
# Comparison of Brooks_pasilla data with DM_0.1.4 using htseq counts with DEXSeq

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

Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM_0.1.4/R/", full.names=TRUE)
for(i in Rfiles) source(i)


##########################################################################
# load metadata
##########################################################################

metadata <- read.table("metadata/Malgorzata Nowicka2014-11-04GSE37764.csv", stringsAsFactors=F, sep=",", header=T) 

metadata <- metadata[metadata$X == "RNA-seq",]

metadata


#######################################################
# merge Diff Expr results 
#######################################################

results <- list()

####################### results produced by Katarina

### DEXSeq
rt <- read.table("DEXSeq_1.10.8/diff_out_Full_tumorVSnormal/DEXSeq_gene_results.txt", header = T, as.is = TRUE)
head(rt)

rt <- rt[,c("geneID", "pvalue")]
colnames(rt) <- c("Gene",  "adjPValue_htseq_dexseq")

results[["htseq_dexseq"]] <- rt

sum(results[["htseq_dexseq"]]$adjPValue_htseq_dexseq < 0.05)



####################### DM_0.1.4 results on htseq counts and bitseq 


filter.methodList <- c("Filtering_DEXSeq", "Filtering0", "Filtering001") 
filter.method <- filter.methodList[1]


count.methodList <- c("htseq", "bitseq")[1]

for(count.method in count.methodList){
  # count.method = count.methodList[1]; filter.method = filter.methodList[1]
  
  res.path <- paste0("DM_0_1_4/", count.method, "/", filter.method, "/")
  
  files <- list.files(path = res.path, pattern = "_results.xls" )
  res  <- gsub(pattern = "_results.xls", replacement = "", x = files)
  
  for(i in seq_len(length(files))){
    # i = 1
    rt <- read.table(paste0(res.path, files[i]), header = T, stringsAsFactors = F)
    head(rt)
    rt <- rt[,c("GeneID","PValue" ,"FDR")]
    colnames(rt) <- c("Gene", paste0("PValue_",res[i]), paste0("adjPValue_",res[i]))
    
    results[[res[i]]] <- rt  
    
  }
  
}


names(results)




#######################################################
# merge Diff Expr results into one table
#######################################################


out.dir <- paste0("PLOTS_DM_0_1_4/", filter.method,"/")
dir.create(out.dir, showWarnings=F, recursive=T)

table <- results[[1]]

for(i in 2:length(results)){
  print(dim(results[[i]]))
  table <- merge(table, results[[i]], by="Gene",  all = TRUE)  
  table <- unique(table)
  print(dim(table))
}

names(table)[1] <- "gene_id"

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


filter.methodList <- c("Filtering_DEXSeq", "Filtering0", "Filtering001") 
filter.method <- filter.methodList[1]

out.dir <- paste0("PLOTS_DM_0_1_4/", filter.method,"/")



### load table with all results

table <- read.table(paste0(out.dir,"/Table_all_results.xls"), header = T, stringsAsFactors = F)
tableNames <- colnames(table)


### assign colors to different methods

allMethods <- data.frame(names = gsub("adjPValue_", "",tableNames[grep(pattern = "adjPValue_", tableNames)]), stringsAsFactors = FALSE)
rownames(allMethods) <- allMethods$names

### generate ggplot hue colors 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

ramp <- gg_color_hue(nrow(allMethods))

allMethods$colors <- ramp

name <- ""





#######################################################
# histograms of p-values
#######################################################


colors <- allMethods[grepl("DM",allMethods$names),]


pdf(paste0(out.dir, "/", name, "hist_pvalues.pdf"), 5, 5)

for(i in seq_len(nrow(colors))){
  # i=1
  cat(i, "\n")
  PValues = table[, paste0("PValue_", colors$names[i]) ] 
  hist(PValues, breaks = 100, main = paste0(colors$names[i]), col = colors$colors[i], cex.main = 1)
  
  
}

dev.off()



#######################################################
# generate venn diagrams 
#######################################################


source("/home/gosia/R/R_Multinomial_project/Plot_Functions/plotVenn_0_1_4.R")

name <- ""

colors <- allMethods

FDR.cutoff <- 0.05



## TP as a separate circle
venn.list <- list()

for(i in colors$names){ 
  venn.list[[i]] <- na.omit(table[ table[,paste0("adjPValue_", i)] < FDR.cutoff, "gene_id"]) 
}

names(venn.list)




count.method <- "htseq"
# count.method <- "bitseq"

colors <- allMethods[grepl(count.method, allMethods$names) & grepl("DM", allMethods$names), ]

methodsListLeg <- methodsList <- list()

for(i in colors$names){
  methodsListLeg[[i]] <- c("DEXSeq", "DM")
    methodsList[[i]] <- c(paste0(count.method, "_dexseq"), i)
}


colors <- allMethods

for(j in names(methodsList)){
  # j = names(methodsList)[1]
  methods <- methodsList[[j]]
  methodsLeg <- methodsListLeg[[j]]
  
  plotPath =  paste0(out.dir, "/", name, "Venn_Diagram_",j,".pdf")
  
  venn.colors <- colors[methods, "colors"]
  names(venn.colors) <- methods
  
  venn.list_tmp <- venn.list[methods]
  names(venn.list_tmp) <- methodsLeg
  names(venn.colors) <- methodsLeg
  
  plotVenn2(venn.list = venn.list_tmp, venn.colors = venn.colors, venn.subset = methodsLeg, margin = 0, cat.cex=2, cex=1.7, plotPath = plotPath)
  
  
}




##################################################################################
################################# plot proportions
##################################################################################
library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


plotProportions <- function(dgeDM, genes2plot, plotPath){
  
  
  pdf(plotPath, width = 12, height = 7)
  
  
  for(g in 1:length(genes2plot)){
    # g = 1
    
    gene <- genes2plot[g]
    # print(gene)
    oo <- order(dgeDM$samples$group)
    Condition <- dgeDM$samples$group[oo]

    expr <- dgeDM$counts[[gene]][,oo]
    
    labels <- strsplit2(rownames(expr), ":")[, 2]
    
    prop.smp <- data.frame( ete_id =  labels, prop.table(expr, 2))  
    prop.est <- data.frame(ete_id = labels, dgeDM$fit[[gene]]$piH)
    
    prop.est.null <- data.frame(ete_id = labels, dgeDM$fit.null[[gene]]$piH)
    
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
    
    main <- paste0(gene, "\n Mean Expression = ", round(dgeDM$meanExpr[gene]), " / Dispersion = ", round(dgeDM$fit.null[[gene]]$gamma0, 2), "\n LR = ", round(dgeDM$table[gene, "LR"], 4) , " / P-value = ", sprintf("%.02e", dgeDM$table[gene, "PValue"]), " / FDR = ", sprintf("%.02e", dgeDM$table[gene, "FDR"]))
    
    
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



FDR.cutoff <- 0.05


### generate list of significant genes for each method
venne.genes <- list()
for(i in allMethods$names){ 
  venne.genes[[i]] <- na.omit(table[ table[,paste0("adjPValue_", i)] < FDR.cutoff, "gene_id"]) 
}


for(count.method in count.methodList){
  # count.method = "htseq"
  res.path <- paste0("DM_0_1_4/", count.method, "/", filter.method, "/")
  
  DM.methods <- allMethods[grepl(count.method, allMethods$names) & grepl("DM", allMethods$names), "names"]
  
  
  for(m in 1:length(DM.methods)){
    # m = 3
    
    method <- DM.methods[m]
    method
    method.ref <- paste0(count.method, "_dexseq")
    method.ref
    
    ### list of FP
    FP <- setdiff(venne.genes[[method]], venne.genes[[method.ref]])
    # FP2plot <- setdiff(venne.genes[[method]], venne.genes[["True"]])
    cat("Number of FP: ",length(FP), "\n")
    
    ### order FP genes by significance
    table.FP <- table[table$gene_id %in% FP, c("gene_id", paste0("adjPValue_", method))]
    table.FP <- table.FP[order(table.FP[, 2], decreasing = FALSE), ]
    
    genes2plot <- table.FP[1:min(20, nrow(table.FP)), "gene_id"]
    
    ### load DM pipeline results
    dgeDM.rdata <- file.path(res.path, paste0(gsub( "\\.", "-",method), "_dgeDM.RData"))
    print(dgeDM.rdata)
    load(dgeDM.rdata)
    # dgeDM
    
    rownames(dgeDM$table) <- dgeDM$table$GeneID
    
    plotPath <- file.path(out.dir, paste0("ExonProportions_topFP_", gsub( "\\.", "_",method),".pdf"))
    
    plotProportions(dgeDM, genes2plot, plotPath)
    
    
    ### list of FN
    FN <- setdiff(venne.genes[[method.ref]], venne.genes[[method]])
    cat("Number of FN: ",length(FN), "\n")
    
    ### order FN genes by significance
    table.FN <- table[table$gene_id %in% FN, c("gene_id", paste0("adjPValue_", method))]   
    nas <- table(is.na(table.FN[, 2]))
    print(nas) 
    if(! "FALSE" %in% names(nas))
      next
    
    table.FN <- table.FN[order(table.FN[, 2], decreasing = TRUE), ]
    
    genes2plot <- table.FN[1:min(20, nas["FALSE"]), "gene_id"]
    
    plotPath <- file.path(out.dir, paste0("ExonProportions_topFN_", gsub( "\\.", "_",method),".pdf"))
    
    plotProportions(dgeDM, genes2plot, plotPath)
    
    
  }
  
}








