#######################################################

# Created 22 June 2015
# BioC 3.1

# Comprare DM_0.1.5 (different dispersion estimators) with DEXSeq
# Use htseq and kallisto counts
# Plots per count method and per filtering


#######################################################


# setwd("/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull")

# setwd("/home/gosia/Multinomial_project/Simulations_hsapiens_Sim5_noDE_noNull")

# setwd("/home/gosia/Multinomial_project/Simulations_hsapiens_Sim5_withDE_noNull")



library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)


### Source all R files in DM package / DM_0.1.5
Rfiles <- list.files("/home/gosia/R/R_Multinomial_project/DM_package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])

DM_out <- "DM_0_1_5/"
DM_plots_out <- "DM_0_1_5_PLOTS/"



##############################################################################################################
# load metadata
##############################################################################################################

# create metadata file
metadata <- data.frame(sampleName = paste0("sample_",1:6), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata



#######################################################
# load information about simulation 
#######################################################


simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, as.is = TRUE)

simulation_details <- simulation_details[, c("gene_id", "transcript_id", "expected_count", "diff_IsoPct", "gene_ds_status", "transcript_ds_status", "gene_de_status")]



#######################################################
# merge Diff Expr results 
#######################################################


# count.methodList <- c("htseq", "kallisto")
# count.method <- count.methodList[1]
# 
# filter.methodList <- c("Filtering_DEXSeq", "Filtering_min3prop0_01min6cpm1maxNA", "Filtering_min3prop0_05min6cpm1maxNA") 
# filter.method <- filter.methodList[1]



results <- list()

####################### results produced by Charlotte


rt <- read.table(paste0("4_results/dexseq_", count.method, ".txt"), header = T, as.is = TRUE)
head(rt)

colnames(rt) <- c("gene_id",  paste0("adjPValue.",count.method,"_dexseq"))

head(rt)
results[[paste0( count.method, "_dexseq" )]] <- rt




####################### DM_0.1.5 results on htseq counts and bitseq 



res.path <- paste0(DM_out, count.method, "/", filter.method, "/")

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
# merge Diff Expr results with simu_info into one table
#######################################################

plots.path <- paste0(DM_plots_out, count.method, "/", filter.method, "/")
dir.create(plots.path, showWarnings=F, recursive=T)


simulation_details <- unique(simulation_details[, c("gene_id", "gene_ds_status")])

table <- merge(simulation_details, results, by = "gene_id", all.x = TRUE)  

write.table(table, paste0(plots.path,"/table_all_results.xls"),  quote = F, sep = "\t", row.names = F, col.names = T)



##############################################################################################################
# ->>>>>> load results
##############################################################################################################

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


# count.methodList <- c("htseq", "kallisto")
# count.method <- count.methodList[1]
# 
# filter.methodList <- c("Filtering_DEXSeq", "Filtering_min3prop0_01min6cpm1maxNA", "Filtering_min3prop0_05min6cpm1maxNA") 
# filter.method <- filter.methodList[1]



plots.path <- paste0(DM_plots_out, count.method, "/", filter.method, "/")

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
# generate ROCx plots 
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/ROCx_0_1_4.R")

colors <- allMethods[grepl("DM",allMethods$names),]


methodsListLeg <- methodsList <- list(
  DM = colors[grepl("DM", colors$names), ])


for(j in names(methodsList)){
  # j = names(methodsList)[1]
  
  methods <- methodsList[[j]]
  
  scores <- lapply(methods$names, function(i){
    print(i)
    df <- data.frame(pvs = table[, paste0("PValue.", i)], apvs = table[ ,paste0("adjPValue.", i)], status = table$gene_ds_status)
    return(df)
    })
names(scores) <- methods$names
  colours <- methods$colors
  names(colours) <- methods$names
  
  plotPath =  paste0(plots.path, "/", name, "rocX_notNorm_", j ,".pdf")

  ROCx_notNorm(scores = scores, colours = colours, plotPath = plotPath, xlim=c(0,1), ylim=c(0,1))
  
}








#######################################################
# generate TPR vs achieved FDR plots
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/TPRvsFDR_0_1_4.R")

colors <- allMethods

methodsListLeg <- methodsList <- list(
  DM = colors)

lim <- list(DM = data.frame(x = c(0, 0.8), y = c(0, 1)))


for(j in names(methodsList)){
  # j = names(methodsList)[1]
  
  methods <- methodsList[[j]]
  
  scores <- list()
  
  for(i in methods$names)
    scores[[i]] <- data.frame(apvs = table[ ,paste0("adjPValue.", i)], status = table$gene_ds_status)
  
  names(scores) <- methods$names
  colours <- methods$colors
  names(colours) <- methods$names
  
  plotPath =  paste0(plots.path, "/", name, "TPRachievedFDR_",j,".pdf")
  
  TPRvsFDR(scores = scores, colours = colours, plotPath = plotPath, xlim = lim[[j]][, "x"], ylim = lim[[j]][, "y"])
  
  
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
    
    main <- paste0(gene, "\n Mean Expression = ", round(dgeDM$meanExpr[gene]), " / Dispersion = ", round(dgeDM$fitNull[[gene]]$gamma0, 2), "\n LR = ", round(dgeDM$table[gene, "LR"], 4) , " / P-value = ", sprintf("%.02e", dgeDM$table[gene, "pValue"]), " / FDR = ", sprintf("%.02e", dgeDM$table[gene, "FDR"]))
    
    
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
  venne.genes[[i]] <- na.omit(table[ table[,paste0("adjPValue.", i)] < FDR.cutoff, "gene_id"]) 
}
venne.genes$True <- na.omit(table[table$gene_ds_status == 1, "gene_id"])


res.path <- paste0(DM_out, count.method, "/", filter.method, "/")

DM.methods <- allMethods[grepl(count.method, allMethods$names) & grepl("DM", allMethods$names) & !grepl("_commonDispersion", allMethods$names), "names"]
DM.methods


for(m in 1:length(DM.methods)){
  # m = 1
  
  method <- DM.methods[m]
  method
  method.ref <- paste0(count.method, "_dexseq")
  method.ref
  
  ### list of FP
  FP <- setdiff(setdiff(venne.genes[[method]], venne.genes[["True"]]), venne.genes[[method.ref]])
  # FP2plot <- setdiff(venne.genes[[method]], venne.genes[["True"]])
  cat("Number of FP: ",length(FP), "\n")
  
  ### order FP genes by significance
  table.FP <- table[table$gene_id %in% FP, c("gene_id", paste0("adjPValue.", method))]
  table.FP <- table.FP[order(table.FP[, 2], decreasing = FALSE), ]
  
  genes2plot <- table.FP[1:20, "gene_id"]
  
  ### load DM pipeline results
  dgeDM.rdata <- file.path(res.path, paste0(gsub( "\\.", "-",method), "_dgeDM.RData"))
  print(dgeDM.rdata)
  load(dgeDM.rdata)
  # dgeDM
  
  
  plotPath <- file.path(out.dir.prop, paste0("ExonProportions_topFP_", gsub( "\\.", "_",method),".pdf"))
  
  plotProportions(dgeDM, genes2plot, plotPath)
  
  
  ### list of FN
  FN <- intersect(setdiff(venne.genes[["True"]], venne.genes[[method]]), venne.genes[[method.ref]])
  cat("Number of FN: ",length(FN), "\n")
  
  ### order FN genes by significance
  table.FN <- table[table$gene_id %in% FN, c("gene_id", paste0("adjPValue.", method))]   
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
venn.list$True <- na.omit(table[table$gene_ds_status == 1, "gene_id"])

names(venn.list)


colors <- allMethods[grepl("DM", allMethods$names), ]


methodsListLeg <- methodsList <- list()
for(i in colors$names){
  methodsListLeg[[i]] <- methodsList[[i]] <- c(paste0(count.method, "_dexseq"), i)
}


colors <- allMethods

for(j in names(methodsList)){
  # j = names(methodsList)[1]
  methods <- methodsList[[j]]
  methodsLeg <- methodsListLeg[[j]]

  plotPath =  paste0(plots.path, "/", name, "Venn_Diagram_", gsub( "\\.", "_",j) ,".pdf")
  
  venn.colors <- colors[methods, "colors"]
  names(venn.colors) <- methods
  venn.colors <- c(True = "grey", venn.colors)
  
  plotVenn2(venn.list, venn.colors = venn.colors, venn.subset = c("True", methods), margin = 0, cat.cex=1.2, cex=1.7, plotPath = plotPath)
  
}





#######################################################
# Run comparison in more automatic way... 
#######################################################


# # setwd("/home/gosia/Multinomial_project/Simulations_drosophila_Sim5_noDE_noNull")
# 
# # setwd("/home/gosia/Multinomial_project/Simulations_hsapiens_Sim5_noDE_noNull")
# 
# setwd("/home/gosia/Multinomial_project/Simulations_hsapiens_Sim5_withDE_noNull")
# 
# 
# 
# 
# count.methodList <- c("htseq", "htseq", "kallisto", "kallisto", "kallisto")
# 
# filter.methodList <- c("Filtering_DEXSeq", "Filtering_min3prop0_01min6cpm1maxNA", "Filtering_DEXSeq", "Filtering_min3prop0_01min6cpm1maxNA", "Filtering_min3prop0_05min6cpm1maxNA") 
# 
# 
# 
# for(i in 3:length(count.methodList)){
#   
#   count.method <- count.methodList[i]
#   filter.method <- filter.methodList[i]
#   
#   source("/home/gosia/R/R_Multinomial_project/Analysis_Sim5_noDE_noNull/Compare_DM_0.1.5.R")
#   
# }




























