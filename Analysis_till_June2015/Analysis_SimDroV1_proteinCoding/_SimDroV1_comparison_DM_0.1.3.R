#######################################################

# Created 26 May 2015
# BioC 3.0

# comprare the DM_0.1.3: tagwise Dipsersion and filtering; using htseq counts and bitseq counts
# For DM: Filtering_DEXSeq


#######################################################


library(DM)

library(limma)
source("/home/gosia/R/R_Multinomial_project/DM_package_devel/0_my_printHead.R")

library(edgeR)





setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V1_proteinCoding")

# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata




#######################################################
# load information about simulation 
#######################################################

simu_info.e <- read.table("Simu_info/true_exons_simulation.txt", header=T)
simu_info.g <- read.table("Simu_info/true_genes_simulation.txt", header=T)

head(simu_info.g)
head(simu_info.e)

simu_info1 <- simu_info.g[!duplicated(simu_info.g$gene_id),c("gene_id", "status", "num", "rpk.mean")]

### count number of defferentially spliced exons

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




####################### DM_0.1.3 results on htseq counts and bitseq 


filter.methodList <- c("Filtering_DEXSeq", "Filtering0", "Filtering001") 
filter.method <- filter.methodList[1]


count.methodList <- c("htseq", "bitseq")

for(count.method in count.methodList){
  # count.method = count.methodList[1]; filter.method = filter.methodList[1]
  
  res.path <- paste0("DM_0_1_3/", count.method, "/", filter.method, "/")
  
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
# merge Diff Expr results with simu_info into one table
#######################################################


out.dir <- paste0("PLOTS_DM_0_1_3/", filter.method,"/")
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


filter.methodList <- c("Filtering_DEXSeq", "Filtering0", "Filtering001") 
filter.method <- filter.methodList[1]

out.dir <- paste0("PLOTS_DM_0_1_3/", filter.method,"/")



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

set.seed(10)

ramp <- gg_color_hue(nrow(allMethods))[sample.int(nrow(allMethods))]

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
  
#   df <- data.frame(PValues = table[, paste0("PValue_", colors$names[i]) ])
#   
#   ggp <- ggplot(df, aes(x = PValues)) +
#     ggtitle(colors$names[i]) +
#     geom_histogram(binwidth = 0.01, fill=colors$colors[i]) +
#     theme(text = element_text(size=10))
#    
#   print(ggp)
  
  PValues = table[, paste0("PValue_", colors$names[i]) ] 
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
  
  scores <- list()
  
  for(i in methods$names)
  scores[[i]] <- data.frame(pvs = table[, paste0("PValue_", i)], apvs = table[ ,paste0("adjPValue_", i)], status = table$status)
  
  names(scores) <- methods$names
  colours <- methods$colors
  names(colours) <- methods$names
  
  plotPath =  paste0(out.dir, "/",name,"rocX_notNorm_",j,".pdf")
  
  # ggROCx_notNorm(scores = scores, colours = colours, plotPath = plotPath)
  ROCx_notNorm(scores = scores, colours = colours, plotPath = plotPath, xlim=c(0,0.8), ylim=c(0.8,1))
  
}








#######################################################
# generate TPR vs achieved FDR plots
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/TPRvsFDR_0_1_4.R")


colors <- allMethods


methodsListLeg <- methodsList <- list(
  all = colors, 
  htseq = colors[grepl("htseq", colors$names), ],
  bitseq = colors[grepl("bitseq", colors$names), ])


lim <- list(all = data.frame(x = c(0, 0.3), y = c(0.4, 1)), htseq = data.frame(x = c(0, 0.2), y = c(0.8, 1)), bitseq = data.frame(x = c(0, 0.3), y = c(0.9, 1)))


for(j in names(methodsList)){
  # j = names(methodsList)[1]
  
  methods <- methodsList[[j]]
  
  scores <- list()
  
  for(i in methods$names)
    scores[[i]] <- data.frame(apvs = table[ ,paste0("adjPValue_", i)], status = table$status)
  
  names(scores) <- methods$names
  colours <- methods$colors
  names(colours) <- methods$names
  
  plotPath =  paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf")
  
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
  
  
  pdf(plotPath, width = 12, height = 7)
  
  
  for(g in 1:length(genes2plot)){
    # g = 1
    
    gene <- genes2plot[g]
    # print(gene)
    Condition <- dgeDM$samples$group 
    expr <- dgeDM$counts[[gene]]
    
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
# plot expression ratios of FP called by DM
#######################################################





FDR.cutoff <- 0.05


### generate list of significant genes for each method
venne.genes <- list()
for(i in allMethods$names){ 
  venne.genes[[i]] <- na.omit(table[ table[,paste0("adjPValue_", i)] < FDR.cutoff, "gene_id"]) 
}
venne.genes$True <- na.omit(table[table$status == 1, "gene_id"])


for(count.method in count.methodList){
  # count.method = "htseq"
  res.path <- paste0("DM_0_1_4/", count.method, "/", filter.method, "/")
  
  DM.methods <- allMethods[grepl(count.method, allMethods$names) & grepl("DM", allMethods$names), "names"]
  
  
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
    table.FP <- table[table$gene_id %in% FP, c("gene_id", paste0("adjPValue_", method))]
    table.FP <- table.FP[order(table.FP[, 2], decreasing = FALSE), ]
    
    genes2plot <- table.FP[1:20, "gene_id"]
    
    ### load DM pipeline results
    dgeDM.rdata <- file.path(res.path, paste0(gsub( "\\.", "-",method), "_dgeDM.RData"))
    print(dgeDM.rdata)
    load(dgeDM.rdata)
    # dgeDM
    
    rownames(dgeDM$table) <- dgeDM$table$GeneID

    plotPath <- file.path(out.dir, paste0("ExonProportions_topFP_", gsub( "\\.", "_",method),".pdf"))
    
    plotProportions(dgeDM, genes2plot, plotPath)
    
    
    ### list of FN
    FN <- intersect(setdiff(venne.genes[["True"]], venne.genes[[method]]), venne.genes[[method.ref]])
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




































