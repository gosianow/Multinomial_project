#######################################################
# 
# Created 01 Sep 2014 

# comprare the DM version 3  commonDispersion and adjustement performance with other methods and DM dirmult

# Update 04 Sep 2014:

# + add DM version 4 


#######################################################
# BioC 2.14

setwd("/home/gosia/Multinomial_project/Simulations_drosophila_V2")




#######################################################
# merge Diff Expr results 
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


# problems when adding this results (same genes have multiple p-values)
# take min adj p-value 
rt <- read.table("Results_from_Katarina/cuffdiff_results.txt", header = T, stringsAsFactors = F)
head(rt)

rt <- rt[,c("ensembl_gene_id","p_value" ,"q_value")]
colnames(rt) <- c("Gene", "PValue_cuffdiff", "adjPValue_cuffdiff")

o <- order(rt$adjPValue_cuffdiff, decreasing = FALSE)

rt <- rt[o, ]
rt <- rt[!duplicated(rt$Gene, fromLast = FALSE), ]

results[["cuffdiff"]] <- rt





####################### DM v3 results / on fc


files <- fcFiles <- list.files(path = "DM_v3/fc/", pattern = "_results.xls" )
filtering <- filtering.fc <- gsub(pattern = "_results.xls", replacement = "", x = files)


for(i in seq_len(length(files))){
  
  rt <- read.table(paste0("DM_v3/fc/", files[i]), header = T, stringsAsFactors = F)
  head(rt)
  rt <- rt[,c("GeneID","PValue" ,"FDR")]
  colnames(rt) <- c("Gene", paste0("PValue_",filtering[i]), paste0("adjPValue_",filtering[i]))
  
  results[[filtering[i]]] <- rt  
  
}



####################### DM v3 results / on htseq


files <- htseqFiles <- list.files(path = "DM_v3/htseq/", pattern = "_results.xls" )
filtering <- filtering.htseq <- gsub(pattern = "_results.xls", replacement = "", x = files)

for(i in seq_len(length(files))){
  
  rt <- read.table(paste0("DM_v3/htseq/", files[i]), header = T, stringsAsFactors = F)
  print(head(rt))
  rt <- rt[,c("GeneID","PValue" ,"FDR")]
  colnames(rt) <- c("Gene", paste0("PValue_",filtering[i]), paste0("adjPValue_",filtering[i]))
  
  results[[filtering[i]]] <- rt  
  
}


####################### DM dirmult results / on fc / from durmult analysis


files <- fcDirmultFiles <- list.files(path = "DM_filtering/fc/", pattern = "_results.xls" )
filtering <- filtering.fcDirmult <- gsub(pattern = "_results.xls", replacement = "", x = files)

for(i in seq_len(length(files))){
  
  rt <- read.table(paste0("DM_filtering/fc/", files[i]), header = T, stringsAsFactors = F)
  head(rt)
  rt <- rt[,c("GeneID","PValue" ,"FDR")]
  colnames(rt) <- c("Gene", paste0("PValue_dirmult",filtering[i]), paste0("adjPValue_dirmult",filtering[i]))
  
  results[[paste0("dirmult", filtering[i])]] <- rt  
  
}

names(results)



####################### DM v4 results / on fc


files <- fcFiles <- list.files(path = "DM_v4/fc/", pattern = "_results.xls" )
filtering <- filtering.fc <- gsub(pattern = "_results.xls", replacement = "", x = files)


for(i in seq_len(length(files))){
  
  rt <- read.table(paste0("DM_v4/fc/", files[i]), header = T, stringsAsFactors = F)
  head(rt)
  rt <- rt[,c("GeneID","PValue" ,"FDR")]
  colnames(rt) <- c("Gene", paste0("PValue_",filtering[i]), paste0("adjPValue_",filtering[i]))
  
  results[[filtering[i]]] <- rt  
  
}



####################### DM v4 results / on htseq


files <- htseqFiles <- list.files(path = "DM_v4/htseq/", pattern = "_results.xls" )
filtering <- filtering.htseq <- gsub(pattern = "_results.xls", replacement = "", x = files)

for(i in seq_len(length(files))){
  
  rt <- read.table(paste0("DM_v4/htseq/", files[i]), header = T, stringsAsFactors = F)
  print(head(rt))
  rt <- rt[,c("GeneID","PValue" ,"FDR")]
  colnames(rt) <- c("Gene", paste0("PValue_",filtering[i]), paste0("adjPValue_",filtering[i]))
  
  results[[filtering[i]]] <- rt  
  
}

names(results)


#######################################################
# load information about simulation / add info about DS exons
#######################################################

simu_info.e <- read.table("Simu_info/true_exons_simulation.txt", header=T)
simu_info.g <- read.table("Simu_info/true_genes_simulation.txt", header=T)

head(simu_info.g)
head(simu_info.e)

simu_info.e[simu_info.e$Gene == "FBgn0001099", ]

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



###################### load simu_info

simu_info <- read.table("Simu_info/simu_info.xls", header=T)



#######################################################
# merge Diff Expr results with simu_info into one table
#######################################################


out.dir <- "PLOTS_DM_v4"
dir.create(out.dir, showWarnings=F, recursive=T)

simu_info <- read.table("Simu_info/simu_info.xls", header=T)
dim(simu_info)

table <- simu_info
for(i in 1:length(results)){
  print(dim(results[[i]]))
  table <- merge(table, results[[i]], by="Gene", all.x = TRUE)  
  table <- unique(table)
  print(dim(table))
}


write.table(table, paste0(out.dir,"/Table_all_results.xls"),  quote = F, sep = "\t", row.names = F, col.names = T)


## check table 

dim(table[!complete.cases(table), ])
colnames(table)

    


##############################################################################################################
# ->>>>>> load results
##############################################################################################################

out.dir <- "PLOTS_DM_v4/"

table <- read.table(paste0(out.dir,"/Table_all_results.xls"), header = T, stringsAsFactors = F)
tableNames <- colnames(table)

library(RColorBrewer)


allMethods <- gsub("adjPValue_", "",tableNames[grep(pattern = "adjPValue_", tableNames)])

ramp <- colorRampPalette(brewer.pal(9,"Set1"))(length(allMethods))

set.seed(42)
ord <- sample.int(length(allMethods))
ramp <- ramp[ord]


colors.org <- ramp
n.colors.org <- allMethods
names(colors.org) <- n.colors.org

colors.org[c("htseq_dexseq", "htseq_dexseq_simes", "fc_voomex",  "cuffdiff")] <- c("pink", "magenta", "blue",  "red")

colors <- colors.org
n.colors <- n.colors.org

name <- ""


# nr <- 9
# ramp <- colorRampPalette(brewer.pal(nr,"Set1"))(nr*3)
# plot(1:length(ramp), pch=19, cex=4, col=ramp)
# dev.off()


#######################################################
# histograms of p-values
#######################################################

colors <- colors.org
n.colors <- n.colors.org
colors <- colors.org[-1]
n.colors <- n.colors.org[-1]

pdf(paste0(out.dir, "/", name, "hist_pvalues.pdf"))


for(i in seq_len(length(colors))){
  # i=4
  cat(i, "\n")
  hist(table[, paste0("PValue_", n.colors[i]) ], col=colors[i], breaks=50, main=n.colors[i], cex.main=2, cex.lab=1.5, cex.axis = 1.5, xlab="P-values")
  
}

dev.off()



#######################################################
# generate ROCx plots 
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/ROCx.R")

colors <- colors.org
n.colors <- n.colors.org

n.colors

############ plots


methodsList <- list(
  all = c("htseq_dexseq_simes", "fc_voomex", "cuffdiff", "dirmultDM_fc_keep0s", "dirmultDM_fc_g0_s4_keep0s"),
  DM = c("htseq_dexseq_simes", "dirmultDM_fc_g0_s4_keep0s", "fc_g0_s4_keep0s_subsetInf_DM", "fc_g0_s4_keep0s_subsetInf_DMadj"),
  DMa = c("htseq_dexseq_simes", "fc_g0_s4_keep0s_subsetInf_DM", "fc_g0_s4_keep0s_subsetInf_DMadj", "fc_g0_s3_keep0s_subsetInf_DM", "fc_g0_s3_keep0s_subsetInf_DMadj"),
  DM4a = c("htseq_dexseq_simes", "fc_g0_s4_keep0s_subsetInf_DM4", "fc_g0_s4_keep0s_subsetInf_DM4adj", "fc_g0_s3_keep0s_subsetInf_DM4", "fc_g0_s3_keep0s_subsetInf_DM4adj"),
  htseq = c("htseq_dexseq_simes",  "fc_g0_s3_keep0s_subsetInf_DM", "htseq_g0_s3_keep0s_subsetInf_DMadj",  "fc_g0_s4_keep0s_subsetInf_DM", "htseq_g0_s4_keep0s_subsetInf_DMadj") )

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
  
  legend("bottomright", methods , col=colors[methods], lty=1, lwd=4, cex=1.5)
  dev.off()
  
  
}





#######################################################
# generate TPR vs achieved FDR plots
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/TPRvsFDR.R")


methodsList <- list(
  # all = c("htseq_dexseq", "htseq_dexseq_simes", "fc_voomex", "cuffdiff", "dirmultDM_fc_keep0s", "dirmultDM_fc_g0_s4_keep0s"),
  DM = c("htseq_dexseq_simes", "dirmultDM_fc_g0_s4_keep0s", "fc_g0_s4_keep0s_subsetInf_DM", "fc_g0_s4_keep0s_subsetInf_DMadj"),
  DMa = c("htseq_dexseq_simes", "fc_g0_s4_keep0s_subsetInf_DM", "fc_g0_s4_keep0s_subsetInf_DMadj", "fc_g0_s3_keep0s_subsetInf_DM", "fc_g0_s3_keep0s_subsetInf_DMadj"),
  DM4a = c("htseq_dexseq_simes", "fc_g0_s4_keep0s_subsetInf_DM4", "fc_g0_s4_keep0s_subsetInf_DM4adj", "fc_g0_s3_keep0s_subsetInf_DM4", "fc_g0_s3_keep0s_subsetInf_DM4adj"),
  htseq = c("htseq_dexseq_simes",  "fc_g0_s3_keep0s_subsetInf_DM", "htseq_g0_s3_keep0s_subsetInf_DM"),
  htseqDMa = c("htseq_dexseq_simes", "htseq_g0_s4_keep0s_subsetInf_DM",  "htseq_g0_s4_keep0s_subsetInf_DMadj", "htseq_g0_s3_keep0s_subsetInf_DM",  "htseq_g0_s3_keep0s_subsetInf_DMadj"),
  htseqDM4a = c("htseq_dexseq_simes", "htseq_g0_s4_keep0s_subsetInf_DM4",  "htseq_g0_s4_keep0s_subsetInf_DM4adj", "htseq_g0_s3_keep0s_subsetInf_DM4",  "htseq_g0_s3_keep0s_subsetInf_DM4adj") 
  )


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  pdf(paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf"))
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.2), ylim=c(0.5,0.8), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
  }
  
  legend("bottomright", methods, col=colors[methods], lty=1, lwd=4, cex=1.5, bty="o", bg="white", box.col="white")
  dev.off()
    
}



#######################################################
# generate venn diagrams 
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/plotVenn.R")

## TP as a separate circle
venne.genes <- list()
for(i in n.colors){ 
  venne.genes[[i]] <- na.omit(table[ table[,paste0("adjPValue_", i)] < FDR.cutoff, "Gene"]) 
}
venne.genes$True <- table[table$status == 1, "Gene"]


methodsList <- list(
  DM1 = c("htseq_dexseq_simes", "fc_voomex", "cuffdiff", "dirmultDM_fc_g0_s4_keep0s"),
  DM2 = c("htseq_dexseq_simes", "fc_voomex", "cuffdiff", "fc_g0_s4_keep0s_subsetInf_DM"),
  DM3 = c("htseq_dexseq_simes", "fc_voomex", "cuffdiff", "fc_g0_s4_keep0s_subsetInf_DMadj"),
  htseq1 = c("htseq_dexseq_simes", "fc_voomex", "cuffdiff", "htseq_g0_s4_keep0s_subsetInf_DM"),
  htseq2 = c("htseq_dexseq_simes", "fc_voomex", "cuffdiff",  "htseq_g0_s4_keep0s_subsetInf_DMadj") )


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  plotVenn(venne.genes, colors=c(colors[methods], True="grey"), venn.methods = c(methods, "True"), name1=name, name2=j, margin=0.1, cat.cex=1.4, cex=1.7)
  
}












































































