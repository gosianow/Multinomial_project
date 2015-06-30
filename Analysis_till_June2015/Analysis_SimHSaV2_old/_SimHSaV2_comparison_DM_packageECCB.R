#######################################################
# 
# Created 01 Sep 2014 

# comprare the DM version 3  commonDispersion and adjustement performance with other methods and DM dirmult

# Update 04 Sep 2014:

# + add DM version 4 


#######################################################
# BioC 2.14

setwd("/home/gosia/Multinomial_project/Simulations_hsapiens_V2")



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



# ####################### DM v4 results / on htseq
# 
# 
# files <- htseqFiles <- list.files(path = "DM_v4/htseq/", pattern = "_results.xls" )
# filtering <- filtering.htseq <- gsub(pattern = "_results.xls", replacement = "", x = files)
# 
# for(i in seq_len(length(files))){
#   
#   rt <- read.table(paste0("DM_v4/htseq/", files[i]), header = T, stringsAsFactors = F)
#   print(head(rt))
#   rt <- rt[,c("GeneID","PValue" ,"FDR")]
#   colnames(rt) <- c("Gene", paste0("PValue_",filtering[i]), paste0("adjPValue_",filtering[i]))
#   
#   results[[filtering[i]]] <- rt  
#   
# }
# 
names(results)




#######################################################
# merge Diff Expr results with simu_info into one table
#######################################################


out.dir <- "PLOTS_DM_v4/ECCB/"
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

out.dir <- "PLOTS_DM_v4/ECCB/"
dir.create(out.dir, recursive = T, showWarnings = FALSE)


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

colors.org[c("htseq_dexseq", "htseq_dexseq_simes", "fc_voomex",  "cuffdiff", "dirmultDM_fc_g0_s4_keep0s", "fc_g0_s4_keep0s_subsetInf_DM4", "fc_g0_s4_keep0s_subsetInf_DM4adj")] <- c("pink", "magenta", "blue",  "red", "slateblue", "chartreuse3", "orange")


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

pdf(paste0(out.dir, "/", name, "hist_pvalues-DEXSeqExon.pdf"))

rt <- read.table("Results_from_Katarina/dexseq_exon_results.txt", header = T, stringsAsFactors = F)
head(rt)

hist(rt[, paste0("pvalue") ], col="pink", breaks=50, main="DEXSeq exon level", cex.main=2, cex.lab=1.5, cex.axis = 1.5, xlab="P-values")


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
  all = c("htseq_dexseq_simes", "cuffdiff", "dirmultDM_fc_g0_s4_keep0s", "fc_g0_s4_keep0s_subsetInf_DM4", "fc_g0_s4_keep0s_subsetInf_DM4adj") )

methodsListLeg <- list(
  all = c("DEXSeq (simes)", "Cuffdiff", "DM tagwise", "DM CPL", "DM ACPL") )

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
  
  legend("bottomright", methodsLeg , col=colors[methods], lty=1, lwd=4, cex=1.5)
  dev.off()
  
  
}





#######################################################
# generate TPR vs achieved FDR plots
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/TPRvsFDR.R")


methodsList <- list(
  all = c("htseq_dexseq" ,"htseq_dexseq_simes", "cuffdiff", "dirmultDM_fc_g0_s4_keep0s", "fc_g0_s4_keep0s_subsetInf_DM4", "fc_g0_s4_keep0s_subsetInf_DM4adj") )

methodsListLeg <- list(
  all = c("DEXSeq", "DEXSeq (simes)", "Cuffdiff", "DM tagwise", "DM CPL", "DM ACPL") )


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  
  pdf(paste0(out.dir, "/", name, "TPRachievedFDR_",j,".pdf"))
  
  apvs <- table[,paste0("adjPValue_", methods[1])]
  status <- table$status
  TPRvsFDR(status, apvs, col=colors[methods[1]], xlim=c(0,0.6), ylim=c(0.8,0.95), lwd=4, cex=2.5, cex.lab=1.45, cex.axis=1.5)
  
  for( i in methods[-1]){
    apvs <- table[,paste0("adjPValue_", i)]
    status <- table$status
    TPRvsFDR(status, apvs, col=colors[i], add=TRUE, lwd=4, cex=2.5)
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
  venne.genes[[i]] <- na.omit(table[ table[,paste0("adjPValue_", i)] < FDR.cutoff, "Gene"]) 
}
venne.genes$True <- table[table$status == 1, "Gene"]


methodsList <- list(
  all = c("htseq_dexseq", "cuffdiff", "fc_g0_s4_keep0s_subsetInf_DM4adj") )

methodsListLeg <- list(
  all = c("DEXSeq", "Cuffdiff", "DM ACPL") )


for(j in names(methodsList)){
  
  methods <- methodsList[[j]]
  methodsLeg <- methodsListLeg[[j]]
  
  plotVenn(venne.genes, colors=c(colors[methods], True="grey"), venn.methods = c(methods, "True"), methodsLeg = c(methodsLeg, "TRUE"),  name1=name, name2=j, margin=0.1, cat.cex=2, cex=1.7)
  
}












































































