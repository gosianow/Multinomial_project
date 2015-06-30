# BioC 3.0

# Created 8 Jan 2014
# Modyfied 27 Jan 2014


##############################################################################################################

# compare DM SQTL with sQTLseekeR results

##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")


resSeeker <- read.table("sQTLseekeR-GEUVADIS-FDR10-annotation/sQTLs-FDR10-CEU.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

res.seeker <- resSeeker[resSeeker$FDR < 0.05, ] ## 2869


resDM <- read.table("DMv5_sQTL_analysis/dgeSQTL_results.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

res.dm <- resDM[resDM$FDR < 0.05 & !is.na(resDM$SNP_id), ] ## 205198


#######################################################
# generate venn diagrams 
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/plotVenn.R")

## TP as a separate circle
venne.list <- list()
venne.list[["Seeker"]] <- resSeekerSign[,"snpId"]
# venne.list[["Seeker"]] <- intersect(resSeekerSign[,"snpId"], resDM[, "SNP_id"])
venne.list[["DM"]] <- na.omit(resDMSign[,"SNP_id"])


  plotVenn(venne.list, colors=c(Seeker="blue", DM="orange"), venn.methods = names(venne.list), margin=0.1, cat.cex=0.8, cex=1.7, out.dir="", name2="")
  


#######################################################
# generate hist of p-values
#######################################################

out.dir <- "DMv5_sQTL_analysis/Comparison/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


pdf(paste0(out.dir, "hist_pvalues.pdf"))

hist(resDM[, "PValue"], col = "#FF7F00", breaks = 100, cex.lab=1.5, cex.axis = 1.5, xlab="P-values", main = "DM", cex.main = 3)

dev.off()



##############################################################################################################

# Plot expression for most significant snp-gene 

##############################################################################################################


#######################################################
# load data
#######################################################

library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v5.R"))

library(ggplot2)
library(reshape2)
library(gridExtra)

#### load dgeSQTL
load("DMv5_sQTL_analysis/dgeSQTL.RData")

out.dir <- "DMv5_sQTL_analysis/Comparison/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)



#######################################################
#### Find which snp-gene are FP
#######################################################

res.seeker$SNPgene <- paste0(res.seeker[,"snpId"], "-", res.seeker[,"geneId"])
res.dm$SNPgene <- paste0(res.dm[,"SNP_id"], "-", res.dm[, "gene_id"])


fp.dm <- setdiff(res.dm$SNPgene, res.seeker$SNPgene)

res.dm.o <- res.dm[order(res.dm$PValue, decreasing = FALSE), ]

res.dm.fp <- res.dm.o[res.dm.o$SNPgene %in% fp.dm, ]

### number of genes with x number of snps

res.dm.t <- table(res.dm.fp$gene_id)
res.dm.tt <- table(res.dm.t)
res.dm.tt <- data.frame(Number_of_SNPs = names(res.dm.tt), Number_of_genes = as.numeric(res.dm.tt))

write.table(res.dm.tt, paste0(out.dir, "Number_of_SNPs_per_gene_DM.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


res.dm.fp <- split(res.dm.fp, f = factor(res.dm.fp$gene_id, levels = unique(res.dm.fp$gene_id) ))

### keep genes with large number of snps > 100
gene <- names(which(res.dm.t > 100))

res.dm.fp <- res.dm.fp[gene] 


### keep top snp from each gene
res.dm.fp.top <-  lapply(seq(length(res.dm.fp)), function(g){ res.dm.fp[[g]][1, , drop = FALSE] })

res.dm.fp.top <- do.call(rbind, res.dm.fp.top)
res.dm.fp.top <- res.dm.fp.top[order(res.dm.fp.top$PValue, decreasing = FALSE), ]

res.dm.fp.top$rank <- 1:nrow(res.dm.fp.top)
res.dm.fp.top$nrSNPs <- as.numeric(res.dm.t[res.dm.fp.top$gene_id])




##############################################################
# Plot raw expression for FP called by DM
##############################################################

out.dir <- "DMv5_sQTL_analysis/Comparison/DM_FP/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)


plot.snps <- res.dm.fp.top[, c("gene_id", "SNP_id", colnames(res.dm.fp.top)[-c(1, 2)])]
dim(plot.snps)

plot.names <- paste0("_rank", res.dm.fp.top$rank, "_nrSNPs", res.dm.fp.top$nrSNPs, "_df", res.dm.fp.top$df)
plot.main <- paste0(plot.names)



for(i in 1:20){
  # i = 1
  cat(paste0("SNPgene ", i, "\n"))
  
  gene <- plot.snps[i, 1]
  snp <- plot.snps[i, 2]
  
  expr <- data.frame(dgeSQTL$counts[[gene]])
  
  snpGene <- which(dgeSQTL$SNPs$SNP_id == snp & dgeSQTL$SNPs$gene_id == gene )[1]
  
  geno <- dgeSQTL$genotypes[ snpGene , ]
  
  samps.keep <- !is.na(geno) & !is.na(expr[1,])
  
  expr <- expr[,samps.keep]
  geno <- geno[samps.keep]
 
  
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  gg_color_variants <- function(n) {
    c("black", "grey50", "grey90")[1:n]
  }
  
  min.mean.sd.max <- function(x) {
    r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  
  
  ######## Plot Expression 
  tom <- cbind(expr, Transcript=rownames(expr))  
  
  m <- melt(tom, id.vars = "Transcript" )
  
  m$genotype <- NA
  
  geno.val <- sort(unique(geno))  
  var.counts <- rep(0, length(geno.val))
  names(var.counts) <- paste0("variant ", geno.val)
  
  for(j in 1:length(geno.val)){    
    m$genotype[m$variable %in% names(geno[geno == geno.val[j]] )] <- paste0("variant ", geno.val[j])
    var.counts[j] <- length(geno[geno == geno.val[j]])   
  }
  

  ggp.expr <- ggplot(data = m, aes(x = genotype, y = value, fill = Transcript)) + 
    geom_boxplot( width = 1) +
    # stat_summary(fun.y= min.mean.sd.max, geom='point', colour = "red", position = position_dodge(width = 1)) + 
    ggtitle(paste0(snp, "-", gene, "\n", plot.main[i])) +
    scale_x_discrete(labels = paste0(names(var.counts), " (",var.counts, ")" ), name="") + 
    ylab("Expression (counts)") + 
    theme(panel.grid.major=element_blank(), axis.text=element_text(size=13),axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(colour = gg_color_variants(length(var.counts)))) + 
    geom_vline(xintercept=c(1.5,2.5),color="white")
  
    ggsave(paste0(out.dir, "DM_expr",plot.names[i],"_", snp, "-", gene, ".pdf"), width = 15, height = 7, units = "in", plot = ggp.expr)
  
  
  ##### Plot violin

# ggp.expr.v <- ggplot(data = m, aes(x = genotype, y = value, fill = Transcript)) + 
#   geom_violin(width = 1) +
#   stat_summary(fun.y= mean, geom='point', position = position_dodge(width = 1)) +
#   ggtitle(paste0(snp, "-", gene, "\n", plot.main[i])) +
#   scale_x_discrete(labels = paste0(names(var.counts), " (",var.counts, ")" ), name="") + 
#   ylab("Expression (counts)") + 
#   theme(panel.grid.major=element_blank(), axis.text=element_text(size=13),axis.title=element_text(size=14,face="bold")) + 
#   geom_vline(xintercept=c(1.5,2.5),color="white")
# 
# ggsave(paste0(out.dir, "DM_expr",plot.names[i],"_", snp, "-", gene, ".pdf"), width = 15, height = 7, units = "in", plot = ggp.expr.v)




######## Plot proportions per sample

tot <- colSums(expr)
prop.smp <- data.frame(t(apply(expr, 1, function(t){ t / tot })))

tom <- cbind(prop.smp, Transcript=rownames(expr))  

m <- melt(tom, id.vars = "Transcript" )

m$genotype <- NA

geno.val <- sort(unique(geno))  
var.counts <- rep(0, length(geno.val))
names(var.counts) <- paste0("variant ", geno.val)

for(j in 1:length(geno.val)){    
  m$genotype[m$variable %in% names(geno[geno == geno.val[j]] )] <- paste0("variant ", geno.val[j])
  var.counts[j] <- length(geno[geno == geno.val[j]])   
}


ggp.prop.smp <- ggplot(data = m, aes(x = genotype, y = value, fill = Transcript)) + 
  geom_boxplot( width = 1) +
  ggtitle(paste0(snp, "-", gene, "\n", plot.main[i])) +
  scale_x_discrete(labels = paste0(names(var.counts), " (",var.counts, ")" ), name="") + 
  ylab("Proportions") + 
  theme(panel.grid.major=element_blank(), axis.text=element_text(size=13),axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(colour = gg_color_variants(length(var.counts)))) + 
  geom_vline(xintercept=c(1.5,2.5),color="white")


# ggsave(paste0(out.dir, "DM_prop_smp",plot.names[i],"_", snp, "-", gene, ".pdf"), width = 15, height = 7, units = "in", plot = ggp.prop.smp)


  ######## Plot Proportions
 
  prop <- dgeSQTL$fit[[snpGene]]$piH
  rownames(prop) <- rownames(expr)
  m.prop <- melt(prop)
  colnames(m.prop) <- c("Transcript", "Variant", "Proportion")
  m.prop$Transcript <- factor(m.prop$Transcript, levels = sort(as.character(unique(m.prop$Transcript))))


  ggp.prop <- ggplot(data = m.prop, aes(x = Transcript, y = Proportion, group = factor(Variant), colour = factor(Variant))) +
  theme_bw() +
    geom_line(size=1.5) +
    geom_point() +
  scale_color_manual(values=gg_color_variants(length(var.counts))) +
    theme(axis.text.x  = element_text(angle=80, vjust=0.5, size=12, colour = gg_color_hue(nlevels(m.prop$Transcript)), face="bold"), panel.background = element_blank(),  axis.line = element_line(colour = "grey")) 

  
  l <- mget(c("ggp.prop.smp", "ggp.prop")) 

  ggsave(paste0(out.dir, "DM_prop",plot.names[i],"_", snp, "-", gene, ".pdf"), width = 13, height = 10, units = "in", plot = do.call(arrangeGrob, l) )
  



}





























