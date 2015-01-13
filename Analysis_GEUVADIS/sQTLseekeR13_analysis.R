
# Created 9 Jan 2014
# Modyfied 9 Jan 2014


# Run sQTLseekeR analysis on 



##################################################################################
### load dgeSQTL produced with dmSNPsFiltering.R
##################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")

library(edgeR)
library(parallel)
library(dirmult)

Rdir <- "/home/gosia/R/R_Multinomial_project/DM_package/version5_sQTL/"

source("/home/gosia/R/R_Multinomial_project/DM_commonDisp/dirmult_code.R")
source(paste0(Rdir, "dmFunctions_v5.R"))

load("DM_genotypes/dgeSQTL.RData")


##################################################################################
# BioC 2.14

### Run sQTLseekeR 1.3 analysis on my filtering dmSNPsFiltering.R
# Too slow 
##################################################################################


############################################
### example
############################################
library(sQTLseekeR)

data(trans.exp.CEU.3genes)
data(genotype.CEU.3genes)

dim(trans.exp.CEU.3genes)
dim(genotype.CEU.3genes)

head(trans.exp.CEU.3genes)
head(genotype.CEU.3genes)

table( c(as.matrix(genotype.CEU.3genes[, -c(1,2)])) , useNA = "always")

genotype.CEU.3genes[2, 4] <- -1 ### does not work when missing values :/

# out <- sQTLseekeR(trans.exp.CEU.3genes, genotype.CEU.3genes, outfile.prefix=NULL, nb.perm.max=10000)
# head(out)

out <- sQTLseekeR(trans.exp = trans.exp.CEU.3genes, genotype = genotype.CEU.3genes, genes.to.test = NULL, outfile.prefix = "sQTLseekeR13_results", verbose = TRUE, min.nb.ext.scores = 1000, nb.perm.max = 3e+06, svQTL = FALSE, approx = TRUE, min.trans.exp = 0.01, min.gene.exp = 0.01, out.headers=TRUE)


out <- sQTLseekeR(trans.exp = trans.exp.CEU.3genes, genotype = genotype.CEU.3genes, genes.to.test = NULL, outfile.prefix = "sQTLseekeR13_resultsPermutations", verbose = TRUE, min.nb.ext.scores = 1000, nb.perm.max = 3e+06, svQTL = FALSE, approx = FALSE, min.trans.exp = 0.01, min.gene.exp = 0.01, out.headers=TRUE) ### veeery slow


out <- read.table("sQTLseekeR13_results.tsv", header = TRUE)

outP <- read.table("sQTLseekeR13_resultsPermutations.tsv", header = TRUE)



############################################
### prepare data
############################################


trans.exp <- cbind(dgeSQTL$genes[,c("ete_id", "gene_id")], do.call(rbind, dgeSQTL$counts))
colnames(trans.exp) <- c("TargetID", "Gene_Symbol", colnames(trans.exp)[-c(1,2)])

genotype <- cbind(dgeSQTL$SNPs, dgeSQTL$genotypes)  
colnames(genotype) <- c("gene", "ID", colnames(genotype)[-c(1,2)])

dim(genotype)

sum(grepl("snp", genotype$ID))


































