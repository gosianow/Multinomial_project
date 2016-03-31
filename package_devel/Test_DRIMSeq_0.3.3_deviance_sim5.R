
# R32dev

library(devtools)

load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")

##############################################################################
# Simulation sim5
##############################################################################


library(DRIMSeq)
library(limma)
library(ggplot2)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_sim5'
simulation='hsapiens_node_nonull'
workers=10
count_method=c('htseq','kallisto','htseqprefiltered15','htseqprefiltered5','kallistofiltered5','kallistoprefiltered5')[2]
filter_method=c("filter0", "filter2")[1]
dispersion_common=TRUE
results_common=TRUE
disp_mode=c('grid','grid','optimize','optim','constrOptim')[1]
disp_moderation=c('none','common','none','none','none')[1]


##############################################################################

setwd(paste0(rwd, "/", simulation))
method_out <- "drimseq_0_3_3_f"

########################################################
# create metadata file
########################################################

metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata


##########################################################################
# 
##########################################################################

out_dir <- paste0(method_out, "/", count_method, "/", filter_method, "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



##########################################################################
# DRIMSeq results
##########################################################################


### Load object d

common_disp <- as.numeric(read.table(paste0(out_dir, "common_dispersion.txt")))
common_disp

disp <- "genewise"

if(disp_mode == "grid"){
  out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_", disp_moderation, "_")
}else{
  out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_")
}


load(paste0(out_name, "d.Rdata"))



### New out directory 
out_dir <- paste0(method_out, "/", count_method, "/", filter_method, "/", "test_deviance/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if(disp_mode == "grid"){
  out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_", disp_moderation, "_")
}else{
  out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_")
}




##########################################################################
# Use different tests
##########################################################################


## LR test
dlr <- dmTest(d, test = "lr", BPPARAM = BiocParallel::MulticoreParam(workers = workers))

reslr <- results(dlr)

plotTest(dlr, out_dir = paste0(out_name, "lr_"))


## Fql test
dfql <- dmTest(d, test = "fql", BPPARAM = BiocParallel::MulticoreParam(workers = workers))

resfql <- results(dfql)

plotTest(dfql, out_dir = paste0(out_name, "fql_"))


## F test
df <- dmTest(d, test = "f", BPPARAM = BiocParallel::MulticoreParam(workers = workers))

resf <- results(df)

plotTest(df, out_dir = paste0(out_name, "f_"))



##########################################################################
### Scatter plots of p-values

resm <- merge(reslr[, c("gene_id", "pvalue")], resf[, c("gene_id", "pvalue")], by = c("gene_id"), suffixes = c(".lr",".f"))

pdf(paste0(out_name, "scatter_pvalues_f.pdf"))

plot(resm[, "pvalue.lr"], resm[, "pvalue.f"], xlim = c(0, 1), ylim = c(0, 1), xlab = "LR", ylab = "F")
abline(a = 0, b = 1, col = "blue")

dev.off()




resm <- merge(reslr[, c("gene_id", "pvalue")], resfql[, c("gene_id", "pvalue")], by = c("gene_id"), suffixes = c(".lr",".fql"))

pdf(paste0(out_name, "scatter_pvalues_fql.pdf"))

plot(resm[, "pvalue.lr"], resm[, "pvalue.fql"], xlim = c(0, 1), ylim = c(0, 1), xlab = "LR", ylab = "Fql")
abline(a = 0, b = 1, col = "blue")

dev.off()


##########################################################################
### Recalculate deviance for negative cases

stats <- statistics(d)


stats <- stats[order(stats$dev_null, decreasing = FALSE), ]
head(stats)

gene <- stats[2, "gene_id"]


plotFit(d, gene_id = gene, out_dir = paste0(out_name, "d", "_neg_dev_"))


pi <- d@fit_null[[gene]][, "null"]
pi
gamma0 <- d@genewise_dispersion[gene]
gamma0
y <- d@counts[[gene]]
y



DRIMSeq:::dm_devG(pi = pi[-length(pi)], gamma0 = gamma0, y = y)

### Save an example for Mark
save("pi", "gamma0", "y", file = "/home/gosia/deviance.Rdata")


ll_mod <- sum(lgamma(y + pi * gamma0) - lgamma(pi * gamma0) , na.rm = TRUE )

pi_sat <- y/matrix(colSums(y), nrow(y), ncol(y), byrow = TRUE)

ll_sat <- sum(lgamma(y + pi_sat * gamma0) - lgamma(pi_sat * gamma0) , na.rm = TRUE) # Inf for y = 0

D <- 2 * (ll_sat - ll_mod)







### Check whether the LR statistics are equal

statsf <- statistics(df)

dev_full <- statsf[, c("dev_C1", "dev_C2")]
dev_null <- statsf[, c("dev_null")]

lrd <- dev_null - rowSums(dev_full)

all.equal(lrd, results(dlr)$lr)


##########################################################################
### Find when negative test statistics occure


## negative LR when transcripts have the same proportions in null and full models 

table(reslr$lr < 0) 
table(resf$f < 0)

genes <- reslr[which(reslr$lr < 0), "gene_id"]

gene <- "ENSG00000198003"

dlr@counts[[gene]]


dlr@fit_full[[gene]]
dlr@fit_null[[gene]]






## negative Fql

statsf <- statistics(df)

table(resfql$f < 0)

table((statsf$dev_C1 + statsf$dev_C2) < 0)

table(statsf$dev_C1 < 0)

table(statsf$dev_C2 < 0)

zero_genes <- sapply(1:length(d@counts), function(g){
  
  any(d@counts[[g]] == 0)
  
})

names(zero_genes) <- names(d@counts)



neg_genes <- statsf[statsf$dev_null < 0, "gene_id"]

table(zero_genes[neg_genes])



statsf <- statsf[order(statsf$dev_null, decreasing = FALSE), ]

rownames(statsf) <- statsf$gene_id


genes <- names(which(zero_genes[neg_genes] == FALSE))
gene <- genes[1]



plotFit(dfql, gene_id = gene, out_dir = paste0(out_name, "dfql", "_neg_dev_no_zeros_1_"))

plotFit(dlr, gene_id = gene, out_dir = paste0(out_name, "dlr", "_neg_dev_no_zeros_1_"))


for(i in 1:length(genes)){
  
  plotFit(dfql, gene_id = genes[i], out_dir = paste0(out_name, "dfql", "_neg_dev_no_zeros_", i, "_"))
  
}



##########################################################################
### Do fitting just for 1 gene

genenr <- 9

gene <- genes[genenr] ## "ENSG00000136026"
gene


dsub <- dfql[gene, ]
disp_init <- common_dispersion(dfql)
disp_init


dsub <- dmDispersion(dsub, mean_expression = TRUE,
  common_dispersion = FALSE, genewise_dispersion = TRUE,
  disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
  disp_tol = 1e-12, disp_init = disp_init, disp_init_weirMoM = TRUE,
  disp_grid_length = 21, disp_grid_range = c(-10, 10),
  disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

genewise_dispersion(dsub)


dsub <- dmFit(dsub, dispersion = "genewise_dispersion",
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotFit(dsub, gene_id = gene, out_dir = paste0(out_name, "dsub_genewise", "_neg_dev_no_zeros_", genenr ,"_"))



## Set dispersion to 100
common_dispersion(dsub) <- 100

dsub <- dmFit(dsub, dispersion = "common_dispersion",
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotFit(dsub, gene_id = gene, out_dir = paste0(out_name, "dsub_set", "_neg_dev_no_zeros_", genenr ,"_"))



## Do not use CR adjustement
dsub <- dmDispersion(dsub, mean_expression = TRUE,
  common_dispersion = FALSE, genewise_dispersion = TRUE,
  disp_adjust = FALSE, disp_mode = "grid", disp_interval = c(0, 1e+05),
  disp_tol = 1e-12, disp_init = disp_init, disp_init_weirMoM = TRUE,
  disp_grid_length = 21, disp_grid_range = c(-10, 10),
  disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

genewise_dispersion(dsub)


dsub <- dmFit(dsub, dispersion = "genewise_dispersion",
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotFit(dsub, gene_id = gene, out_dir = paste0(out_name, "dsub_noadj", "_neg_dev_no_zeros_", genenr ,"_"))








# ### Run version with minmax constrains on proportions
# dsub <- dmFit(dsub, dispersion = "genewise_dispersion",
#   prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
#   BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# plotFit(dsub, gene_id = gene, out_dir = paste0(out_name, "dsub_minmax", "_neg_dev_no_zeros_", genenr ,"_"))






pi <- d@fit_null[[gene]][, "null"]
pi
gamma0 <- d@genewise_dispersion[gene]
gamma0
y <- d@counts[[gene]]
y

statsf[gene, ]

DRIMSeq:::dm_devG(pi = pi[-length(pi)], gamma0 = gamma0, y = y)


DRIMSeq:::dm_devG(pi = pi[-length(pi)], gamma0 = 100, y = y)




##########################################################################
### Estimate dispersion per group 


disp_init <- common_dispersion(dfql)
disp_init

dsub1 <- dfql[gene, 1:3]

dsub1 <- dmDispersion(dsub1, mean_expression = TRUE,
  common_dispersion = FALSE, genewise_dispersion = TRUE,
  disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
  disp_tol = 1e-12, disp_init = disp_init, disp_init_weirMoM = TRUE,
  disp_grid_length = 21, disp_grid_range = c(-10, 10),
  disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

genewise_dispersion(dsub1)


dsub1 <- dmFit(dsub1, dispersion = "genewise_dispersion",
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotFit(dsub1, gene_id = gene, out_dir = paste0(out_name, "dsub1", "_neg_dev_no_zeros_", genenr ,"_"))



dsub2 <- dfql[gene, 4:6]
samples(dsub2)


dsub2 <- dmDispersion(dsub2, mean_expression = TRUE,
  common_dispersion = FALSE, genewise_dispersion = TRUE,
  disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
  disp_tol = 1e-12, disp_init = disp_init, disp_init_weirMoM = TRUE,
  disp_grid_length = 21, disp_grid_range = c(-10, 10),
  disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

genewise_dispersion(dsub2)


dsub2 <- dmFit(dsub2, dispersion = "genewise_dispersion",
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotFit(dsub2, gene_id = gene, out_dir = paste0(out_name, "dsub2", "_neg_dev_no_zeros_", genenr ,"_"))






##########################################################################
### Estimate proportions with smaller tol_prop


d2 <- dmDispersion(d, mean_expression = TRUE,
  common_dispersion = TRUE, genewise_dispersion = TRUE,
  disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
  disp_tol = 1e-12, disp_init = 100, disp_init_weirMoM = TRUE,
  disp_grid_length = 21, disp_grid_range = c(-10, 10),
  disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotDispersion(d2, out_dir = paste0(out_name, "smaller_tol_"))

d2 <- dmFit(d2, dispersion = "genewise_dispersion",
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))


plotFit(d2, gene_id = gene, out_dir = paste0(out_name, "d2", "_neg_dev_no_zeros_", genenr ,"_"))


load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")


##########################################################################
### Estimation with new minmax constrain in constrOptim

disp_init <- common_dispersion(d)

d3 <- d[88, ]

d3 <- dmDispersion(d, mean_expression = TRUE,
  common_dispersion = FALSE, genewise_dispersion = TRUE,
  disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
  disp_tol = 1e-08, disp_init = disp_init, disp_init_weirMoM = TRUE,
  disp_grid_length = 21, disp_grid_range = c(-10, 10),
  disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
  prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0,
  BPPARAM = BiocParallel::SerialParam())



d3 <- dmDispersion(d, mean_expression = TRUE,
  common_dispersion = FALSE, genewise_dispersion = TRUE,
  disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
  disp_tol = 1e-08, disp_init = disp_init, disp_init_weirMoM = TRUE,
  disp_grid_length = 21, disp_grid_range = c(-10, 10),
  disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
  prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))


plotDispersion(d3, out_dir = paste0(out_name, "ci_minmax_"))

d3 <- dmFit(d3, dispersion = "genewise_dispersion",
  prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))


plotFit(d3, gene_id = gene, out_dir = paste0(out_name, "d3", "_neg_dev_no_zeros_", genenr ,"_"))






##########################################################################
### Workflow from the original DRIMSeq analysis


# # genewise dispersion
# d <- dmDispersion(d, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_mode = disp_mode, disp_init = common_disp, disp_moderation = disp_moderation, disp_prior_df = 0.1, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# common_dispersion(d) <- common_disp
# 
# plotDispersion(d, out_dir = out_name)
# 
# d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# ## LR test
# d <- dmTest(d, test = "lr", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# plotTest(d, out_dir = paste0(out_name, "lr_"))
# 
# res <- results(d)
# 
# write.table(res, paste0(out_name, "results_lr.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# ## F test
# d <- dmTest(d, test = "f", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# plotTest(d, out_dir = paste0(out_name, "f_"))
# 
# res <- results(d)
# 
# write.table(res, paste0(out_name, "results_f.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# save(d, file = paste0(out_name, "d.Rdata"))







##########################################################################
### Plot of how using F statistics gives heavy tails

pdf("/home/gosia/F.pdf")

x <- seq(1, 20, 0.1)

plot(x, 4 * dchisq(x * 4, df = 4), type = "l", lwd = 2)
lines(x, df(x, df1 = 4, df2 = 3), type = "l", col = "blue", lwd = 2)
lines(x, df(x, df1 = 4, df2 = 10), type = "l", col = "orange", lwd = 2)
lines(x, df(x, df1 = 4, df2 = 50), type = "l", col = "salmon", lwd = 2)
legend("topright", legend = c("4X ~ Chi^2(4)", "X ~ F(4, 3)", "X ~ F(4, 10)", "X ~ F(4, 50)"), lty = 1, col = c("black", "blue", "orange", "salmon"))


dev.off()












