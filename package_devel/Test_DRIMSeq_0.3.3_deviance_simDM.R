
# R32dev

library(devtools)

load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")

##############################################################################
# Simulations from DM - proportions_f
##############################################################################


library(BiocParallel)
library(pryr)
library(plyr)
library(dirmult)
library(limma)
library(DRIMSeq)
library(ggplot2)
library(reshape2)
library(tools)


##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_dm/drimseq/'
simulation_script='/home/gosia/R/drimseq_paper/simulations_dm/dm_simulate.R'
workers=8
sim_name=''
run='run1'
m=1000
n=3 # Number of samples
nm=1000
nd=0
nr_features=c(3,5,7,10,13,15,17)
param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters_drimseq_0_3_3/kim_kallisto/disp_common_kim_kallisto.txt'
save=TRUE
out_suffix='proportions_decay'

##############################################################################
# Read in the arguments
##############################################################################


source(simulation_script)

### Dispersion
params <- read.table(param_gamma_path, header = FALSE, sep = "\t")

sim_disp_common <- FALSE
sim_disp_genewise <- FALSE

# Common dispersion
if(ncol(params) == 1){
  g0 <- as.numeric(params)
  print(g0)
  sim_disp_common <- TRUE
}

# Genewise dispersion from lognormal distribution
if(ncol(params) == 2){
  g0_meanlog <- params[1, 2]
  g0_sdlog <- params[2, 2]
  print(g0_meanlog)
  print(g0_sdlog)
  sim_disp_genewise <- TRUE
}



##############################################################################

dir.create(rwd, recursive = T, showWarnings = FALSE)
setwd(rwd)

out_dir <- "proportions_f/run/"
dir.create(out_dir, recursive = T, showWarnings = FALSE)

out_name <- paste0(sim_name, "n", n, "_nm", nm, "_nd", nd, "_",  basename(file_path_sans_ext(param_gamma_path)), "_")

out_name

if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}


plot_dir <- "proportions_f/test_deviance/"
dir.create(plot_dir, recursive = T, showWarnings = FALSE)

plot_name <- paste0(plot_dir, out_name)


##############################################################################
### Load simulation results 
##############################################################################



for(j in 1:length(nr_features)){
  # j = 5
  print(nr_features[j])
  
  out_name <- paste0(sim_name, "n", n, "_nm", nm, "_nd", nd, "_",  basename(file_path_sans_ext(param_gamma_path)), "_")
  
  load(paste0(out_dir, out_name, "d_", out_suffix, "_", run, "_", nr_features[j], "_moderation_none",".Rdata"))
  
  out_name <- paste0(plot_dir, out_name, nr_features[j], "_")
  
  
  
  
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
  
  
  
  statsf <- statistics(df)
  
  neg_genes <- statsf[statsf$dev_null < 0, "gene_id"]
  
  zero_genes <- sapply(1:length(d@counts), function(g){
    
    any(d@counts[[g]] == 0)
    
  })
  names(zero_genes) <- names(d@counts)
  
  
  
  genes <- names(which(zero_genes[neg_genes] == FALSE))
  length(genes)
  
  
  genes <- neg_genes
  length(genes)
  
  gene <- genes[1]
  
  
  for(i in 1:10){
    
    plotFit(dfql, gene_id = genes[i], out_dir = paste0(out_name, "dfql", "_neg_dev_no_zeros_", i, "_"))
    
  }
  
  
  
  disp_init <- common_dispersion(d)
  
  d3 <- dmDispersion(d, mean_expression = TRUE,
                     common_dispersion = FALSE, genewise_dispersion = TRUE,
                     disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
                     disp_tol = 1e-08, disp_init = disp_init, disp_init_weirMoM = TRUE,
                     disp_grid_length = 21, disp_grid_range = c(-10, 10),
                     disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
                     prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0,
                     BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  
  load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")
  
  d3 <- d[88, ]
  
  d3 <- dmDispersion(d3, mean_expression = TRUE,
                     common_dispersion = FALSE, genewise_dispersion = TRUE,
                     disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
                     disp_tol = 1e-08, disp_init = disp_init, disp_init_weirMoM = TRUE,
                     disp_grid_length = 21, disp_grid_range = c(-10, 10),
                     disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
                     prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0,
                     BPPARAM = BiocParallel::SerialParam())
  
  
  
  plotDispersion(d3, out_dir = paste0(out_name, "ci_minmax_"))
  
  d3 <- dmFit(d3, dispersion = "genewise_dispersion",
              prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0,
              BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  

}


























