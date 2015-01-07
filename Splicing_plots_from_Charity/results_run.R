# 18 June 2014 (Last updated 23 June 2014)

# null simulations for drosophila
sim <- "drosophila_null"
path <- paste0("./simulation_", sim, "/")
source("results_loadfiles.R")
source("results_makeplots.R")
rm(list=ls(all=TRUE))

# null simulations for hsapiens
sim <- "hsapiens_null"
path <- paste0("./simulation_", sim, "/")
source("results_loadfiles.R")
source("results_makeplots.R")
rm(list=ls(all=TRUE))

# nonnull simulations for drosophila
sim <- "drosophila"
path <- paste0("./simulation_", sim, "/")
source("results_loadfiles.R")
source("results_makeplots.R")
rm(list=ls(all=TRUE))

# nonnull simulations for hsapiens
sim <- "hsapiens"
path <- paste0("./simulation_", sim, "/")
source("results_loadfiles.R")
source("results_makeplots.R")
rm(list=ls(all=TRUE))

# nonnull simulations for drosophila (V2)
sim <- "drosophila_V2"
path <- paste0("./simulation_", sim, "/")
source("results_loadfiles.R")
source("results_makeplots.R")
rm(list=ls(all=TRUE))

# null simulations for drosophila (V2)
sim <- "drosophila_V2_null"
path <- paste0("./simulation_", sim, "/")
source("results_loadfiles.R")
source("results_makeplots.R")
rm(list=ls(all=TRUE))

# nonnull simulations for hsapiens (V2)
sim <- "hsapiens_V2"
path <- paste0("./simulation_", sim, "/")
source("results_loadfiles.R")
source("results_makeplots.R")
rm(list=ls(all=TRUE))
# group specific filter
# sim <- "hsapiens_V2"
# path <- paste0("./simulation_", sim, "/")
# voom.exon <- read.table("~/voom_exon_grpspfilter_results.txt", header=TRUE)[,-3]
# colnames(voom.exon)[1:2] <- c("GeneID", "ExonID")
# voom.gene <- read.table("~/voom_gene_grpspfilter_results.txt", header=TRUE)[,-1]




# null simulations for hsapiens (V2)
sim <- "hsapiens_V2_null"
path <- paste0("./simulation_", sim, "/")
source("results_loadfiles.R")
source("results_makeplots.R")
rm(list=ls(all=TRUE))