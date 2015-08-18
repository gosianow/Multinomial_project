# BioC 3.1

# Created 3 June 2015

##############################################################################################################
# Test DM package
##############################################################################################################

library(DM)


data <- dataDS_dmDSdata
dmDSplotData(data)


data <- dmDSfilter(data)
dmDSplotData(data)


data <- dmDSdispersion(data)


dmDSplotDispersion(data)


data_fit <- dmDSfit(data, dispersion = "tagwise_dispersion", prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 1))



dmDSplotFit(data_fit, gene_id = "FBgn0001316", plot_type = "barplot", order = TRUE, plot_full = FALSE, plot_nunll = FALSE)



table <- dmDStest(data_fit)


dmDSplotTest(table)









