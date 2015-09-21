##############################################################################

# BioC 3.1
# Created 10 Sep 2015:

# Run DM_0.2.3 on 
# - htseq couts
# - kallisto counts

# Updated 11 Sep 2015
# - Use htseq counts from INCOMPLETE_KALLISTOEST

##############################################################################################################

### run DM with different modes 

##############################################################################################################

setwd("/home/gosia/multinomial_project/simulations_sim5_drosophila_noDE_noNull/")

# setwd("/home/gosia/multinomial_project/simulations_sim5_hsapiens_noDE_noNull/")


library(ggplot2)
library(DM)




count_method_list <- c("htseq", "kallisto", "htseq_kallisto_filter")
count_method <- count_method_list[3]

filter_method_list <- c("filtering_dexseq", "filtering_min3prop0_01min6cpm1maxInf") 
filter_method <- filter_method_list[1]

out_dir <- paste0("dm_0_2_3/", count_method, "/", filter_method, "/")



##### run DM pipelines : common_dispersion

load(paste0(out_dir, "d.Rdata"))

disp <- "common_dispersion"
out_name <- paste0(out_dir, "/DM_", disp, "_")

d <- dmDispersion(d, genewise_dispersion = FALSE, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

d <- dmFit(d, dispersion = disp, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

d <- dmLRT(d, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

plotLRT(d, out_dir = out_name)

results <- DM::results(d)

write.table(results, paste0(out_name, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
save(d, file = paste0(out_name, "d.Rdata"))


common_dispersion <- DM::common_dispersion(d)
write.table(common_dispersion, paste0(out_dir, "common_dispersion.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)




##### run DM pipelines : genewise_dispersion

disp_mode_list = c("grid", "grid", "grid", "grid", "grid")
disp_moderation_list <- c("none", "common", "trended", "common", "trended")
disp_prior_df_list <- c(10, 10, 10, 4, 4)


for(i in 2:3){
  # i = 3
  print(i)
  
  load(paste0(out_dir, "d.Rdata"))
  common_dispersion <- as.numeric(read.table(paste0(out_dir, "common_dispersion.txt")))
  
  disp <- "genewise_dispersion"
  disp_mode <- disp_mode_list[i]
  disp_moderation <- disp_moderation_list[i]
  disp_prior_df <- disp_prior_df_list[i]
  
  
  if(disp_mode == "grid")
    out_name <- paste0(out_dir, "/DM_", disp, "_", disp_mode, "_", disp_moderation, "_")
  
  if(disp_prior_df != 10)
    out_name <- paste0(out_name, disp_prior_df, "_")
  
  
  d <- dmDispersion(d, common_dispersion = FALSE, genewise_dispersion = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 5))
  common_dispersion(d) <- common_dispersion
  
  d <- dmDispersion(d, genewise_dispersion = TRUE, disp_mode = disp_mode, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 5))
  
  plotDispersion(d, out_dir = out_name)
  
  d <- dmFit(d, dispersion = disp, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

  d <- dmLRT(d, BPPARAM = BiocParallel::MulticoreParam(workers = 5))

  plotLRT(d, out_dir = out_name)

  results <- DM::results(d)

  write.table(results, paste0(out_name, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  save(d, file = paste0(out_name, "d.Rdata"))

}



























