##############################################################################
# Fix the issue with trended dispersion by not considering the genes with 
# dispersion that is on the border of the grid
# Simulation sim5

# Created on 24 Mar 2016
# Modiefied on 30 Mar 2016
##############################################################################

# R32

library(DRIMSeq)
library(limma)
library(ggplot2)
library(reshape2)

library(devtools)
load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")

##############################################################################
# Test arguments
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_sim5'
simulation='drosophila_node_nonull'
workers=10
count_method=c('htseq','kallisto','htseqprefiltered15','htseqprefiltered5','kallistofiltered5','kallistoprefiltered5')[1]
filter_method=c("filter0", "filter2")[1]
dispersion_common=TRUE
results_common=TRUE
disp_mode=c('grid','grid','optimize','optim','constrOptim')[1]
disp_moderation=c('none','common','none','none','none')[1]
lr_contribution_function_path <- "/home/gosia/R/drimseq_paper/help_functions/dmDS_lr_contribution.R"

##############################################################################

setwd(paste0(rwd, "/", simulation))
method_out <- "drimseq_0_3_3"

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
out_dir <- paste0(method_out, "/", count_method, "/", filter_method, "/", "trended_dispersion/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if(disp_mode == "grid"){
  out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_", disp_moderation, "_")
}else{
  out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_")
}



res <- results(d)

plotTest(d, out_dir = paste0(out_name))
plotDispersion(d, out_dir = paste0(out_name))


#######################################################
# load simulation info
#######################################################


simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, as.is = TRUE, sep = "\t")

truth_file <- list.files("3_truth/", pattern = "truth")
truth_file

truth <- read.table(paste0("3_truth/", truth_file), header = TRUE, as.is = TRUE, sep = "\t")

truth <- truth[, c("gene", "ds_status", "de_status", "TPM", "nbr_isoforms", "diff_IsoPct", "nbrexonbins")]
rownames(truth) <- truth$gene

colnames(truth)[1] <- "gene_id"


##########################################################################
### merge results with truth


resm <- merge(res, truth, by = c("gene_id"), sort = FALSE)

resm <- unique(resm)


##########################################################################
### Plot dispersion versus mean with marked FP ----

common_disp <- common_dispersion(d)

ggp <- plotDispersion(d)

df <- ggp$data

res_fp <- resm[resm$ds_status == 0 & resm$adj_pvalue < 0.05, ]

df_fp <- df[res_fp[, "gene_id"], ]

ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion")) +
  theme_bw() +
  xlab("Log10 of mean expression") +
  ylab("Log10 of gamma_+") +
  geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") + 
  geom_hline(yintercept = log10(common_disp), colour = "black", linetype = "dashed", size =  0.5) +
  geom_point(data = df_fp, aes_string(x = "mean_expression", y = "dispersion"), color = "orange")


pdf(paste0(out_name, "dispersion_vs_mean_fp.pdf"))
print(ggp)
dev.off()


##########################################################################
### Plot dispersion versus mean with marked TP ####

common_disp <- common_dispersion(d)

ggp <- plotDispersion(d)

df <- ggp$data

res_tp <- resm[resm$ds_status == 1 & resm$adj_pvalue < 0.05, ]

df_tp <- df[res_tp[, "gene_id"], ]

ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion")) +
  theme_bw() +
  xlab("Log10 of mean expression") +
  ylab("Log10 of gamma_+") +
  geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") + 
  geom_hline(yintercept = log10(common_disp), colour = "black", linetype = "dashed", size =  0.5) +
  geom_point(data = df_tp, aes_string(x = "mean_expression", y = "dispersion"), color = "blue")


pdf(paste0(out_name, "dispersion_vs_mean_tp.pdf"))
print(ggp)
dev.off()



##########################################################################
### Plot dispersion versus mean with marked P ####

common_disp <- common_dispersion(d)

ggp <- plotDispersion(d)

df <- ggp$data

res_p <- resm[resm$ds_status == 1, ]

df_p <- df[res_p[, "gene_id"], ]

ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion")) +
  theme_bw() +
  xlab("Log10 of mean expression") +
  ylab("Log10 of gamma_+") +
  geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") + 
  geom_hline(yintercept = log10(common_disp), colour = "black", linetype = "dashed", size =  0.5) +
  geom_point(data = df_p, aes_string(x = "mean_expression", y = "dispersion"), color = "dodgerblue")


pdf(paste0(out_name, "dispersion_vs_mean_p.pdf"))
print(ggp)
dev.off()







##########################################################################
### Fix the trended dispersion method by not taking into account the outliers

x <- d

disp_init <- common_dispersion(x)

counts = x@counts; samples = x@samples; mean_expression = x@mean_expression; disp_adjust = TRUE; disp_mode = c("optimize", "optim", "constrOptim", "grid")[4]; disp_interval = c(0, 1e+5); disp_tol = 1e-08; disp_init = disp_init; disp_init_weirMoM = TRUE; disp_grid_length = 21; disp_grid_range = c(-10, 10); disp_moderation = c("none", "common", "trended")[2]; disp_prior_df = 1; disp_span = 0.2; prop_mode = c( "constrOptim", "constrOptimG", "FisherScoring")[2]; prop_tol = 1e-12; verbose = FALSE; BPPARAM = BiocParallel::MulticoreParam(workers = 10)



inds <- 1:length(counts)

group <- samples$group
ngroups <- nlevels(group)
lgroups <- levels(group)
nlibs <- length(group)

igroups <- lapply(lgroups, function(gr){which(group == gr)})
names(igroups) <- lgroups




#####################################
### Plot dispersion versus mean with marked grid points

splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)
splineDisp <- disp_init * 2^splinePts

common_disp <- common_dispersion(d)

ggp <- plotDispersion(d)

df <- ggp$data

ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion")) +
  theme_bw() +
  xlab("Log10 of mean expression") +
  ylab("Log10 of gamma_+") +
  geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") + 
  geom_hline(yintercept = log10(common_disp), colour = "black", linetype = "dashed", size =  0.5) +
  geom_hline(yintercept = log10(splineDisp), colour = "orange", linetype = "dashed", size =  0.3)

pdf(paste0(out_name, "dispersion_vs_mean_grid.pdf"))
print(ggp)
dev.off()





#####################################
### Plot dispersion versus mean with marked grid points


# ### Twice more dense grid toward the common dispersion
# splinePts_uni <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)
# 
# signPts <- sign(splinePts_uni)
# 
# splinePts <- 2^abs(splinePts_uni) * signPts
# 
# splinePts <- min(splinePts_uni) + (splinePts - min(splinePts)) * (max(splinePts_uni) - min(splinePts_uni)) / (max(splinePts) - min(splinePts))



# ## Half of the grid with twice density
# splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = ceiling(4/3 * disp_grid_length))
# 
# grid_exclude_length <- floor((ceiling(4/3 * disp_grid_length) - disp_grid_length) / 2)
# exclude_left <- (1:grid_exclude_length) * 2
# exclude_right <- length(splinePts) + 1 - (1:grid_exclude_length) * 2
# 
# splinePts <- splinePts[-c(exclude_left, exclude_right)]



### More dense grid toward the common dispersion 
splinePts_uni <- sort(unique(c(0, seq(from = -10, to = 10, length = 11))))

nr_positive_splitting <- sum(sign(splinePts_uni) == 1)
nr_negative_splitting <- sum(sign(splinePts_uni) == -1)

max_splitting <- max(nr_positive_splitting, nr_negative_splitting)
min_splitting <- min(nr_positive_splitting, nr_negative_splitting)

if(nr_positive_splitting == max_splitting)
  nr_splitting <- c((max_splitting - min_splitting + 1):max_splitting, max_splitting:1) + 2

if(nr_negative_splitting == max_splitting)
  nr_splitting <- c(1:max_splitting, max_splitting:(max_splitting - min_splitting + 1)) + 2



splinePts <- lapply(1:(length(splinePts_uni) - 1), function(i){
  
  seq(from = splinePts_uni[i], to = splinePts_uni[i + 1], length = nr_splitting[i])
  
})

splinePts <- sort(unique(unlist(splinePts)))

disp_grid_length <- length(splinePts)


### Standard
# splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)


splineDisp <- disp_init * 2^splinePts


common_disp <- common_dispersion(d)

ggp <- plotDispersion(d)

df <- ggp$data

ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion")) +
  theme_bw() +
  xlab("Log10 of mean expression") +
  ylab("Log10 of gamma_+") +
  geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") + 
  geom_hline(yintercept = log10(common_disp), colour = "black", linetype = "dashed", size =  0.5) +
  geom_hline(yintercept = log10(splineDisp), colour = "orange", linetype = "dashed", size =  0.3)

pdf(paste0(out_name, "dispersion_vs_mean_grid_test1.pdf"))
print(ggp)
dev.off()




#####################################
### Continue with fixing the trended moderation

splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)
splineDisp <- disp_init * 2^splinePts


### calculate the likelihood for each gene at the spline dispersion points
seq_disp_grid_length <- seq(disp_grid_length)


loglikL <- BiocParallel::bplapply(inds, function(g){
  # g = 1237
  # print(g)
  
  ll <- numeric(disp_grid_length)
  
  for(i in seq_disp_grid_length){
    # i = 1
    
    out <- dm_profileLikTagwise(gamma0 = splineDisp[i], y = counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
    
    if(is.na(out)){
      ll <- rep(NA, disp_grid_length)
      break
    }
    
    ll[i] <- out
    
  }
  
  return(ll)
  
}, BPPARAM = BPPARAM)


loglik <- do.call(rbind, loglikL)

### See how many NAs for each gene
nas <- is.na(loglik)
table(rowSums(nas))





not_nas <- complete.cases(loglik)        


### Plot mean expression for genes with NAs 
ggdf <- data.frame(mean_expression = mean_expression, not_nas = not_nas)

ggp <- ggplot(ggdf, aes(x = log10(mean_expression), fill = not_nas)) +
  geom_density(alpha = 0.6)


pdf(paste0(out_name, "mean_expression_nas.pdf"))
print(ggp)
dev.off()






loglik <- loglik[not_nas, , drop = FALSE]
mean_expression <- mean_expression[not_nas]




#####################################
### Plot profile likelihoods on grid for interesting genes
ggdf <- melt(data.frame(loglik, gene_id = names(counts)[not_nas]), id.vars = "gene_id", variable.name = "splineId", value.name = "profile_loglik", stringsAsFactors = FALSE)

spline_match <- data.frame(splineDisp = splineDisp, splineId = paste0("X", 1:length(splineDisp)), stringsAsFactors = FALSE)

splinem <- match(ggdf$splineId, spline_match$splineId)

ggdf$splineDisp <- spline_match$splineDisp[splinem]




### Gene with boudry dispersion 
gene <- "ENSG00000005059"

ggp <- ggplot(ggdf[ggdf$gene_id == gene, ], aes(x = log10(splineDisp), y = profile_loglik)) +
  geom_line() + 
  geom_vline(xintercept = log10(splineDisp), colour = "orange", linetype = "dashed", size =  0.3) +
  geom_vline(xintercept = log10(disp_init), colour = "red", linetype = "dashed", size =  0.3)


pdf(paste0(out_name, "profile_loglik_", gene,".pdf"))
print(ggp)
dev.off()

ggp <- ggplot(ggdf[ggdf$gene_id == gene, ], aes(x = splineDisp, y = profile_loglik)) +
  geom_line() + 
  geom_vline(xintercept = splineDisp, colour = "orange", linetype = "dashed", size =  0.3) +
  geom_vline(xintercept = disp_init, colour = "red", linetype = "dashed", size =  0.3)

pdf(paste0(out_name, "profile_loglik_", gene,"_notlog.pdf"))
print(ggp)
dev.off()

plotFit(d, gene_id = gene, out_dir = paste0(out_name, "profile_loglik_"))



### Gene with "good" dispersion
gene <- "ENSG00000000003"

ggp <- ggplot(ggdf[ggdf$gene_id == gene, ], aes(x = log10(splineDisp), y = profile_loglik)) +
  geom_line() + 
  geom_vline(xintercept = log10(splineDisp), colour = "orange", linetype = "dashed", size =  0.3) +
  geom_vline(xintercept = log10(disp_init), colour = "red", linetype = "dashed", size =  0.3)

pdf(paste0(out_name, "profile_loglik_", gene,".pdf"))
print(ggp)
dev.off()

ggp <- ggplot(ggdf[ggdf$gene_id == gene, ], aes(x = splineDisp, y = profile_loglik)) +
  geom_line() + 
  geom_vline(xintercept = splineDisp, colour = "orange", linetype = "dashed", size =  0.3) +
  geom_vline(xintercept = disp_init, colour = "red", linetype = "dashed", size =  0.3)

pdf(paste0(out_name, "profile_loglik_", gene,"_notlog.pdf"))
print(ggp)
dev.off()


ggp <- ggplot(ggdf[ggdf$gene_id == gene, ], aes(x = splineDisp, y = profile_loglik)) +
  geom_line() + 
  geom_vline(xintercept = splineDisp, colour = "orange", linetype = "dashed", size =  0.3) +
  geom_vline(xintercept = disp_init, colour = "red", linetype = "dashed", size =  0.3) +
  coord_cartesian(xlim = c(0, 5000))

pdf(paste0(out_name, "profile_loglik_", gene,"_notlog_xlim.pdf"))
print(ggp)
dev.off()


plotFit(d, gene_id = gene, out_dir = paste0(out_name, "profile_loglik_"))




### Gene with very high dispersion = low concentration gamma_+

gene <- "ENSG00000205835"

ggp <- ggplot(ggdf[ggdf$gene_id == gene, ], aes(x = log10(splineDisp), y = profile_loglik)) +
  geom_line() + 
  geom_vline(xintercept = log10(splineDisp), colour = "orange", linetype = "dashed", size =  0.3) +
  geom_vline(xintercept = log10(disp_init), colour = "red", linetype = "dashed", size =  0.3)

pdf(paste0(out_name, "profile_loglik_", gene,".pdf"))
print(ggp)
dev.off()

ggp <- ggplot(ggdf[ggdf$gene_id == gene, ], aes(x = splineDisp, y = profile_loglik)) +
  geom_line() + 
  geom_vline(xintercept = splineDisp, colour = "orange", linetype = "dashed", size =  0.3) +
  geom_vline(xintercept = disp_init, colour = "red", linetype = "dashed", size =  0.3)

pdf(paste0(out_name, "profile_loglik_", gene,"_notlog.pdf"))
print(ggp)
dev.off()


ggp <- ggplot(ggdf[ggdf$gene_id == gene, ], aes(x = splineDisp, y = profile_loglik)) +
  geom_line() + 
  geom_vline(xintercept = splineDisp, colour = "orange", linetype = "dashed", size =  0.3) +
  geom_vline(xintercept = disp_init, colour = "red", linetype = "dashed", size =  0.3) +
  coord_cartesian(xlim = c(0, 5000))

pdf(paste0(out_name, "profile_loglik_", gene,"_notlog_xlim.pdf"))
print(ggp)
dev.off()


plotFit(d, gene_id = gene, out_dir = paste0(out_name, "profile_loglik_"))



#####################################


if(nrow(loglik) == 0){
  dispersion <- rep(NA, length(inds))
  names(dispersion) <- names(counts)
  if(verbose) cat("*** Genewise dispersion: ", head(dispersion), "... \n")
  return(dispersion)
}




### Check where the grid is maximized 
grid_max <- apply(loglik, 1, which.max)
table(grid_max)


### In the calculation of moderation, do not take into account genes that have dispersion on the top boundry of the grid (3 last grid points)
not_boundry <- grid_max < (disp_grid_length - 3)
table(not_boundry)

not_boundry[7351] <- FALSE


if(disp_moderation != "none"){
  
  # nlibs <- length(group)
  # priorN <- disp_prior_df/(nlibs - ngroups) ### analogy to edgeR
  
  priorN <- disp_prior_df
  
  switch(disp_moderation, 
    
    common = {
      
      moderation <- colMeans(loglik[not_boundry, , drop = FALSE])
      
      loglik <- sweep(loglik, 2, priorN * moderation, FUN = "+")
      
    },
    
    trended = {
      
      ### Use non boundry genes for calculating the moderation
      mean_expression_not_boundry <- mean_expression[not_boundry]
      loglik_not_boundry <- loglik[not_boundry, , drop = FALSE]
      
      o <- order(mean_expression_not_boundry)
      oo <- order(o)
      
      width <- floor(disp_span * nrow(loglik_not_boundry))
      
      moderation_not_boundry <- edgeR::movingAverageByCol(loglik_not_boundry[o, , drop = FALSE], width = width)[oo, , drop = FALSE]
      
      ### Fill in moderation values for the boundy genes
      moderation <- matrix(NA, nrow = nrow(loglik), ncol = ncol(loglik))
      
      moderation[not_boundry, ] <- moderation_not_boundry
      
      o <- order(mean_expression)
      oo <- order(o)
      
      moderation <- moderation[o, , drop = FALSE]
      not_boundry <- not_boundry[o]
      
      ### Last value in not_boundry must be TRUE
      if(not_boundry[length(not_boundry)] == FALSE){
        
        last_true <- max(which(not_boundry))
        moderation[length(not_boundry), ] <- moderation[last_true, ]
        
        not_boundry[length(not_boundry)] <- TRUE

      }
      
      not_boundry_diff <- diff(not_boundry, lag = 1)
      
      not_boundry_cumsum <- cumsum(not_boundry)
      
      # df <- data.frame(not_boundry = not_boundry, not_boundry_diff = c(not_boundry_diff, NA), not_boundry_cumsum = not_boundry_cumsum)
      
      ### Values used for filling in the boundry NAs - swith from FALSE to TRUE
      replacement_indx <- which(not_boundry_diff == 1) + 1
      
      replaced_indx <- which(!not_boundry)
      
      replaced_freq <- as.numeric(table(not_boundry_cumsum[replaced_indx]))
      
      moderation_boundry  <- moderation[rep(replacement_indx, times = replaced_freq), , drop = FALSE]
      
      moderation[!not_boundry, ] <- moderation_boundry
      
      moderation <- moderation[oo, , drop = FALSE]
      
      loglik <- loglik + priorN * moderation ### like in egdeR estimateTagwiseDisp
      # loglik <- (loglik + priorN * moderation)/(1 + priorN) ### like in edgeR dispCoxReidInterpolateTagwise
      
    }
  )
  
}


out <- edgeR::maximizeInterpolant(splinePts, loglik)


#### set NA for genes that tagwise disp could not be calculated 
dispersion <- rep(NA, length(inds))
names(dispersion) <- names(counts)
dispersion[not_nas] <- disp_init * 2^out








































