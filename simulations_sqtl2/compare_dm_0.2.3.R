#######################################################

# Created 11 Sep 2015
# BioC 3.1

#######################################################


setwd("/home/gosia/multinomial_project/simulations_sqtl2_hsapiens_noDE_noNull")


library(ggplot2)
library(DM)


DM_out <- DM_plots_out <- "dm_0_2_3_plots/"
dir.create(DM_plots_out, showWarnings=F, recursive=T)

results_dm <- "dm_0_2_3/"


##############################################################################################################
# ->>>>>> load results
##############################################################################################################


results <- list()


### sqtlseeker

rt <- read.table(paste0("sqtlseeker_2_0/results_min_dispersion001/results.txt"), header = T, as.is = TRUE)
head(rt)

rt <- rt[, c("geneId", "pv", "qv")]

colnames(rt) <- c("gene_id", "pvalue", "adj_pvalue")

results[[1]] <- rt



### DM

rt <- read.table(paste0("dm_0_2_3/results.txt"), header = TRUE, as.is = TRUE)

head(rt)

rt <- rt[,c("gene_id", "pvalue" ,"adj_pvalue")]

colnames(rt) <- c("gene_id", "pvalue", "adj_pvalue")
head(rt)

results[[2]] <- rt  



split_levels <- data.frame(method = c("sqtlseeker", "DM"))

#######################################################
# load simulation info
#######################################################


# simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, as.is = TRUE)

load("3_truth/simulation_details.Rdata")

status <- unique(simulation_details[, c("gene_id", "gene_ds_status")])

colnames(status) <- c("gene_id", "status")


#######################################################
# plot_TPRFDR function
#######################################################




plot_TPRFDR <- function(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = NULL, xylim_one = FALSE){
  
  split_levels <- split_levels[, c(plot_levels, facet_levels), drop = FALSE]
  
  
  TPRFDRlist <- lapply(1:nrow(split_levels), function(r){
    # r = 1
    
    TPRFDR <- data_TPRFDR[[r]]
    
    TPRFDR <- cbind(TPRFDR, split_levels[rep(r, nrow(TPRFDR)), , drop = FALSE])
    
    })
  

  TPRFDR <- do.call(rbind, TPRFDRlist)
  
  TPRFDR$white <- ifelse(TPRFDR$FDR <= TPRFDR$threshold, NA, TPRFDR$TPR)
  
  pointsize <- 2.5

  ggp <- ggplot(data = TPRFDR, aes_string(x = "FDR", y = "TPR", group = plot_levels, colour = plot_levels)) +
  theme_bw() +
  xlab("FDR") +
  geom_line(size = 1.5, na.rm=TRUE) +
  geom_vline(aes(xintercept = threshold), linetype = "dashed") + 
  geom_point(size = pointsize + 1, shape = 19, na.rm=TRUE) + 
  geom_point(aes_string(y = "white"), size = pointsize, shape = 21, fill = "white", na.rm=TRUE) + 
  theme(axis.text=element_text(size = 16), axis.title = element_text(size = 18, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12), strip.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, shape = NA), nrow = 3)) 
  
  if(xylim_one)
  ggp <- ggp + coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
  
  if(!is.null(plot_colors) && nlevels(TPRFDR[, plot_levels]) == length(plot_colors))   
  ggp <- ggp + scale_color_manual(values = plot_colors)
  else
  ggp <- ggp + scale_color_manual(values = colorb(nlevels(TPRFDR[, plot_levels])))
  
  if(length(facet_levels) == 1)
  ggp <- ggp + facet_wrap(reformulate(facet_levels[1]))
  
  if(length(facet_levels) == 2)
  ggp <- ggp + facet_grid(reformulate(facet_levels[1], facet_levels[2]))

  
  # pdf("./TPRFDR.pdf")
  
  print(ggp)

  # dev.off()
  

}


#######################################################
# generate TPR vs achieved FDR plots
#######################################################



data_TPRFDR <- DM::calculate_TPRFDR(results, status, thresholds = c(0.01, 0.05, 0.1))


levels(split_levels$method)

plot_levels <- "method"
facet_levels <- numeric()

plot_colors <- c("dodgerblue3", "orange")




pdf(paste0(DM_plots_out,"TPRFDR.pdf"), width = 7, height = 7)

plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = FALSE)

dev.off()



pdf(paste0(DM_plots_out,"TPRFDR1.pdf"), width = 7, height = 7)

plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = TRUE)

dev.off()





#######################################################
# generate ROC plots
#######################################################



data_ROCx <- DM::calculate_ROCx(results, status)


levels(split_levels$method)

plot_levels <- "method"
facet_levels <- numeric()

plot_colors <- c("dodgerblue3", "orange")




pdf(paste0(DM_plots_out,"ROC.pdf"), width = 7, height = 7)

plot_ROCx(data_ROCx, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = FALSE)

dev.off()



pdf(paste0(DM_plots_out,"ROC1.pdf"), width = 7, height = 7)

plot_ROCx(data_ROCx, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = TRUE)

dev.off()







##############################################################################################################
# ->>>>>> load results
##############################################################################################################

gene_expr_st1 <- unique(simulation_details[simulation_details$gene_ds_status == 1 , c("gene_id", "expected_gene_count_gr1")])


gene_expr_st1$gene_expr <- Hmisc::cut2(gene_expr_st1[, 2], g = 3)

levels(gene_expr_st1$gene_expr)
table(gene_expr_st1$gene_expr)


gene_expr_st1$gene_expr <- Hmisc::cut2(gene_expr_st1[, 2], cuts = c(1000, 2500))

levels(gene_expr_st1$gene_expr)
table(gene_expr_st1$gene_expr)





gene_expr <- unique(simulation_details[, c("gene_id", "expected_gene_count_gr1")])
rownames(gene_expr) <- gene_expr$gene_id

gene_expr$gene_expr <- Hmisc::cut2(gene_expr[, 2], cuts = c(1000, 2500))

levels(gene_expr$gene_expr)
table(gene_expr$gene_expr)





results <- list()


### sqtlseeker

rt <- read.table(paste0("sqtlseeker_2_0/results_min_dispersion001/results.txt"), header = T, as.is = TRUE)
head(rt)

rt <- rt[, c("geneId", "pv", "qv")]

colnames(rt) <- c("gene_id", "pvalue", "adj_pvalue")

results[[1]] <- split(rt, gene_expr[rt$gene_id, "gene_expr"]) 



### DM

rt <- read.table(paste0("dm_0_2_3/results.txt"), header = TRUE, as.is = TRUE)

head(rt)

rt <- rt[,c("gene_id", "pvalue" ,"adj_pvalue")]

colnames(rt) <- c("gene_id", "pvalue", "adj_pvalue")
head(rt)

results[[2]] <- split(rt, gene_expr[rt$gene_id, "gene_expr"])  



split_levels <- data.frame(method = c("sqtlseeker", "DM"))




  #######################################################
  # load simulation info
  #######################################################


  load("3_truth/simulation_details.Rdata")

  status <- unique(simulation_details[, c("gene_id", "gene_ds_status")])

  colnames(status) <- c("gene_id", "status")


  status <- split(status, gene_expr[status$gene_id, "gene_expr"])


  #######################################################
  # generate TPR vs achieved FDR plots
  #######################################################


  data_TPRFDR_list <- list()

  for(i in 1:3){
    
    data_TPRFDR_list[[i]] <- DM::calculate_TPRFDR(results = lapply(results, function(r) r[[i]]), status[[i]], thresholds = c(0.01, 0.05, 0.1))

    
  }


  data_TPRFDR <- unlist(data_TPRFDR_list, recursive = FALSE)

  nr_method <- nrow(split_levels)

  split_levels <- split_levels[rep(1:nrow(split_levels), 3), , drop = FALSE]
  
  split_levels$gene_expr <- factor(rep(levels(gene_expr$gene_expr), each = nr_method), levels = levels(gene_expr$gene_expr))


  levels(split_levels$method)

  plot_levels <- "method"
  facet_levels <- c("gene_expr")

  plot_colors <- c("dodgerblue3", "orange")



  pdf(paste0(DM_plots_out,"TPRFDR_gene_expr.pdf"), width = 10.5, height = 5)

  plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = FALSE)

  dev.off()



  pdf(paste0(DM_plots_out,"TPRFDR_gene_expr1.pdf"), width = 10.5, height = 5)

  plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = TRUE)

  dev.off()





















