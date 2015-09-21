#######################################################
# BioC 3.1
# Created 11 Sep 2015
# Updated 21 Sep 2015
# - add plots of p-values DM vs sqtlseeker

#######################################################


setwd("/home/gosia/multinomial_project/simulations_sqtl3_hsapiens_noDE_noNull")


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





results_all <- merge(results[[1]], results[[2]], by = "gene_id", suffixes = c("_sqtlseeker", "_DM"))


ggp <- ggplot(data = results_all, aes(x = pvalue_sqtlseeker, y = pvalue_DM)) +
geom_point()


ggp2 <- ggplot(data = results_all, aes(x = pvalue_sqtlseeker, y = pvalue_DM)) +
geom_point() + scale_x_log10() + scale_y_log10()


pdf(paste0(DM_plots_out,"scatter_pvalues.pdf"), width = 7, height = 7)

print(ggp)
print(ggp2)

dev.off()





#######################################################
# load simulation info
#######################################################


# simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, as.is = TRUE)

load("3_truth/simulation_details.Rdata")

status <- unique(simulation_details[, c("gene_id", "gene_ds_status")])

colnames(status) <- c("gene_id", "status")



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




#######################################################
# generate venne diagrams
#######################################################

library(VennDiagram)

calculate_venn <- function(results, status, threshold = 0.05){
  
  
  venn <- lapply(results, function(r){
    
    r[which(r$adj_pvalue < threshold), "gene_id"]
    
    })
  
  venn[["status"]] <- status[which(status$status == 1), "gene_id"]
  
  return(venn)
}



plot_venn <- function(data_venn, plot_results, split_levels, plot_levels, plot_colors = NULL, plot_status = TRUE){
  
  
  if(plot_status && "status" %in% names(data_venn)){
    stopifnot(length(plot_results) <= 4)
    
    data_venn_tmp <- data_venn[plot_results]
    names(data_venn_tmp) <- split_levels[plot_results, plot_levels]
    data_venn_tmp[["status"]] <- data_venn[["status"]]
    
    if(!is.null(plot_colors))
    colors <- c(plot_colors, "grey")
    else
    colors <- c(colorb(length(plot_results)), "grey")
    
    }else{
      stopifnot(length(plot_results) <= 5)
      
      data_venn_tmp <- data_venn[plot_results]
      names(data_venn_tmp) <- split_levels[plot_results, plot_levels]
      
      if(!is.null(plot_colors))
      colors <- plot_colors
      else
      colors <- colorb(length(plot_results))
    }

    
    
    cex = 1.5
    cat.cex=1.5
    lwd=2
    lty=1
    alpha=0.5
    margin=0.1
    
    if(length(colors) == 4)
    colors <- colors[c(1, 3, 4 ,2)]
    
    venn <- VennDiagram::venn.diagram(data_venn_tmp, filename = NULL, fill = colors, col = colors, cex=cex, cat.cex=cat.cex, lwd=lwd, lty=lty, alpha=alpha, margin = margin, scaled = FALSE)
    
    
    return(grid.draw(venn))
    
    
  }





  data_venn <- calculate_venn(results, status, threshold = 0.05)


  plot_results <- 1:2
  plot_levels <- "method"
  plot_colors <- c("dodgerblue3", "orange")


  pdf(paste0(DM_plots_out,"venn.pdf"), width = 7, height = 7)

  plot_venn(data_venn, plot_results, split_levels, plot_levels, plot_colors = plot_colors, plot_status = TRUE)

  dev.off()



##############################################################################################################
# stratified plots by gene expression 
##############################################################################################################


### decide the cuts: equal ammount of ds genes in each strata

gene_expr_st1 <- unique(simulation_details[simulation_details$gene_ds_status == 1 , c("gene_id", "expected_gene_count_gr1")])


gene_expr_st1$gene_expr <- Hmisc::cut2(gene_expr_st1[, "expected_gene_count_gr1"], g = 3)

levels(gene_expr_st1$gene_expr)
table(gene_expr_st1$gene_expr)


gene_expr_st1$gene_expr <- Hmisc::cut2(gene_expr_st1[, "expected_gene_count_gr1"], cuts = c(1000, 2500))

levels(gene_expr_st1$gene_expr)
table(gene_expr_st1$gene_expr)



### stratification


gene_expr <- unique(simulation_details[, c("gene_id", "gene_ds_status", "expected_gene_count_gr1")])
rownames(gene_expr) <- gene_expr$gene_id

gene_expr$gene_expr <- Hmisc::cut2(gene_expr[, "expected_gene_count_gr1"], cuts = c(1000, 2500))

levels(gene_expr$gene_expr)
table(gene_expr$gene_expr)
table(gene_expr$gene_expr[gene_expr$gene_ds_status == 1])


# names with counts of n and n.ds

new_levels <- paste0(gsub(" ", "", levels(gene_expr$gene_expr)), " (n = ", table(gene_expr$gene_expr), ", n.ds = ", table(gene_expr$gene_expr[gene_expr$gene_ds_status == 1]), ")")

levels(gene_expr$gene_expr) <- new_levels


### results

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





##############################################################################################################
# stratified plots by diff_IsoPct
##############################################################################################################

strat_var <- "diff_IsoPct"

### decide the cuts: equal ammount of ds genes in each strata

strata <- unique(simulation_details[simulation_details$gene_ds_status == 1 , c("gene_id", strat_var)])


strata$strata <- Hmisc::cut2(strata[, strat_var], g = 3)

levels(strata$strata)
table(strata$strata)


strata$strata <- Hmisc::cut2(strata[, strat_var], cuts = c(1/3, 2/3, 1))

levels(strata$strata)
table(strata$strata)



### stratification


strata <- unique(simulation_details[, c("gene_id", "gene_ds_status", strat_var)])
rownames(strata) <- strata$gene_id

strata$strata <- Hmisc::cut2(strata[, strat_var], cuts = c(1/3, 2/3, 1))

levels(strata$strata)
table(strata$strata)
table(strata$strata[strata$gene_ds_status == 1])


# names with counts of n and n.ds

new_levels <- paste0(gsub(" ", "", levels(strata$strata)), " (n = ", table(strata$strata), ", n.ds = ", table(strata$strata[strata$gene_ds_status == 1]), ")")

levels(strata$strata) <- new_levels

levels(strata$strata)


### results

results <- list()


### sqtlseeker

rt <- read.table(paste0("sqtlseeker_2_0/results_min_dispersion001/results.txt"), header = T, as.is = TRUE)
head(rt)

rt <- rt[, c("geneId", "pv", "qv")]

colnames(rt) <- c("gene_id", "pvalue", "adj_pvalue")

results[[1]] <- split(rt, strata[rt$gene_id, "strata"]) 



### DM

rt <- read.table(paste0("dm_0_2_3/results.txt"), header = TRUE, as.is = TRUE)

head(rt)

rt <- rt[,c("gene_id", "pvalue" ,"adj_pvalue")]

colnames(rt) <- c("gene_id", "pvalue", "adj_pvalue")
head(rt)

results[[2]] <- split(rt, strata[rt$gene_id, "strata"])  



split_levels <- data.frame(method = c("sqtlseeker", "DM"))




#######################################################
# load simulation info
#######################################################


load("3_truth/simulation_details.Rdata")

status <- unique(simulation_details[, c("gene_id", "gene_ds_status")])

colnames(status) <- c("gene_id", "status")


status <- split(status, strata[status$gene_id, "strata"])


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

split_levels$strata <- factor(rep(levels(strata$strata), each = nr_method), levels = levels(strata$strata))


levels(split_levels$method)

plot_levels <- "method"
facet_levels <- c("strata")

plot_colors <- c("dodgerblue3", "orange")



pdf(paste0(DM_plots_out,"TPRFDR_",strat_var,".pdf"), width = 10.5, height = 5)

plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = FALSE)

dev.off()



pdf(paste0(DM_plots_out,"TPRFDR_", strat_var,"1.pdf"), width = 10.5, height = 5)

plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = TRUE)

dev.off()
















