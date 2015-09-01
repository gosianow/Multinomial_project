#######################################################

# Created 28 Aug 2015
# BioC 3.1

# Comprare DM_0.1.5 (different dispersion estimators) with DEXSeq
# Use htseq and kallisto counts
# Plots per count method and per filtering

# Create new ggplot versions of ROC and FDRTPR plots

#######################################################


setwd("/home/gosia/multinomial_project/simulations_sim5_drosophila_noDE_noNull")


library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)


DM_out <- DM_plots_out <- "dm_0_1_5_new_ggplot_ROC_FDR/"
dir.create(DM_plots_out, showWarnings=F, recursive=T)

results_dm <- "dm_0_1_5/"


##############################################################################################################
# ->>>>>> load results
##############################################################################################################


count_methodList <- c("htseq", "kallisto")

filter_methodList <- c("Filtering_DEXSeq", "Filtering_min3prop0_01min6cpm1maxNA") 

split_levels <- list()

results <- list()

counter <- 1

for(i in 1:length(count_methodList)){
	# i = 1
  
	for(j in 1:length(filter_methodList)){
		# j = 1
    
		count_method <- count_methodList[i]
		filter_method <- filter_methodList[j]

		####################### results produced by Charlotte


		rt <- read.table(paste0("4_results/dexseq_", count_method, ".txt"), header = T, as.is = TRUE)
		head(rt)

		colnames(rt) <- c("gene_id", "adj_pvalue")

		results[[counter]] <- rt
		split_levels[[counter]] <- data.frame(method = "dexseq", counting = count_method, filtering = filter_method)
    
		counter <- counter + 1


		####################### DM_0.1.5 results on htseq counts and bitseq 

		res_path <- paste0(results_dm, count_method, "/", filter_method, "/")

		files <- list.files(path = res_path, pattern = "_results.xls", full.names = TRUE)
    
		print(files)
    
		res  <- gsub(pattern = "_results.xls", replacement = "", x = basename(files))
		res  <- gsub(pattern = paste0(count_method, "_"), replacement = "", res)
		res  <- gsub(pattern = "constrOptim2G_", replacement = "", res)
    
		res 

		for(k in 1:length(files)){
			# k = 1
      
			# rt <- read.table(paste0(res_path, files[k]), header = TRUE, as.is = TRUE)
			rt <- read.table(files[k], header = TRUE, as.is = TRUE)
      
			head(rt)
      
			rt <- rt[,c("geneID", "pValue" ,"FDR")]
      
			colnames(rt) <- c("gene_id", "pvalue", "adj_pvalue")
			head(rt)
      
			results[[counter]] <- rt  
			split_levels[[counter]] <- data.frame(method = res[k], counting = count_method, filtering = filter_method)
      
			counter <- counter + 1
      
      
		}
    
	}
  
}


split_levels <- do.call(rbind, split_levels)

#######################################################
# load simulation info
#######################################################


simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, as.is = TRUE)

status <- unique(simulation_details[, c("gene_id", "gene_ds_status")])

colnames(status) <- c("gene_id", "status")





#######################################################
# generate ROCx plots 
#######################################################



calculate_ROCx <- function(results, status){
	# results: list of data.frames with "gene_id", "pvalue" and "adj_pvalue" columns
	# status: data.frame with "gene_id" and "status" columns. "status" must consists of 0 and 1 (or FALSE, TRUE)
   
	status <- status[complete.cases(status[, c("gene_id", "status")]), , drop = FALSE]
  
	P <- sum(status$status == 1)
	N <- sum(status$status == 0)
  
	Xlist <- list()
	ROClist <- list()
  
	for(m in 1:length(results)){
		# m = 2
		print(m)
    
		if(! "pvalue" %in% colnames(results[[m]])){
      
			message("No 'pvalue' column in results ", m)
			ROClist[[m]] <- data.frame(FPR = NA, TPR = NA)
			Xlist[[m]] <- data.frame(FPR = NA, TPR = NA) 
      
		}else{
      
			sc <- merge(status, results[[m]], by = "gene_id", all.x = TRUE)
      
			NAs <- is.na(sc[, "pvalue"]) | is.na(sc[, "adj_pvalue"])
      
			pvs <- sc$pvalue[!NAs]
			apvs <- sc$adj_pvalue[!NAs]
			sts <- sc$status[!NAs]

			ord <- order(pvs, decreasing = FALSE)
			sts <- sts[ord]
			pvs <- pvs[ord]
      
			TPRv <- cumsum(sts) / P
			FPRv <- cumsum(!sts) / N
			ROClist[[m]] <- data.frame(FPR = FPRv, TPR = TPRv)
      
			### what is the TPR if the error was controlled 
			TPR <- sum(apvs[sts == 1] < 0.05) / P    
			Xlist[[m]] <- data.frame(FPR = approx(TPRv, FPRv, xout=TPR)$y, TPR = TPR)    
      
      
		}
    

	}
  
	return(list(ROClist = ROClist, Xlist = Xlist))
   
}



data_ROCx <- calculate_ROCx(results, status)



colorb <- function(n){
  
	clrs <- c("dodgerblue3", "maroon2",  "forestgreen",  "blueviolet", "firebrick3", "deepskyblue",  "orchid2", "chartreuse3", "tomato" , "slateblue1")
  
	nc <- length(clrs)
  
	if(n > nc)
		clrs <- rep(clrs, ceiling(n/nc))
    
	clrs[1:n]
  
	# colorRampPalette(clrs)(n)

}




plot_levels <- "method"
facet_levels <- c("counting", "filtering")
plot_colors <- colorb(4)



plot_ROCx <- function(data_ROCx, split_levels, plot_levels, facet_levels, plot_colors = NULL, xylim_one = TRUE){
  
	split_levels <- split_levels[, c(plot_levels, facet_levels), drop = FALSE]
  
	ROClist <- lapply(1:nrow(split_levels), function(r){
		# r = 1
    
		ROC <- data_ROCx$ROClist[[r]]
    
		ROC <- cbind(ROC, split_levels[rep(r, nrow(ROC)), , drop = FALSE])
    
	})
  
	Xlist <- lapply(1:nrow(split_levels), function(r){
		# r = 1
    
		X <- data_ROCx$Xlist[[r]]
    
		X <- cbind(X, split_levels[rep(r, nrow(X)), , drop = FALSE])
    
	})
  
	X <- do.call(rbind, Xlist)
	ROC <- do.call(rbind, ROClist)
    
    
	ggp <- ggplot(data = ROC, aes_string(x = "FPR", y = "TPR", group = plot_levels, colour = plot_levels)) +
	theme_bw() +
	geom_line(size = 1.5, na.rm=TRUE) +
	theme(axis.text=element_text(size = 16), axis.title = element_text(size = 18, face = "bold"), legend.justification = c(1, 0), legend.position = "bottom", legend.title = element_blank()) +
	guides(colour = guide_legend(override.aes = list(size = 1.5, shape = NA))) +
	geom_point(data = X, aes_string(x = "FPR", y = "TPR", group = plot_levels, colour = plot_levels), size = 8, shape = "X", na.rm=TRUE) 
  
  if(xylim_one)
  ggp <- ggp + coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
      
	if(!is.null(plot_colors) && nlevels(ROC[, plot_levels]) == length(plot_colors))     
		ggp <- ggp + scale_color_manual(values = plot_colors)
  else
    ggp <- ggp + scale_color_manual(values = colorb(nlevels(TPRFDR[, plot_levels])))
      
	if(length(facet_levels) == 1)
		ggp <- ggp + facet_wrap(reformulate(facet_levels[1]))
      
	if(length(facet_levels) == 2)
		ggp <- ggp + facet_grid(reformulate(facet_levels[1], facet_levels[2]))

    
	# pdf("./ROC.pdf")
   
	print(ggp)

	# dev.off()
  
  
}



pdf("./ROC.pdf")
 
plot_ROCx(data_ROCx, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = FALSE)

dev.off()





#######################################################
# generate TPR vs achieved FDR plots
#######################################################



calculate_TPRFDR <- function(results, status, thresholds = c(0.01, 0.05, 0.1)){
	# results: list of data.frames with "gene_id" and "adj_pvalue" columns
	# status: data.frame with "gene_id" and "status" columns. "status" must consists of 0 and 1 (or FALSE, TRUE)
   
	status <- status[complete.cases(status[, c("gene_id", "status")]), , drop = FALSE]
  
	TPRFDRlist <- list()

	for(m in 1:length(results)){
		# m = 1
		# print(m)
    
		if(! "adj_pvalue" %in% colnames(results[[m]])){
      
			message("No 'adj_pvalue' column in results ", m)
			TPRFDRlist[[m]] <- data.frame(threshold = NA, FDR = NA, TPR = NA)

      
		}else{
      
			sc <- merge(status, results[[m]], by = "gene_id", all.x = TRUE)
      
			# mm <- match(status$gene_id, results[[m]]$gene_id)
       
			apvs <- sc$adj_pvalue
			apvs[is.na(apvs)] <- 1
			sts <- sc$status
      
			q <- length(thresholds)
			TPR <- rep(0, q)
			FDR <- rep(0, q)
      
			for(i in 1:q){
				# i=1
				sts_est <- as.numeric(apvs < thresholds[i])
        
				TP <- sum(sts==1 & sts_est==1)
				FP <- sum(sts==0 & sts_est==1)
				FN <- sum(sts==1 & sts_est==0)
        
				TPR[i] <- TP/(TP+FN)
				FDR[i] <- FP/(FP+TP)
        
			}  
      
			TPRFDRlist[[m]] <- data.frame(threshold = thresholds, FDR = FDR, TPR = TPR)
      
		}
    

	}
  
	return(TPRFDRlist)
  
  
}



data_TPRFDR <- calculate_TPRFDR(results, status, thresholds = c(0.01, 0.05, 0.1))



plot_TPRFDR <- function(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = NULL, xylim_one = TRUE){
  
	split_levels <- split_levels[, c(plot_levels, facet_levels), drop = FALSE]
  
  
	TPRFDRlist <- lapply(1:nrow(split_levels), function(r){
		# r = 1
    
		TPRFDR <- data_TPRFDR[[r]]
    
		TPRFDR <- cbind(TPRFDR, split_levels[rep(r, nrow(TPRFDR)), , drop = FALSE])
    
	})
  

	TPRFDR <- do.call(rbind, TPRFDRlist)
  
  TPRFDR$white <- ifelse(TPRFDR$FDR <= TPRFDR$threshold, NA, TPRFDR$TPR)
  
  pointsize <- 4

	ggp <- ggplot(data = TPRFDR, aes_string(x = "FDR", y = "TPR", group = plot_levels, colour = plot_levels)) +
	theme_bw() +
	xlab("Achieved FDR") +
	geom_line(size = 1.5, na.rm=TRUE) +
	geom_vline(aes(xintercept = threshold), linetype = "dashed") + 
	geom_point(size = pointsize + 1, shape = 19, na.rm=TRUE) + 
	geom_point(aes_string(y = "white"), size = pointsize, shape = 21, fill = "white", na.rm=TRUE) + 
	theme(axis.text=element_text(size = 16), axis.title = element_text(size = 18, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12), strip.text = element_text(size = 12)) +
	guides(colour = guide_legend(override.aes = list(size = 1.5, shape = NA), nrow = 2)) 
  
  if(xylim_one)
  ggp <- ggp + coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1))
      
	if(!is.null(plot_colors) && nlevels(TPRFDR[, plot_levels]) == length(plot_colors)){
		ggp <- ggp + scale_color_manual(values = plot_colors)
  }else{
    ggp <- ggp + scale_color_manual(values = colorb(nlevels(TPRFDR[, plot_levels])))
  }
      
	if(length(facet_levels) == 1)
		ggp <- ggp + facet_wrap(reformulate(facet_levels[1]))
      
	if(length(facet_levels) == 2)
		ggp <- ggp + facet_grid(reformulate(facet_levels[1], facet_levels[2]))

    
	pdf("./TPRFDR.pdf")
   
	print(ggp)

	dev.off()
  

}







pdf("./TPRFDR.pdf")
  
plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = FALSE)

dev.off()








#######################################################
# generate venn diagrams 
#######################################################



























