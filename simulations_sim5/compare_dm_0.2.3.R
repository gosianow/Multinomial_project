#######################################################

# Created 11 Sep 2015
# BioC 3.1

#######################################################


setwd("/home/gosia/multinomial_project/simulations_sim5_drosophila_noDE_noNull")

# setwd("/home/gosia/multinomial_project/simulations_sim5_hsapiens_noDE_noNull")

# setwd("/home/gosia/multinomial_project/simulations_sim5_hsapiens_withDE_noNull")


library(ggplot2)
library(DM)


DM_out <- DM_plots_out <- "dm_0_2_3_plots/"
dir.create(DM_plots_out, showWarnings=F, recursive=T)

results_dm <- "dm_0_2_3/"


##############################################################################################################
# ->>>>>> load results
##############################################################################################################


count_methodList <- c("htseq", "kallisto", "htseq_kallisto_filter")
dexseq_resultsList <- c("dexseq_htseq_nomerge", "dexseq_kallisto", "INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast5")

filter_methodList <- c("filtering_dexseq", "filtering_min3prop0_01min6cpm1maxInf") 


split_levels <- list()

results <- list()

counter <- 1

for(i in 1:2){
	# i = 3
  
	for(j in 1:length(filter_methodList)){
		# j = 1
    
		count_method <- count_methodList[i]
		filter_method <- filter_methodList[j]

		####################### results produced by Charlotte


		rt <- read.table(paste0("4_results/", dexseq_resultsList[i], ".txt"), header = T, as.is = TRUE)
		head(rt)

		colnames(rt) <- c("gene_id", "adj_pvalue")

		results[[counter]] <- rt
		split_levels[[counter]] <- data.frame(method = "dexseq", counting = count_method, filtering = filter_method)
    
		counter <- counter + 1


		####################### DM results  

		res_path <- paste0(results_dm, count_method, "/", filter_method, "/")

		files <- list.files(path = res_path, pattern = "_results.txt", full.names = TRUE)
    
		print(files)
    
		res  <- gsub(pattern = "_results.txt", replacement = "", x = basename(files))
		# res  <- gsub(pattern = paste0(count_method, "_"), replacement = "", res)
		# res  <- gsub(pattern = "constrOptim2G_", replacement = "", res)
		# res 

		for(k in 1:length(files)){
			# k = 1
      
			# rt <- read.table(paste0(res_path, files[k]), header = TRUE, as.is = TRUE)
			rt <- read.table(files[k], header = TRUE, as.is = TRUE)
      
			head(rt)
      
			rt <- rt[,c("gene_id", "pvalue" ,"adj_pvalue")]
      
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
# generate TPR vs achieved FDR plots
#######################################################



data_TPRFDR <- DM::calculate_TPRFDR(results, status, thresholds = c(0.01, 0.05, 0.1))


levels(split_levels$method)

plot_levels <- "method"
facet_levels <- c("counting", "filtering")


plot_colors <- c("orange", "chartreuse4", "lightblue", "palevioletred2", "firebrick1", "dodgerblue3", "firebrick4")[1:nlevels(split_levels$method)]




pdf(paste0(DM_plots_out,"TPRFDR.pdf"), width = 7, height = 7)

plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = FALSE)

dev.off()



pdf(paste0(DM_plots_out,"TPRFDR1.pdf"), width = 7, height = 7)

plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = TRUE)

dev.off()





#######################################################
# generate venn diagrams 
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
  guides(colour = guide_legend(override.aes = list(size = 1.5, shape = NA), nrow = 1)) 
  
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





##############################################################################################################
# ->>>>>> load results
##############################################################################################################


count_methodList <- c("htseq", "kallisto", "htseq_kallisto_filter")
dexseq_resultsList <- c("dexseq_htseq_nomerge", "dexseq_kallisto", "INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast5")

### for Human
cutsList <- list(c(10, 25), c(4, 7), c(10, 20))

### for Drosophila
# cutsList <- list(c(5, 15), c(4, 7), c(5, 15))


filter_methodList <- c("filtering_dexseq") 



for(i in 1:2){
  # i = 1
  
  
split_levels <- list()

results <- list()

counter <- 1


  cuts <- cutsList[[i]]


  j = 1

  count_method <- count_methodList[i]
  filter_method <- filter_methodList[j]


  load(paste0(results_dm, count_method, "/", filter_method, "/d.Rdata"))

  nr_iso <- Hmisc::cut2(width(d@counts), cuts = cuts)
  names(nr_iso) <- names(d@counts)

  levels(nr_iso)

  table(nr_iso)


  # nr_iso <- Hmisc::cut2(width(d@counts), g = 3)
  # names(nr_iso) <- names(d@counts)

  # levels(nr_iso)

  # table(nr_iso)



  ####################### results produced by Charlotte


  rt <- read.table(paste0("4_results/", dexseq_resultsList[i], ".txt"), header = T, as.is = TRUE)
  head(rt)

  colnames(rt) <- c("gene_id", "adj_pvalue")

  rt <- split(rt, nr_iso)

  results[[counter]] <- rt

  split_levels[[counter]] <- data.frame(method = rep("dexseq", 1), counting = rep(count_method, 1), filtering = rep(filter_method, 1))

  counter <- counter + 1


  ####################### DM results  

  res_path <- paste0(results_dm, count_method, "/", filter_method, "/")

  files <- list.files(path = res_path, pattern = "_results.txt", full.names = TRUE)

  print(files)

  res  <- gsub(pattern = "_results.txt", replacement = "", x = basename(files))

  for(k in 1:length(files)){
    # k = 1
    
    # rt <- read.table(paste0(res_path, files[k]), header = TRUE, as.is = TRUE)
    rt <- read.table(files[k], header = TRUE, as.is = TRUE)
    
    head(rt)
    
    rt <- rt[,c("gene_id", "pvalue" ,"adj_pvalue")]
    
    colnames(rt) <- c("gene_id", "pvalue", "adj_pvalue")
    head(rt)
    
    rt <- split(rt, nr_iso)
    
    results[[counter]] <- rt  
    split_levels[[counter]] <- data.frame(method = rep(res[k], 1), counting = rep(count_method, 1), filtering = rep(filter_method, 1))
    
    counter <- counter + 1
    
  }



  split_levels <- do.call(rbind, split_levels)

  #######################################################
  # load simulation info
  #######################################################


  simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, as.is = TRUE)

  status <- unique(simulation_details[, c("gene_id", "gene_ds_status")])

  colnames(status) <- c("gene_id", "status")

  status <- split(status, nr_iso[status$gene_id])


  #######################################################
  # generate TPR vs achieved FDR plots
  #######################################################

  data_TPRFDR_list <- list()

  for(i in 1:3){
    
    data_TPRFDR_list[[i]] <- DM::calculate_TPRFDR(results = lapply(results, function(r) r[[i]]), status[[i]], thresholds = c(0.01, 0.05, 0.1))

    
  }


  data_TPRFDR <- unlist(data_TPRFDR_list, recursive = FALSE)

  nr_method <- nrow(split_levels)

  split_levels <- split_levels[rep(1:nrow(split_levels), 3), ]
  split_levels$nr_iso <- factor(rep(levels(nr_iso), each = nr_method), levels = levels(nr_iso))


  levels(split_levels$method)

  plot_levels <- "method"
  facet_levels <- c("nr_iso")


  plot_colors <- c("orange", "chartreuse4", "lightblue", "palevioletred2", "firebrick1", "dodgerblue3", "firebrick4")[1:nlevels(split_levels$method)]




  pdf(paste0(DM_plots_out,"TPRFDR_nr_iso_", count_method,".pdf"), width = 10.5, height = 5)

  plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = FALSE)

  dev.off()



  pdf(paste0(DM_plots_out,"TPRFDR_nr_iso1_",count_method,".pdf"), width = 10.5, height = 5)

  plot_TPRFDR(data_TPRFDR, split_levels, plot_levels, facet_levels, plot_colors = plot_colors, xylim_one = TRUE)

  dev.off()




  
}




















