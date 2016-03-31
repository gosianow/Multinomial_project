library(VennDiagram)

colorb <- function(n){
  
  clrs <- c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1" ,"blueviolet", "firebrick2", "deepskyblue",  "orchid2", "chartreuse3", "gold", "slateblue1", "tomato" , "blue", "magenta", "green3", "yellow", "purple3", "red" ,"darkslategray1", "lightpink1", "lightgreen", "khaki1", "plum3", "salmon")
  
  nc <- length(clrs)
  
  if(n > nc)
    clrs <- rep(clrs, ceiling(n/nc))
  
  clrs[1:n]
  
  # colorRampPalette(clrs)(n)
  
}


# nb <- 24
# barplot(rep(1, nb), col = colorb(nb))
# dev.off()


#' Venn diagram
#' 
#' Venn diagram overlaps for a list of result data frames.
#'  
#' @inheritParams plotROCx
#' @param threshold Threshold for the FDR.
#' @return 
#' \code{calculateVenn} returns a list of vectors with names of the genes that are significant \code{adj_pvalue < threshold} for different methods. If status is not \code{NULL}, the last element of the output is a list of genes with status equal 1.
#' 
#' \code{plotVenn} returns a venn diagram overlaps of significant genes for different methods. \code{plotVenn} is a wrapper function that uses \code{\link[VennDiagram]{venn.diagram}} for plotting.
#' 
#'  @author Malgorzata Nowicka
#'  @name plotVenn
#' @export
calculateVenn <- function(results, status = NULL, threshold = 0.05){
  
  venn <- lapply(results, function(r){
    
    if(!all(c("gene_id", "adj_pvalue") %in% colnames(r))){
      
      message("Some columns 'gene_id', 'adj_pvalue' are missing in one of the results")
      character()
      
    }else{
      
      as.character(r[which(r$adj_pvalue < threshold), "gene_id"])
      
    }

    })
  
  names(venn) <- NULL
  
  if(!is.null(status)){
    stopifnot(all(c("gene_id", "status") %in% colnames(status)))
    venn[["status"]] <- as.character(status[which(status$status == 1), "gene_id"])
  }

  
  return(venn)
}



################################################################################

#' @param data_venn List structured like the output of \code{calculateVenn}.
#' @param label_var Name of one of the variables in \code{metadata}. Levels of this variable will be used as category names in venn diagram.
#' @param plot_results Numeric vector of maximum length 4 (if \code{plot_status = TRUE}) or 5 (if \code{plot_status = FALSE}) indicating which results should be plotted. 
#' @param plot_status Logical. Whether to plot a circle that corresponds to truly significant genes.
#' @examples 
#' 
#' status <- dataDS_status
#' d <- dataDS_dmDStest
#' 
#' results <- list()
#' results[[1]] <- results(d)
#' metadata <- data.frame(method = "DM")
#' 
#' data_venn <- calculateVenn(results, status = status, threshold = 0.05)
#' 
#' plotVenn(data_venn, plot_results = 1, metadata, label_var = "method", 
#'    label_colors = NULL, plot_status = TRUE)
#' 
#' @seealso \code{\link{dataDS_status}}, \code{\link{dataDS_dmDStest}}, \code{\link{plotROCx}}, \code{\link{plotTPRFDR}}
#' @rdname plotVenn
#' @export
plotVenn <- function(data_venn, plot_results, metadata, label_var, label_colors = NULL, plot_status = TRUE, cex = 1.5, cat.cex=1.5, lwd=2, lty=1, alpha=0.5, margin=0.1){
  
  # cex = 1.5
  # cat.cex=1.5
  # lwd=2
  # lty=1
  # alpha=0.5
  # margin=0.1
  
  stopifnot(is.logical(plot_status))
  if("status" %in% names(data_venn))
    stopifnot(length(data_venn) == nrow(metadata) + 1)
  else
    stopifnot(length(data_venn) == nrow(metadata))
  stopifnot(label_var %in% colnames(metadata))
  
  if(!is.null(label_colors))
    stopifnot(nlevels(metadata[, label_var]) == length(label_colors))
  
  if(plot_status && "status" %in% names(data_venn)){
    stopifnot(length(plot_results) <= 4)
    
    data_venn_tmp <- data_venn[plot_results]
    names(data_venn_tmp) <- metadata[plot_results, label_var]
    data_venn_tmp[["status"]] <- data_venn[["status"]]
    
    if(!is.null(label_colors)){
      colors <- metadata[, label_var]
      levels(colors) <- label_colors
      colors <- c(colors[plot_results], "grey")
    }else{
      colors <- c(colorb(length(plot_results)), "grey")
    }

  }else{
    stopifnot(length(plot_results) <= 5)
    
    data_venn_tmp <- data_venn[plot_results]
    names(data_venn_tmp) <- metadata[plot_results, label_var]
    
    if(!is.null(label_colors)){
      colors <- metadata[, label_var]
      levels(colors) <- label_colors
      colors <- c(colors[plot_results])
    }else{
      colors <- c(colorb(length(plot_results)))
    }
    
  }

  colors_col <- colors
  if(length(colors) == 4)
    colors_col <- colors[c(1, 3, 4, 2)]
    
  cat_pos = rep(0, length(colors))
  if(length(colors) == 3)
  cat_pos = c(0, 0, 180)
  if(length(colors) == 5)
  cat_pos = c(0, 0, 180, 180, 0)
  
  venn <- VennDiagram::venn.diagram(data_venn_tmp, filename = NULL, fill = colors, col = colors_col, cex=cex, cat.cex = cat.cex, lwd=lwd, lty=lty, alpha=alpha, margin = margin, scaled = FALSE, cat.pos = cat_pos)
  
  
  return(grid::grid.draw(venn))
  
  
}

