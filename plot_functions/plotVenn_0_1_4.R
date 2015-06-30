#######################################################
# venn diagrams --- my function
#######################################################
library(VennDiagram)


plotVenn <- function(venn.list, venn.colors, venn.subset = names(venne.list), plotPath = "Venn_Diagram.pdf", cex=1, cat.cex=1.4, lwd=2, lty=1, alpha=0.5, margin=0.1){
  
  colors <- venn.colors[venn.subset]
  
  if(length(venn.subset)==4)
    colors <- colors[c(1, 3, 4 ,2)]
  
  venn.d <- venn.diagram(venn.list[venn.subset], filename = NULL, fill = venn.colors[venn.subset], col = colors, cex=cex, cat.cex=cat.cex, lwd=lwd, lty=lty, alpha=alpha, margin = margin)
  
  pdf(plotPath)
  grid.draw(venn.d)
  dev.off()
  
  
}



library(limma)

plotVenn2 <- function(venn.list, venn.colors, venn.subset = names(venne.list), plotPath = "Venn_Diagram.pdf", cex=1, cat.cex=1.4, lwd=2, margin=0.1){
	
	venn.list <- venn.list[venn.subset]
	
	elem <- unique(unlist(venn.list))
	
	res <- lapply(venn.list, function(x){
		
		elem %in% x
		
	})
	
	res <- do.call(cbind, res)
	
	vc <- vennCounts(res)
	
  pdf(plotPath)
	vennDiagram(vc, mar = rep(margin,4), lwd = lwd, circle.col = venn.colors[venn.subset], cex = c(cat.cex, cex, cex))
  dev.off()
	
}












