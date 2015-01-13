#######################################################
# venn diagrams --- my function
#######################################################
library(VennDiagram)


plotVenn <- function(venne.genes, colors, venn.methods, name1="", name2="", cex=1, cat.cex=1.4, lwd=2, lty=1, alpha=0.5, margin=0.05, out.dir){
  
  colors.col <- colors[venn.methods]
  
  if(length(venn.methods)==4)
    colors.col <- colors.col[c(1, 3, 4 ,2)]
  
  venn.d <- venn.diagram(venne.genes[venn.methods], filename=NULL, fill = colors[venn.methods], col=colors.col, cex=cex, cat.cex=cat.cex, lwd=lwd, lty=lty, alpha=alpha, margin=margin)
  
  pdf(paste0(out.dir, name1, "Venn_Diagr_", name2,".pdf"))
  grid.draw(venn.d)
  dev.off()
  
  
}
