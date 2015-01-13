# BioC 214

# Created 8 Jan 2014
# Modyfied 8 Jan 2014


##############################################################################################################

# compare DM SQTL with sQTLseekeR results

##############################################################################################################

setwd("/home/Shared/data/seq/GEUVADIS/")


resSeeker <- read.table("sQTLs-GEUVADIS-FDR10-annotation/sQTLs-FDR10-CEU.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

resSeekerSign <- resSeeker[resSeeker$FDR < 0.05, ] ## 2869


resDM <- read.table("DM_sQTLanalysis/dgeSQTL_results.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

resDMSign <- resDM[resDM$FDR < 0.05 & !is.na(resDM$SNP_id), ] ## 205198


#######################################################
# generate venn diagrams 
#######################################################

source("/home/gosia/R/R_Multinomial_project/Plot_Functions/plotVenn.R")

## TP as a separate circle
venne.list <- list()
venne.list[["Seeker"]] <- resSeekerSign[,"snpId"]
# venne.list[["Seeker"]] <- intersect(resSeekerSign[,"snpId"], resDM[, "SNP_id"])
venne.list[["DM"]] <- na.omit(resDMSign[,"SNP_id"])




  plotVenn(venne.list, colors=c(Seeker="blue", DM="orange"), venn.methods = names(venne.list), margin=0.1, cat.cex=0.8, cex=1.7, out.dir="", name2="")
  












































