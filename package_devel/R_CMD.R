################################################################################
### R CMD build & R CMD check from R
################################################################################

setwd("/home/gosia/R/multinomial_project/package_devel/")

R <- shQuote(file.path(R.home(component="bin"), "R"))

Sys.setenv("R_TESTS"="") # needed for R CMD check; thanks for the tip, Hadley


## R31 CMD build DM
system(paste(R, "CMD build DM"))


## R31 CMD check DM 
# system(paste(R, "CMD check DM_1.0.tar.gz"))


## R31 CMD INSTALL DM 
system(paste(R, "CMD INSTALL DM_0.1.5.tar.gz"))



################################################################################
### R CMD build & R CMD check from terminal 
################################################################################

# cd R/multinomial_project/package_devel

# R CMD build --no-build-vignettes DRIMSeq

# R CMD INSTALL DRIMSeq_0.3.1.tar.gz


### Inportant: use tar.gz thwn it does not contail files from .buildignore

# R CMD check --no-build-vignettes DRIMSeq_0.3.1.tar.gz

# R CMD BiocCheck DRIMSeq_0.3.1.tar.gz




# R CMD build . && R CMD check *tar.gz && R CMD BiocCheck *tar.gz

################################################################################
### compile vignette
################################################################################
# cd DRIMSeq/vignettes

# R CMD Sweave --engine=knitr::knitr --pdf DRIMSeq.Rnw





