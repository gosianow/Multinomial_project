### R CMD build & R CMD check

setwd("/home/gosia/R/R_Multinomial_project/DM_package_devel/")

R <- shQuote(file.path(R.home(component="bin"), "R"))

Sys.setenv("R_TESTS"="") # needed for R CMD check; thanks for the tip, Hadley




## R30 CMD build DM
system(paste(R, "CMD build DM"))


## R30 CMD check DM DM_1.0.tar.gz
# system(paste(R, "CMD check DM_1.0.tar.gz"))


## R30 CMD INSTALL DM DM_1.0.tar.gz
system(paste(R, "CMD INSTALL DM_0.1.5.tar.gz"))





