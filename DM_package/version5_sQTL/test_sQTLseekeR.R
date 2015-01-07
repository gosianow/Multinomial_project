################################################################################

# Created 15 Sep 2014
# Updated 03 Oct 2014
# BioC 14

################################################################################

### install vegan
# install.packages("/home/gosia/R/packages/permute_0.8-3.tar.gz", lib="/home/gosia/R/libraries/3.1.0/")
# install.packages("/home/gosia/R/packages/vegan_2.0-10.tar.gz", lib="/home/gosia/R/libraries/3.1.0/")

### installation in terminal
# R14 CMD INSTALL sQTLseekeR-1-2_0.gz

# R14 CMD INSTALL sQTLseekeR-1-3.gz



library(sQTLseekeR)

data(trans.exp.CEU.3genes)
data(genotype.CEU.3genes)


out <- sQTLseekeR(trans.exp.CEU.3genes, genotype.CEU.3genes, outfile.prefix=NULL, nb.perm.max=10000)

head(out)













