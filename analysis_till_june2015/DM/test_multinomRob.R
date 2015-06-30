#######################################################
# run example
#######################################################


#install.packages("multinomRob", dependencies=TRUE, lib="/home/gosia/R/libraries/3.0.2")

library(multinomRob)


# make some multinomial data
x1 <- rnorm(50);
x2 <- rnorm(50);
p1 <- exp(2*x1)/(1+exp(2*x1)+exp(x2));
p2 <- exp(x2)/(1+exp(2*x1)+exp(x2));
p3 <- 1 - (p1 + p2);

cbind(p1, p2, p3)

y <- matrix(0, 50, 3);
for (i in 1:50) {
  y[i,] <- rmultinomial(1000, c(p1[i], p2[i], p3[i]));
}

# perturb the first 5 observations
y[1:5,c(1,2,3)] <- y[1:5,c(3,1,2)];
y1 <- y[,1];
y2 <- y[,2];
y3 <- y[,3];

# put data into a dataframe
dtf <- data.frame(x1, x2, y1, y2, y3);

head(dtf)

## Set parameters for Genoud
zz.genoud.parms <- list( pop.size             = 1000,
                         wait.generations      = 10,
                         max.generations       = 100,
                         scale.domains         = 5,
                         print.level = 0
)

# estimate a model, with "y3" being the reference category
# true coefficient values are:  (Intercept) = 0, x = 1
# impose an equality constraint
# equality constraint:  coefficients of x1 and x2 are equal
mulrobE <- multinomRob(list(y1 ~ x1, y2 ~ x2, y3 ~ 0),
                       dtf,
                       genoud.parms = zz.genoud.parms,
                       print.level = 3, iter=FALSE);


summary(mulrobE, weights=TRUE);

#Do only MLE estimation.  The following model is NOT identified if we
#try to estimate the overdispersed MNL.
dtf <- data.frame(y1=c(1,1),y2=c(2,1),y3=c(1,2),x=c(0,1))
dtf

summary(multinomRob(list(y1 ~ 0, y2 ~ x, y3 ~ x), data=dtf, MLEonly=TRUE))


#######################################################

#######################################################

setwd("/home/gosia/Multinomial/Simulations")

# load information about simulation 
simu_info <- read.table("Simu_info/simu12_reads_final.txt", header=T)
# write.table(simu_info, "Simu_info/simu12_reads_final.xls", sep="\t", quote=F, row.names=F, col.names=T)
head(simu_info)


# create metadata file
metadata <- data.frame(SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "R",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))


library("edgeR")
library("DEXSeq")

# run read.HTSeqCounts to create ExonCountSet object
metadata$condition <- as.factor(metadata$condition)


ecs <- read.HTSeqCounts(countfiles = paste0("MISO_counts/", metadata$SampleName, "_miso_sim12.txt"), design = metadata[,c("condition")])
sampleNames(ecs) = metadata$SampleName


counts <- counts(ecs)
gene.id <- fData(ecs)$geneID
ete.id <- fData(ecs)$exonID


gene.id.unique <- unique(gene.id)

#g = gene.id.unique[1]
g="FBgn0261573"

counts.tmp <- counts[gene.id == g, ] + 1
dim(counts.tmp)

design <- model.matrix(~condition, data=metadata)

# put data into a dataframe
dtf <- data.frame(design, t(counts.tmp));
# dtf <- rbind(dtf, dtf)
dim(dtf)

head(dtf)

regressors <- names(dtf)[1:2]
events <- names(dtf)[-(1:2)]

formulas1 <- paste0(events, "~ -1 + ", regressors[1], "+", regressors[2])
formulas0 <- paste0(events, "~ -1 + ", regressors[1])

formulas.list1 <- lapply(formulas1, as.formula)
formulas.list1[[length(formulas.list1)]] <- as.formula(paste0(events[length(events)], "~ 0"))

formulas.list0 <- lapply(formulas0, as.formula)
formulas.list0[[length(formulas.list0)]] <- as.formula(paste0(events[length(events)], "~ 0"))


mulrobE1 <- multinomRob(formulas.list1, dtf, MLEonly=TRUE)
mulrobE0 <- multinomRob(formulas.list0, dtf, MLEonly=TRUE)


summary(mulrobE1)
summary(mulrobE0)


names(mulrobE1)

















