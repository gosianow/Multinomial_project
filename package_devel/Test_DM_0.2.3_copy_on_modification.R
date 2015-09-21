### Run on Mac

setwd("Desktop/")

install.packages("DM_0.2.3.tar.gz", type = "source")

library(pryr)
library(DM)


counts <- as.matrix(dataDS_counts[,-1])
group_id <- dataDS_counts[,1]
group_split <- limma::strsplit2(group_id, ":")
gene_id <- group_split[, 1]
feature_id <- group_split[, 2]
sample_id = dataDS_metadata$sample_id
group = dataDS_metadata$group


mem_used()

d <- dmDSdata(counts = counts, gene_id = gene_id, feature_id = feature_id, sample_id = sample_id, group = group)
print(c(address(d), refs(d)))


tracemem(d)
mem_change(d <- dmFilter(d))


d <- dmFilter(d)
print(c(address(d), refs(d)))



tracemem(d)

d <- dmDispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 4))
print(c(address(d), refs(d)))



plotDispersion(d)


d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 1))


d <- dmLRT(d)

res <- results(d)

plotLRT(d)
plot(d)


gene_id <- res$gene_id[1]
plotFit(d, gene_id = gene_id)




