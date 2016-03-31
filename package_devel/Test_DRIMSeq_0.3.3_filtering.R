
library(devtools)

load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")

samples <- data.frame(sample_id = factor(c("c1", "c2")), group = factor(c("c1", "c2")))

mm <- matrix(c(1, rep(0, 5)), nrow = 3)
mm <- matrix(c(1, 1, rep(0, 4)), nrow = 3, byrow = FALSE)

rownames(mm) <- paste0("r", 1:nrow(mm))


counts <- DRIMSeq:::MatrixList(mm = mm)
counts


d <- DRIMSeq:::dmDS_filter(counts = counts, samples = samples, min_samps_gene_expr = 0, min_gene_expr = 0, min_samps_feature_expr = 0, min_feature_expr = 0, min_samps_feature_prop = 0, min_feature_prop = 0, max_features = Inf)

d@counts









