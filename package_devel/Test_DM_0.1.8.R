# BioC Devel

# Created 19 Aug 2015


##############################################################################################################
# Test DM package on data from Mark
##############################################################################################################


setwd("/home/Shared/tmp/gosia/")

load("first_hurdle.Rdata")


o <- order(lookup$gene_id[m], lookup$transcript_id[m]) 

library(DM)


d <- dmDSdata(counts = as.matrix(countsw[o,]), gene_id_counts = lookup$gene_id[m][o], 
 feature_id_counts = rownames(countsw)[o], sample_id = colnames(countsw), 
 group = as.factor(df$condition))

# system.time(d <- dmDSfilter(d))


myFilter <- function(counts, samples, min_samps_gene_expr = 3, min_gene_expr = 1,
 min_samps_feature_prop = 3, min_feature_prop = 0.01, max_features = Inf)
{
  
 counts_cpm <- new("MatrixList", unlistData = DM:::dm_cpm(counts@unlistData),
   partitioning = counts@partitioning)
 gene_list <- names(counts)


# clunky
z <- as.data.frame(counts@partitioning)

inds <- vector("list",nrow(z))
for(i in 1:nrow(z))
inds[[i]] <- z$start[i]:z$end[i]
names(inds) <- z$names


n <- sapply(inds,length)

counts_new <- lapply(inds[n>1], function(g) {
  
 expr_cpm_gene <- counts_cpm@unlistData[g,,drop=FALSE]
 expr_gene <- counts@unlistData[g,,drop=FALSE]
 
 
 if (dim(expr_gene)[1] == 1)
 return(NULL)
 if (!sum(colSums(expr_cpm_gene) > min_gene_expr) >= min_samps_gene_expr)
 return(NULL)
 samps2keep <- colSums(expr_cpm_gene) != 0 & !is.na(expr_cpm_gene[1,
   ])
 prop <- prop.table(expr_gene[, samps2keep], 2)
 trans2keep <- rowSums(prop > min_feature_prop) >= min_samps_feature_prop
 if (sum(trans2keep) <= 1)
 return(NULL)
 if (!max_features == Inf) {
   if (sum(trans2keep) > max_features) {
     tr_order <- order(apply(aggregate(t(prop[trans2keep,
       ]), by = list(group = samples$group[samps2keep]),
     median)[, -1], 2, max), decreasing = TRUE)
     trans2keep <- trans2keep[trans2keep]
     trans2keep <- names(trans2keep[sort(tr_order[1:max_features])])
   }
 }
 expr <- expr_gene[trans2keep, ]
 return(expr)
 })

   #names(counts_new) <- names(counts)
   #return(counts_new)
   counts_new <- DM:::MatrixList(counts_new)
   data <- new("dmDSdata", counts = counts_new, samples = samples)
   return(data)
   
 }


system.time(d <- myFilter(d@counts,d@samples))


d_disp <- dmDSdispersion(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))








################################################################################
### Test new MatrixListFast
################################################################################

setClass("MatrixListFast",
  representation(unlistData = "matrix", 
    partitioning = "list"))


mms <- d@counts



# w <- width(d@counts@partitioning)
# inds <- split( 1:sum(w), rep(names(d@counts@partitioning),w) )


# clunky
z <- as.data.frame(d@counts@partitioning)
inds <- vector("list",nrow(z))
for(i in 1:nrow(z))
inds[[i]] <- z$start[i]:z$end[i]
names(inds) <- z$names


mmf <- new("MatrixListFast", unlistData = d@counts@unlistData, partitioning = inds)


object.size(mms)

object.size(mmf)



system.time(x <- lapply(names(mmf@partitioning), function(i){
  
  mmf@unlistData[mmf@partitioning[[i]], , drop = FALSE]
  
  }))


system.time(x <- lapply(1:length(mmf@partitioning), function(i){
  
  mmf@unlistData[mmf@partitioning[[i]], , drop = FALSE]
  
  }))


system.time(x <- lapply(mmf@partitioning, function(i){
  
  mmf@unlistData[i, , drop = FALSE]
  
  }))






access_mmf <- function(x, i){
  
  x@unlistData[x@partitioning[[i]], , drop = FALSE]
  
}


access_mms <- function(x, i){
  
  x@unlistData[x@partitioning[[i]], , drop = FALSE]
  
}





system.time(x <- lapply(names(mmf@partitioning), function(i){
  
  access_mmf(mmf, i)
  
  }))


system.time(x <- lapply(1:length(mmf@partitioning), function(i){
  
  access_mmf(mmf, i)
  
  }))



system.time(x <- lapply(names(mms@partitioning), function(i){
  
  access_mms(mms, i)
  
  }))


system.time(x <- lapply(1:length(mms@partitioning), function(i){
  
  access_mms(mms, i)
  
  }))








