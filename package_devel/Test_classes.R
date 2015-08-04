
# BioC 3.1
# Created 21 July 2015

##################################################################################################
### dm_fit, dmMatrixList, dmList classes 
##################################################################################################

########################################################
# on rum_dm_0_1_6 DS
########################################################

setwd("/home/gosia/multinomial_project/simulations_sim5_drosophila_noDE_noNull/")


library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

library(limma)
library(edgeR)

library(DM)


### Source all R files in DM package
Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])




count_method_list <- c("htseq", "kallisto")
count_method <- count_method_list[1]

filter_method_list <- c("filtering_dexseq", "filtering_min3prop0_01min6cpm1maxInf", "filtering_min3prop0_05min6cpm1maxInf") 
filter_method <- filter_method_list[1]

out_dir <- paste0("dm_0_1_6/", count_method, "/", filter_method, "/")

load(paste0(out_dir, "/data.RData"))

BPPARAM <- MulticoreParam(workers = 5)


# object.size(data@counts)

# library(IRanges)

# counts_partitioning <- PartitioningByEnd(cumsum(table(factor(gene_id, levels = unique(gene_id)))))

# object.size(counts_partitioning)




##### run DM pipeline with common dispersion 

mean_expression <- dm_estimateMeanExpression(data, BPPARAM = MulticoreParam(workers = 2))


dispersion <- 4000


fit <- dmDS_fit(data, dispersion, prop_mode = c("constrOptim", "constrOptimG", "FisherScoring")[2], prop_tol = 1e-12, verbose = FALSE, BPPARAM = MulticoreParam(workers = 10))


setClass("dm_fit", representation(pi = "matrix", gamma0 = "numeric", df = "numeric", lik = "numeric"))

new("dm_fit")



fit_full <- fit$fit_full
fit_null <- fit$fit_null



fit_full_s4 <- lapply(fit_full, function(f){

  if(is.null(f))
  return(NULL)
  
  new("dm_fit", pi = f$pi, gamma0 = f$gamma0, df = f$df, lik = f$lik)
  
  })

fit_null_s4 <- lapply(fit_null, function(f){

  if(is.null(f))
  return(NULL)

  new("dm_fit", pi = f$pi, gamma0 = f$gamma0, df = f$df, lik = f$lik)

  })



object.size(fit_full)
object.size(fit_full_s4)

object.size(fit_null)
object.size(fit_null_s4)


setClass("dmDS_fit", representation = representation(fit_full = "list", fit_null = "list"), contains = "dmDSdata")
fit <- new("dmDS_fit", fit_full = fit_full_s4, fit_null = fit_null_s4)




setClass("dmDS_fit2", contains = "dmDSdata", representation = representation(fit_full = "list", fit_null = "list"))

new("dmDS_fit2")








##################################################################################################
### dm_list class 
##################################################################################################


library(GenomicRanges)


x <- relist(matrix(1, 20, 3), PartitioningByEnd(seq(5, 20, 5)))
class(x) ### SimpleList

# setClassUnion("dmMatrix", "matrix")

# setClass("dmMatrix", contains = "matrix")

# setClass("dmMatrixList", contains = "CompressedList", representation(unlistData = "dmMatrix"))
# ll <- new("dmMatrixList")


setClass("dmMatrixList", contains = "CompressedList", representation(unlistData = "matrix"), prototype(elementType = "matrix"))

ll <- new("dmMatrixList")

ll <- new("dmMatrixList", unlistData = matrix(rep(1:4, each = 15), 20, 3, byrow = TRUE), partitioning = PartitioningByEnd(seq(5, 20, 5)))



dmMatrixList <- function(...){

  listData <- list(...)
  if (length(listData) == 1L && is.list(listData[[1L]]))
  listData <- listData[[1L]]
  if (length(listData) == 0L) {
    return(new("dmMatrixList"))
  } 
  else {
    if (!all(sapply(listData, is, "matrixORNULL")))
    stop("all elements in '...' must be matrix or NULL objects")

    unlistData <- do.call(rbind, listData)

    return(new("dmMatrixList", unlistData = unlistData, partitioning = PartitioningByEnd(listData)))		

  }
}

dmMatrixList(pis)




show_matrix <- function(object, nhead = 3, ntail = 3){
  nr <- nrow(object)
  nc <- ncol(object)

  cat(mode(object), " " ,class(object), " with ", nr, ifelse(nr == 1, " row and ", " rows and "), nc, ifelse(nc == 1, " column\n", " columns\n"), sep = "")

  if(nr > 0 && nc > 0){
    nms <- rownames(object)
    if(nr < (nhead + ntail + 1L)) {
      out <- object
      } else {
        out <- do.call(rbind, list(head(object, nhead), rep.int("...", nc), tail(object, ntail)))
        rownames(out) <- rownames_matrix(nms, nr, nhead, ntail) 
      }
      print(out, quote = FALSE, right = TRUE)
    }
  }


  rownames_matrix <- function(nms, nrow, nhead, ntail){
    p1 <- ifelse (nhead == 0, 0L, 1L)
    p2 <- ifelse (ntail == 0, 0L, ntail-1L)
    s1 <- s2 <- character(0)

    if (is.null(nms)) {
      if (nhead > 0) 
      s1 <- paste0("[",as.character(p1:nhead), ",]")
      if (ntail > 0) 
      s2 <- paste0("[",as.character((nrow-p2):nrow), ",]")
      } else { 
        if (nhead > 0) 
        s1 <- paste0(head(nms, nhead))
        if (ntail > 0) 
        s2 <- paste0(tail(nms, ntail))
      }
      c(s1, "...", s2)
    }


    setMethod("show", "dmMatrixList", function(object){

      nn <- length(object)
      cat("dmMatrixList of length", nn,"\n")

      if(nn){

        i <- names(object)
        if(is.null(i)) i <- seq_len(nn)

        what <- i[1]
        cat("$", what, "\n", sep="")
        show_matrix(object@unlistData[object@partitioning[[what]], , drop = FALSE])

        if(nn > 1){
          if(nn > 2)
          cat("$", "...", "\n", sep="")
          what <- i[nn]
          cat("$", what, "\n", sep="")
          show_matrix(object@unlistData[object@partitioning[[what]], , drop = FALSE])
        }

      }

      cat("with mcols \n")
      print(object@elementMetadata)


  # if(!is.null(object@elementMetadata))
  # show_matrix(object@elementMetadata, 1, 1)
  # else 
  # print(NULL)
  
  })






ll[[1]] ### gives only [1] 1 1 1 1 1

ll[1] ### works

attributes(ll)

rr <- relist(ll, PartitioningByEnd(seq(2, 20, 2))) ### does NOT work 


setMethod("[[", "dmMatrixList", function(x, i, ...){

  x@unlistData[x@partitioning[[i]], , drop = FALSE]
  
  })


ll[[1]] <- matrix(8, 8, 3)




names(ll) <- paste("name", 1:length(ll))
mcols(ll) <- DataFrame(info = rep("valid", length(ll)))
ll1 <- ll
mcols(ll1) <- NULL
tls <- unname(list(ll, ll, ll1 ))



setMethod("c", "dmMatrixList", function(x, ..., recursive = FALSE){

  if (!identical(recursive, FALSE))
  stop("\"c\" method for dmMatrixList objects ", "does not support the 'recursive' argument")
  if (missing(x))
  tls <- unname(list(...))
  else tls <- unname(list(x, ...))
  if (!all(sapply(tls, is, "dmMatrixList")))
  stop("all arguments in '...' must be dmMatrixList objects")
  ecs <- sapply(tls, elementType)
  if (!all(sapply(ecs, extends, ecs[[1L]])))
  stop("all arguments in '...' must have an element class that extends that of the first argument")
  
  
  unlistData <- do.call(rbind, lapply(tls, function(tl) tl@unlistData))
  
  ends <- cumsum(do.call(c, lapply(tls, function(tl){
    # tl <- tls[[1]]
    
    w <- width(tl@partitioning)
    if(!is.null(names(tl@partitioning)))
    names(w) <- names(tl@partitioning)
    w
    
    })))
  
  if(all(sapply(tls, function(tl) is.null(tl@elementMetadata))))
  elementMetadata <- NULL
  
  dim_meta <- unique(unlist(lapply(tls, function(tl) ncol(tl@elementMetadata))))
  
  if(!length(dim_meta) == 1){
    elementMetadata <- NULL
    }else{

      names_meta <- do.call(rbind, lapply(tls, function(tl) colnames(tl@elementMetadata)))

      if(!nrow(unique(names_meta)) == 1){
        elementMetadata <- NULL

        }else{
          elementMetadata <- do.call(rbind, lapply(tls, function(tl){
        # tl <- tls[[3]]
        if(is.null(tl@elementMetadata)){
          el_meta <- DataFrame(matrix(NA, length(tl@partitioning), dim_meta))
          colnames(el_meta) <- names_meta[1, ]
          # if(!is.null(names(tl@partitioning)))
          # 	rownames(el_meta) <- names(tl@partitioning)
          return(el_meta)
        }
        tl@elementMetadata
        }))
        }
      }

      new("dmMatrixList", unlistData = unlistData, partitioning = PartitioningByEnd(ends), elementMetadata = elementMetadata)

      })




## Construction with GRangesList():
gr1 <- GRanges(seqnames = "chr2", ranges = IRanges(3, 6),
 strand = "+", score = 5L, GC = 0.45)
gr2 <- GRanges(seqnames = c("chr1", "chr1"),
 ranges = IRanges(c(7,13), width = 3),
 strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
gr3 <- GRanges(seqnames = c("chr1", "chr2"),
 ranges = IRanges(c(1, 4), c(3, 9)),
 strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))
grl <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
grl


relist(grl, PartitioningByEnd(seq(1,5,1))) ### does NOT work eaither

c(grl, grl)

showMethods("c")

getMethod("c", c(x="CompressedList"))




#################### dm_list class with matrix as elementMetadata


# setClass("dm_list", contains = "Vector", representation(unlistData = "matrix", partitioning = "PartitioningByEnd", elementMetadata = "DataFrame"))

# x <- new("dm_list")

# x <- new("dm_list", unlistData = matrix(rep(1:4, each = 15), 20, 3, byrow = TRUE), partitioning = PartitioningByEnd(seq(5, 20, 5), names = paste0("part", 1:4)), elementMetadata = new("DataFrame", nrows=length(ends)))

# methods(class = "dm_list")

setClassUnion("matrixORNULL", c("matrix", "NULL"))

setClass("dm_list", representation(unlistData = "matrix", partitioning = "PartitioningByEnd", elementMetadata = "matrixORNULL"))

x <- new("dm_list")

x <- new("dm_list", unlistData = matrix(rep(1:4, each = 15), 20, 3, byrow = TRUE), partitioning = PartitioningByEnd(seq(5, 20, 5), names = paste0("part", 1:4)), elementMetadata = NULL)

methods(class = "dm_list")


object <- matrix(rep(1:4, each = 15), 20, 3, byrow = TRUE)


show_matrix <- function(object, nhead = 3, ntail = 3){
  nr <- nrow(object)
  nc <- ncol(object)
  
  cat(mode(object), " " ,class(object), " with ", nr, ifelse(nr == 1, " row and ", " rows and "), nc, ifelse(nc == 1, " column\n", " columns\n"), sep = "")
  
  if(nr > 0 && nc > 0){
    nms <- rownames(object)
    if(nr < (nhead + ntail + 1L)) {
      out <- object
      } else {
        out <- do.call(rbind, list(head(object, nhead), rep.int("...", nc), tail(object, ntail)))
        rownames(out) <- rownames_matrix(nms, nr, nhead, ntail) 
      }
      print(out, quote = FALSE, right = TRUE)
    }
  }


  rownames_matrix <- function(nms, nrow, nhead, ntail){
    p1 <- ifelse (nhead == 0, 0L, 1L)
    p2 <- ifelse (ntail == 0, 0L, ntail-1L)
    s1 <- s2 <- character(0)

    if (is.null(nms)) {
      if (nhead > 0) 
      s1 <- paste0("[",as.character(p1:nhead), ",]")
      if (ntail > 0) 
      s2 <- paste0("[",as.character((nrow-p2):nrow), ",]")
      } else { 
        if (nhead > 0) 
        s1 <- paste0(head(nms, nhead))
        if (ntail > 0) 
        s2 <- paste0(tail(nms, ntail))
      }
      c(s1, "...", s2)
    }


    object <- x


    setMethod("show", "dm_list", function(object){

      nn <- length(object)
      cat("dm_list of length", nn,"\n")

      if(nn){

        i <- names(x)
        if(is.null(i)) i <- seq_len(nn)

        what <- i[1]
        cat("$", what, "\n", sep="")
        show_matrix(object@unlistData[object@partitioning[[what]], , drop = FALSE])

        if(nn > 1){
          cat("$", "...", "\n", sep="")
          what <- i[nn]
          cat("$", what, "\n", sep="")
          show_matrix(object@unlistData[object@partitioning[[what]], , drop = FALSE])
        }

      }

      cat("with mcols \n")

      if(!is.null(object@elementMetadata))
      show_matrix(object@elementMetadata, 1, 1)
      else 
      print(NULL)

      })




    setMethod("[[", "dm_list", function(x, i, ...){

      x@unlistData[x@partitioning[[i]], , drop = FALSE]

      })


    setMethod("$", "dm_list", function(x, name){
      x[[name]] 
      })



    setMethod("[", "dm_list", function(x, i, ...){

      subset <- lapply(i, function(ii){
        x@unlistData[x@partitioning[[ii]], , drop = FALSE]
        })

      unlistData <- do.call(rbind, subset)
      partitioning <- PartitioningByEnd(subset, names = names(x@partitioning)[i])

      new("dm_list", unlistData = unlistData, partitioning = partitioning)

      })




    setMethod("elementLengths", "dm_list", function(x){
      ans <- elementLengths(x@partitioning)
      names(ans) <- names(x)
      ans
      })


    setMethod("length", "dm_list", function(x) length(x@partitioning))

    setMethod("names", "dm_list", function(x) names(x@partitioning))

    setReplaceMethod("names", "dm_list", function(x, value){
      names(x@partitioning) <- value
      x
      })


    setMethod("mcols", "dm_list", function (x, use.names = FALSE, ...){
      if (!isTRUEorFALSE(use.names))
      stop("'use.names' must be TRUE or FALSE")
      ans <- x@elementMetadata
      if (use.names && !is.null(ans))
      rownames(ans) <- names(x)
      ans
      })

    setReplaceMethod("mcols", "dm_list", function (x, ..., value){
      if (!is(value, "matrixORNULL"))
      stop("replacement 'elementMetadata' value must be a matrix or NULL")
      if ("elementMetadata" %in% names(attributes(x))) {
        if (!is.null(value) && length(x) != nrow(value))
        stop("supplied metadata columns must have the length of 'x'")
        if (!is.null(value))
        rownames(value) <- NULL
        x@elementMetadata <- value
      }
      x
      })



# setMethod("relist", "dm_list", function(flesh, skeleton))









########################################################
# on rum_dm_0_1_5 SQTL
########################################################



setwd("/home/Shared/data/seq/geuvadis/")

library(DM)
library(BiocParallel)

library(limma)
source("/home/gosia/R/multinomial_project/package_devel/0_my_printHead.R")

library(edgeR)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

library(pryr)

Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
for(i in Rfiles) source(i)



out.dir <- "dm_0_1_5_analysis/Results_Data_DM_0_1_5_TagwiseDisp_gridNone/"

load(paste0(out.dir, "dgeSQTL_1000.RData"))


object.size(dgeSQTL$counts)
object.size(dgeSQTL$genotypes)

object.size(dgeSQTL$fitFull)
object.size(dgeSQTL$fitNull)

object_size(dgeSQTL$counts)
object_size(dgeSQTL$genotypes)

object_size(dgeSQTL$fitFull)
object_size(dgeSQTL$fitNull)




fitFull_pi <- lapply(dgeSQTL$fitFull, function(ffg){
  # ffg <- dgeSQTL$fitFull[[1]]
  pis <- lapply(ffg, function(ffs){

    if(is.null(ffs))
    return(NULL)
    
    new_pi <- matrix(NA, nrow = nrow(ffs$pi), ncol = 3)
    colnames(new_pi) <- c("0", "1", "2")
    rownames(new_pi) <- rownames(ffs$pi)
    
    new_pi[,colnames(ffs$pi)] <- ffs$pi
    
    new_pi
    
    })
  
  do.call(rbind, pis)
  
  })

object_size(dgeSQTL$fitFull)
object_size(fitFull_pi)


fitFull_gamma0 <- dgeSQTL$tagwiseDispersion
object_size(fitFull_gamma0)



fitFull_lik <- lapply(dgeSQTL$fitFull, function(ffg){

  liks <- lapply(ffg, function(ffs){

    if(is.null(ffs))
    return(NA)
    
    new_lik <- sum(ffs$logLik)
    new_lik
    
    })
  
  names(liks) <- names(ffg)
  
  unlist(liks)
  
  })


object_size(fitFull_lik)





fitFull_meta <- lapply(dgeSQTL$fitFull, function(ffg){

  metas <- lapply(ffg, function(ffs){

    if(is.null(ffs))
    return(rep(NA, 2))
    
    new_meta <- c(ffs$gamma0, sum(ffs$logLik))
    new_meta
    
    })
  
  metas <- do.call(rbind, metas)
  colnames(metas) <- c("gamma0", "lik")
  metas
  })


object_size(fitFull_meta)






fitFull<- lapply(dgeSQTL$fitFull, function(ffg){
  # ffg <- dgeSQTL$fitFull[[1]]
  
  pis <- lapply(ffg, function(ffs){

    if(is.null(ffs))
    return(NULL)
    
    new_pi <- matrix(NA, nrow = nrow(ffs$pi), ncol = 3)
    colnames(new_pi) <- c("0", "1", "2")
    rownames(new_pi) <- rownames(ffs$pi)
    
    new_pi[,colnames(ffs$pi)] <- ffs$pi
    
    new_pi
    
    })
  
  dmFit <- dmMatrixList(pis)
  
  metas <- lapply(ffg, function(ffs){

    if(is.null(ffs))
    return(rep(NA, 2))
    
    new_meta <- c(ffs$gamma0, sum(ffs$logLik))
    new_meta
    
    })
  
  metas <- do.call(rbind, metas)
  colnames(metas) <- c("gamma0", "lik")
  
  mcols(dmFit) <- DataFrame(metas)
  
  dmFit
  
  })

object_size(dgeSQTL$fitFull)
object_size(fitFull)



DF <- mcols(fitFull[[1]])
M <- as.matrix(DF)
object.size(DF)
object.size(M)
















