# BioC 3.1
# Created 24 July 2015

##################################################################################################
### MulticoreParam vs SnowParam
##################################################################################################
### Taupo:

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





mean_expression <- dm_estimateMeanExpression(data, BPPARAM = MulticoreParam(workers = 3))

mean_expression <- dm_estimateMeanExpression(data, BPPARAM = SnowParam(workers = 2, log = TRUE))


#############################################
# initialize
#############################################
### Mac:

setClass("DotDash", representation = representation(dot = "raw", dash = "integer"))

aDotDash <- new("DotDash", dot = as.raw(1), dash = 1:1e8)
tracemem(aDotDash)
aDotDash@dot <- as.raw(2)

methods(class = "DotDash")
showMethods("initialize")

getMethod("initialize", "DotDash")


aDotDash <- new("DotDash", dot = as.raw(1), dash = 1:1e8)
tracemem(aDotDash)
aDotDash <- initialize(aDotDash, dot = as.raw(2))


aDotDash <- new("DotDash", dot = as.raw(1), dash = 1:1e8)
tracemem(aDotDash)
aDotDash <- new("DotDash", dot = as.raw(2), dash = 1:1e8)



##################################################################################################
### size of counts and genotypes in dm_0_1_5_data
##################################################################################################
### Mac:

library(DM)

library(BiocParallel)
library(limma)
source("/Users/gosia/R/multinomial_project/package_devel/0_my_printHead.R")
library(edgeR)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

# Rfiles <- list.files("/Users/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
# for(i in Rfiles) source(i)

library(pryr)



setwd("/Users/gosia/geuvadis/")

out.dir <- "dm_0_1_5_analysis/Results_Data_DM_0_1_5_TagwiseDisp_gridNone/"

load(paste0(out.dir, "dgeSQTL.RData"))

length(dgeSQTL$counts)
length(dgeSQTL$fitFull)

object.size(dgeSQTL)
object_size(dgeSQTL)

object.size(dgeSQTL$counts)
object.size(dgeSQTL$genotypes)
object.size(dgeSQTL$fitFull)


data <- new("dmSQTLdata", counts = dgeSQTL$counts, genotypes = dgeSQTL$genotypes, samples = data.frame(sample_id = dgeSQTL$samples$sample_names))


xx <- do.call(rbind, data@counts)
object.size(data@counts)
object.size(xx)


yy <- do.call(rbind, data@genotypes)
object.size(data@genotypes)
object.size(yy)

######## size of fit when "compressed"

load(paste0(out.dir, "dgeSQTL_1000.RData"))

object.size(dgeSQTL$counts)
object.size(dgeSQTL$genotypes)
object.size(dgeSQTL$fitFull)


### pi

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


object_size(fitFull_pi)

### gamma0

fitFull_gamma0 <- dgeSQTL$tagwiseDispersion
object_size(fitFull_gamma0)


### lik

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



### gamma0 + lik

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




object_size(dgeSQTL$fitFull)

object_size(fitFull_pi) + object_size(fitFull_gamma0) + object_size(fitFull_lik)

object_size(fitFull_pi) + object_size(fitFull_meta) 




##################################################################################################
### SummarizedExperiment
##################################################################################################

library(GenomicRanges)

nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowData <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                   IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                   strand=sample(c("+", "-"), 200, TRUE))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
               rowData=rowData, colData=colData)
sset


# sset@


##################################################################################################
### CompressedList - VIRTUAL
##################################################################################################

### IRanges
setClass("CompressedList",
         contains="List",
         representation(
                        "VIRTUAL",
                        partitioning="PartitioningByEnd",
                        unlistData = "ANY"
                       )
         )


### S4Vectors
setClass("List",
    contains="Vector",
    representation(
        "VIRTUAL",
        elementType="character"
    ),
    prototype(elementType="ANY")
)


### S4Vectors
setClassUnion("DataTableORNULL", c("DataTable", "NULL"))

setClass("Vector",
    contains="Annotated",
    representation(
        "VIRTUAL",
        elementMetadata="DataTableORNULL"
    )
)


##################################################################################################
### PartitioningByEnd
##################################################################################################

### IRanges

### Constructor:

PartitioningByEnd(x=integer(), NG=NULL, names=NULL)

# ‘x’ must be either a list-like object or a sorted integer vector.  ‘NG’ must be either ‘NULL’ or a single integer.  ‘names’ must be either ‘NULL’ or a character vector of length ‘NG’ (if supplied) or ‘length(x)’ (if ‘NG’ is not supplied).

partitioning <- PartitioningByEnd(seq(5, 25, 5))
partitioning[[1]]


##################################################################################################
### GRangesList
##################################################################################################

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


### GenomicRanges
setClass("GRangesList",
    contains=c("CompressedList", "GenomicRangesList"),
    representation(
        unlistData="GRanges",
        elementMetadata="DataFrame"
    ),
    prototype(
        elementType="GRanges"
    )
)



# grl@

attributes(grl)

methods(class = "GRangesList")

### shows a method that does not work...
relist(grl, PartitioningByEnd(seq(1,5,1))) ### does NOT work eaither

c(grl, grl)

showMethods("c")

getMethod("c", c(x="CompressedList"))



##################################################################################################
### relist() 
##################################################################################################

showMethods("relist")

x <- relist(rep(1:4, each = 5), PartitioningByEnd(seq(5, 20, 5)))
x

x <- relist(matrix(1, 20, 3), PartitioningByEnd(seq(5, 20, 5)))
x
class(x) ### SimpleList

x <- relist(matrix(1, 22, 3), PartitioningByEnd(seq(5, 20, 5)))
x
class(x) ### SimpleList


x <- relist(data.frame(matrix(1, 20, 3)), PartitioningByEnd(seq(5, 20, 5)))
x
class(x)

x <- relist(DataFrame(matrix(1, 20, 3)), PartitioningByEnd(seq(5, 20, 5)))
x
class(x) 



##################################################################################################
### dmMatrixList 
##################################################################################################
library(GenomicRanges)


setClass("dmMatrixList", contains = "CompressedList", representation(unlistData = "matrix"), prototype(elementType = "matrix"))

ll <- new("dmMatrixList")
ll
ll <- new("dmMatrixList", unlistData = matrix(rep(1:4, each = 15), 20, 3, byrow = TRUE), partitioning = PartitioningByEnd(seq(5, 20, 5)))
ll

ll[[1]] ### gives only [1] 1 1 1 1 1

ll[1] ### works

attributes(ll)

rr <- relist(ll, PartitioningByEnd(seq(2, 20, 2))) ### does NOT work 

### [[

setMethod("[[", "dmMatrixList", function(x, i, ...){

  x@unlistData[x@partitioning[[i]], , drop = FALSE]

})


ll[[1]] <- matrix(8, 8, 3)



ll[[1]] ### works!



### show


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


ll

lll <- c(ll, ll)
class(lll)


### c

names(ll) <- paste("name", 1:length(ll))
mcols(ll) <- DataFrame(info = rep("valid", length(ll)))
ll1 <- ll
mcols(ll1) <- NULL
# tls <- unname(list(ll, ll, ll1 ))



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


c(ll, ll, ll1)


dmMatrixList <- function(...){

  listData <- list(...)
  if (length(listData) == 1L && is.list(listData[[1L]]))
      listData <- listData[[1L]]
  if (length(listData) == 0L) {
      return(new("dmMatrixList"))
  } else {
      if (!all(sapply(listData, is, "matrixORNULL")))
          stop("all elements in '...' must be matrix or NULL objects")
      
      unlistData <- do.call(rbind, listData)
      
      return(new("dmMatrixList", unlistData = unlistData, partitioning = PartitioningByEnd(listData)))    

  }
}


##################################################################
# apply dmMatrixList to dgeSQTL$fitFull
##################################################################



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


##################################################################
# merge multiple data.frames
##################################################################



setwd("/home/gosia/multinomial_project/simulations_sim5_drosophila_noDE_noNull/")

# setwd("/Users/gosia/multinomial_project/simulations_sim5_drosophila_noDE_noNull/")


library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

library(limma)
library(edgeR)

library(DM)


### Source all R files in DM package
Rfiles <- list.files("/home/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
# Rfiles <- list.files("/Users/gosia/R/multinomial_project/package_devel/DM/R/", full.names=TRUE)
for(i in 1:length(Rfiles)) source(Rfiles[i])


########################################################
# load metadata
########################################################

# create metadata file
metadata <- data.frame(sample_id = paste0("sample_",1:6), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


##############################################################################################################
# htseq counts
##############################################################################################################

count_method <- "htseq"

out_dir <- paste0("dm_0_1_6/", count_method, "/")
dir.create(out_dir, showWarnings=F, recursive=T)


### load htseq counts
htseqList <- lapply(1:6, function(i){
  # i = 1
  htseq <- read.table(paste0("2_counts/dexseq_nomerge/dexseq", i, ".txt"), header = FALSE, as.is = TRUE)
  colnames(htseq) <- c("group_id", paste0("sample_", i))  
  return(htseq)
})

htseq <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), htseqList)
























