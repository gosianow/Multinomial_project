#############################################
# initialize
#############################################


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


#############################################
# SummarizedExperiment
#############################################

library(GenomicRanges)

nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

rowData <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
	IRanges(floor(runif(200, 1e5, 1e6)), width=100))

colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
	row.names=LETTERS[1:6])


sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
	rowData=rowData, colData=colData)
sset

sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
  rowData=rowData, colData=colData)
sset



assay(sset[1, ])


rowDataList <- relist(rowData, PartitioningByEnd(seq(5, 200, by = 5)))


sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
	rowData=rowDataList, colData=colData)
sset




















