
## ----echo=FALSE----------------------------------------------------------
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(cgdv17))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(PolyPhen.Hsapiens.dbSNP131))


## ------------------------------------------------------------------------
library(VariantAnnotation)
library(cgdv17)
file <- system.file("vcf", "NA06985_17.vcf.gz", package = "cgdv17")

## Explore the file header with scanVcfHeader
hdr <- scanVcfHeader(file)

info(hdr) 

geno(hdr) 


## ------------------------------------------------------------------------
## get entrez ids from gene symbols
library(org.Hs.eg.db)
genesym <- c("TRPV1", "TRPV2", "TRPV3")
geneid <- select(org.Hs.eg.db, keys=genesym, keytype="SYMBOL",
                 columns="ENTREZID")
geneid


## ------------------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb


## ------------------------------------------------------------------------
txdb <- renameSeqlevels(txdb, gsub("chr", "", seqlevels(txdb)))
txdb <- keepSeqlevels(txdb, "17")


## ------------------------------------------------------------------------
txbygene = transcriptsBy(txdb, "gene")


## ------------------------------------------------------------------------
gnrng <- unlist(range(txbygene[geneid$ENTREZID]), use.names=FALSE)
names(gnrng) <- geneid$SYMBOL


## ------------------------------------------------------------------------
param <- ScanVcfParam(which = gnrng, info = "DP", geno = c("GT", "cPd"))
param

## Extract the TRPV ranges from the VCF file 
vcf <- readVcf(file, "hg19", param)
## Inspect the VCF object with the 'fixed', 'info' and 'geno' accessors
vcf

head(fixed(vcf))
dim(fixed(vcf))
geno(vcf)


## ----, eval=FALSE--------------------------------------------------------
## ## Use the 'region' argument to define the region
## ## of interest. See ?locateVariants for details.
## cds <- locateVariants(vcf, txdb, CodingVariants())
## five <- locateVariants(vcf, txdb, FiveUTRVariants())
## splice <- locateVariants(vcf, txdb, SpliceSiteVariants())
## intron <- locateVariants(vcf, txdb, IntronVariants())


## ------------------------------------------------------------------------
all <- locateVariants(vcf, txdb, AllVariants())


## ------------------------------------------------------------------------
## Did any variants match more than one gene?
table(sapply(split(mcols(all)$GENEID, mcols(all)$QUERYID), 
             function(x) length(unique(x)) > 1))

## Summarize the number of variants by gene:
idx <- sapply(split(mcols(all)$QUERYID, mcols(all)$GENEID), unique)
sapply(idx, length)

## Summarize variant location by gene:
sapply(names(idx), 
       function(nm) {
         d <- all[mcols(all)$GENEID %in% nm, c("QUERYID", "LOCATION")]
         table(mcols(d)$LOCATION[duplicated(d) == FALSE])
       })














## ------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
seqlevelsStyle(vcf) <- "UCSC"
seqlevelsStyle(txdb) <- "UCSC"
aa <- predictCoding(vcf, txdb, Hsapiens)


## ------------------------------------------------------------------------
## Did any variants match more than one gene?
table(sapply(split(mcols(aa)$GENEID, mcols(aa)$QUERYID), 
             function(x) length(unique(x)) > 1))

## Summarize the number of variants by gene:
idx <- sapply(split(mcols(aa)$QUERYID, mcols(aa)$GENEID, drop=TRUE), unique)
sapply(idx, length)

## Summarize variant consequence by gene:
sapply(names(idx), 
       function(nm) {
         d <- aa[mcols(aa)$GENEID %in% nm, c("QUERYID","CONSEQUENCE")]
         table(mcols(d)$CONSEQUENCE[duplicated(d) == FALSE])
       })


## ------------------------------------------------------------------------
## Load the PolyPhen package and explore the available keys and columns
library(PolyPhen.Hsapiens.dbSNP131)
keys <- keys(PolyPhen.Hsapiens.dbSNP131)
cols <- columns(PolyPhen.Hsapiens.dbSNP131)
## column descriptions are found at ?PolyPhenDbColumns
cols

## Get the rsids for the non-synonymous variants from the
## predictCoding results
rsid <- unique(names(aa)[mcols(aa)$CONSEQUENCE == "nonsynonymous"]) 

## Retrieve predictions for non-synonymous variants. Two of the six variants 
## are found in the PolyPhen database. 
select(PolyPhen.Hsapiens.dbSNP131, keys=rsid, 
       columns=c("AA1", "AA2", "PREDICTION"))


## ----eval=FALSE----------------------------------------------------------
## library(BiocInstaller)
## biocLite("VariantAnnotation")


## ----eval=FALSE----------------------------------------------------------
## library(VariantAnnotation)


## ----eval=FALSE----------------------------------------------------------
## browseVignettes(package="VariantAnnotation")


## ----eval=FALSE----------------------------------------------------------
## help.start()


## ------------------------------------------------------------------------
sessionInfo()

