## ---- echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages(library(VariantAnnotation))
suppressMessages(suppressPackageStartupMessages(library(cgdv17)))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(PolyPhen.Hsapiens.dbSNP131))

## ---- eval=FALSE---------------------------------------------------------
## library(VariantAnnotation)
## library(cgdv17)
## library(org.Hs.eg.db)
## library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## library(BSgenome.Hsapiens.UCSC.hg19)
## library(PolyPhen.Hsapiens.dbSNP131)

## ---- eval=FALSE---------------------------------------------------------
## source("http://bioconductor.org/biocLite.R")
## biocLite("mypackage")

## ---- eval=FALSE---------------------------------------------------------
## browseVignettes("cgdv17")

## ------------------------------------------------------------------------
library(VariantAnnotation)
library(cgdv17)
file <- system.file("vcf", "NA06985_17.vcf.gz", package = "cgdv17")

## ------------------------------------------------------------------------
hdr <- scanVcfHeader(file)
 
info(hdr) 
 
geno(hdr) 

## ------------------------------------------------------------------------
meta(hdr)$META

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

geno(vcf)

## ---- eval=FALSE---------------------------------------------------------
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

## ---- echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages(library(ensemblVEP))

## ------------------------------------------------------------------------
library(ensemblVEP)

## ------------------------------------------------------------------------
dest <- tempfile()
writeVcf(vcf, dest)

## ---- eval=FALSE---------------------------------------------------------
## gr <- ensemblVEP(file = dest)

## ---- eval=FALSE---------------------------------------------------------
## 2015-06-09 08:50:22 - Starting...
## 2015-06-09 08:50:24 - Detected format of input file as vcf
## 2015-06-09 08:50:24 - Read 277 variants into buffer
## 2015-06-09 08:50:24 - Reading transcript data from cache and/or database
## [================================================================================================================================================================]  [ 100% ]
## 2015-06-09 08:54:29 - Retrieved 55 transcripts (0 mem, 0 cached, 55 DB, 0 duplicates)
## 2015-06-09 08:54:29 - Analyzing chromosome 17
## 2015-06-09 08:54:29 - Analyzing variants
## [================================================================================================================================================================]  [ 100% ]
## 2015-06-09 08:54:29 - Calculating consequences
## [================================================================================================================================================================]  [ 100% ]
## 2015-06-09 08:54:30 - Processed 277 total variants (1 vars/sec, 1 vars/sec total)
## 2015-06-09 08:54:30 - Wrote stats summary to /tmp/RtmppUYah1/file339f64bf5700_summary.html
## 2015-06-09 08:54:30 - Finished!

## ---- eval=FALSE---------------------------------------------------------
## head(gr, 3)

## ---- eval=FALSE---------------------------------------------------------
## GRanges object with 3 ranges and 22 metadata columns:
##            seqnames             ranges strand |   Allele
##               <Rle>          <IRanges>  <Rle> | <factor>
##   rs402369    chr17 [3469401, 3469401]      * |        G
##   rs402369    chr17 [3469401, 3469401]      * |        G
##   rs402369    chr17 [3469401, 3469401]      * |        G
##                                     Consequence   IMPACT   SYMBOL
##                                        <factor> <factor> <factor>
##   rs402369 splice_region_variant&intron_variant      LOW  SPATA22
##   rs402369 splice_region_variant&intron_variant      LOW  SPATA22
##   rs402369                upstream_gene_variant MODIFIER     ASPA
##                       Gene Feature_type         Feature        BIOTYPE     EXON
##                   <factor>     <factor>        <factor>       <factor> <factor>
##   rs402369 ENSG00000141255   Transcript ENST00000397168 protein_coding     <NA>
##   rs402369 ENSG00000141255   Transcript ENST00000355380 protein_coding     <NA>
##   rs402369 ENSG00000108381   Transcript ENST00000456349 protein_coding     <NA>
##              INTRON    HGVSc    HGVSp cDNA_position CDS_position
##            <factor> <factor> <factor>      <factor>     <factor>
##   rs402369      1/8     <NA>     <NA>          <NA>         <NA>
##   rs402369      1/7     <NA>     <NA>          <NA>         <NA>
##   rs402369     <NA>     <NA>     <NA>          <NA>         <NA>
##            Protein_position Amino_acids   Codons Existing_variation DISTANCE
##                    <factor>    <factor> <factor>           <factor> <factor>
##   rs402369             <NA>        <NA>     <NA>               <NA>     <NA>
##   rs402369             <NA>        <NA>     <NA>               <NA>     <NA>
##   rs402369             <NA>        <NA>     <NA>               <NA>     4652
##              STRAND SYMBOL_SOURCE    HGNC_ID
##            <factor>      <factor>   <factor>
##   rs402369       -1          HGNC HGNC:30705
##   rs402369       -1          HGNC HGNC:30705
##   rs402369        1          HGNC   HGNC:756
##   -------
##   seqinfo: 25 sequences from  genome; no seqlengths

## ------------------------------------------------------------------------
VEPParam()

## ------------------------------------------------------------------------
basic(VEPParam())

## ----eval=FALSE----------------------------------------------------------
## browseVignettes(package="VariantAnnotation")

## ----eval=FALSE----------------------------------------------------------
## help.start()

## ------------------------------------------------------------------------
sessionInfo()
