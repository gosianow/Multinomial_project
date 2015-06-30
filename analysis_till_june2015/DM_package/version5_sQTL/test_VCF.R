setwd("/home/Shared/data/seq/GEUVADIS/")


library(VariantAnnotation)
library(GenomicFeatures)


## Pleas make sure exact data.zip under your current work directory
vcf <- readVcf("genotypes/GEUVADIS.chr19.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz", "hg19")
vcf


#Extract the header information stored in the exptData slot
hdr <-exptData(vcf)[["header"]]
hdr




head(rowData(vcf))


head(fixed(vcf),3)


alternate <- alt(vcf)


## number of ALT values per variant
unique(elementLengths(alternate))

geno(vcf)

geno(vcf)$DP[1:3,]


## variants filtering

index1 <-lapply(as.data.frame(info(vcf),row.names=NULL)[["DP4"]], function(x) (x[3]+x[4])>(x[1]+x[2]) & (x[3]+x[4])>3)
index2 <- as.data.frame(info(vcf),row.names=NULL)[["MQ"]]>10
index <- which(index1==TRUE & index2==TRUE)
filtered_fixed = fixed(vcf)[index,]
length(filtered_fixed)


## subsetting VCF


svp <- ScanVcfParam(fixed="ALT",info="DP4", geno="GT")
vcf_flds <- readVcf("data/var.raw.vcf.gz", "TAIR10", svp)
geno(vcf_flds)





## Subsetting VCF with Ranges

which <- RangesList(Chr1=IRanges(1, 100000),
                      Chr2=IRanges(100, 100000),
                      Chr3=IRanges(100, 100000),
                      Chr4=IRanges(1, 100000),
                      Chr5=IRanges(1, 100000))
vcf <- readVcf("data/var.raw.vcf.gz","TAIR10",ScanVcfParam(which=which))





































