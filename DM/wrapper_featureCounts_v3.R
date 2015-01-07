# BioC14
# Created 15 July 2014

###########################################################################
# flattenGtf (flatten_gtf_v2.R)
###########################################################################

### flatten the annotation with disjoin()
library(GenomicRanges)
library(rtracklayer)


# gtfFile <- "/home/Shared_penticton/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf"
# gtfFile.flat <- "/home/Shared_penticton/data/annotation/Drosophila/Ensembl70/Rsubread/Drosophila_melanogaster.BDGP5.70_flat.gtf"


flattenGtf <- function(gtfFile, gtfFile.flat, ignore.strand=FALSE){ 
  
  gtf0 <- import(gtfFile)
  
  idx <- mcols(gtf0)$type == "exon" 
  gtf0 <- gtf0[idx]
  
  # split by gene
  gtf.g <- split(gtf0, mcols(gtf0)$gene_id)
  
  # flatten per gene
  gtf.g.flat <- disjoin(gtf.g)
  
  gtf.g.flat <- unlist(gtf.g.flat)
  mcols(gtf.g.flat)$gene_id <- names(gtf.g.flat)
  
  # flatten overall
  gtf.o.flat <- disjoin(gtf.g.flat,  ignore.strand=ignore.strand)
  
  # remove overlapping parts
  co <- countOverlaps(gtf.o.flat, gtf.g.flat)
  gtf.o.flat <- gtf.o.flat[co==1]
  
  # find the gene names
  fo <- findOverlaps(gtf.o.flat, gtf.g.flat)
  mcols(gtf.o.flat)$gene_id <- mcols(gtf.g.flat)$gene_id[subjectHits(fo)]
  mcols(gtf.o.flat)$transcript_id <- mcols(gtf.o.flat)$gene_id
  mcols(gtf.o.flat)$type <- "exon"
  
  # add exon id
  #   exon.id <- split(gtf.o.flat, mcols(gtf.o.flat)$gene_id)
  # 
  #   exon.id.new <- unlist(lapply(exon.id, function(g){ 
  #     # g=exon.id[[1]]
  #     gg <- mcols(g)$gene_id
  #     mcols(g)$exon_id <- paste0(gg, ":", sprintf( "%03d", seq(1, length(gg) ) ) )
  #     return(g)       
  #   } ) )
  
  
  gtf.o.flat.sort <- gtf.o.flat[order(mcols(gtf.o.flat)$gene_id, decreasing=FALSE)]
  
  exon.id <- split(mcols(gtf.o.flat.sort)$gene_id, mcols(gtf.o.flat.sort)$gene_id)
  
  exon.id.new <- unlist(lapply(exon.id, function(g){ 
    # g=exon.id[[1]]
    seq(1, length(g) ) 
    
  } ) )
  
  mcols(gtf.o.flat.sort)$exon_number <- exon.id.new
  
  mcols(gtf.o.flat.sort)$exon_id <- paste0(mcols(gtf.o.flat.sort)$gene_id, ":", sprintf( "%03d", exon.id.new ) )
  
  
  export(gtf.o.flat.sort, gtfFile.flat)
  
  invisible(gtf.o.flat.sort)
  
}


###########################################################################
# wrapper_featureCounts (wrapper_featureCounts.R)
###########################################################################


# gtfFile <- "/home/Shared_penticton/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf"
# out.dir <- "Rsubread_1.14.0_test"
# 
# SE_bam_files <- with(metadata[metadata$LibraryLayout == "SINGLE", ], paste0("sam/", SampleName, "_s.bam"))
# SE_bam_files <- SE_bam_files[1]
# 
# PE_bam_files <- with(metadata[metadata$LibraryLayout == "PAIRED", ], paste0("sam/", SampleName, "_s.bam"))
# PE_bam_files <- PE_bam_files[1]

library(Rsubread)

wrapper_featureCounts <- function(gtfFile, SE_bam_files=NULL, PE_bam_files=NULL,  nthreads=10){
  
  dir.create(out.dir, showWarnings=FALSE)
  
  ### flatten GTF fully (between genes from diff strands)
  
  base.gtfFile <- basename(gtfFile)
  dir.gtfFile <- dirname(gtfFile)
  
  gtfFile.flat <- paste0(dir.gtfFile, "/Rsubread_flat_", base.gtfFile)
  gtfFile.flat.fully <- paste0(dir.gtfFile, "/Rsubread_flat_fully_", base.gtfFile)
  
  gtf.ff <- flattenGtf(gtfFile, gtfFile.flat.fully, ignore.strand=TRUE)
  
  
  ### count with featureCounts
  
  if(!is.null(SE_bam_files)){
    fc_SE <- featureCounts(SE_bam_files, annot.ext=gtfFile.flat.fully, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="exon_id", useMetaFeatures=FALSE, allowMultiOverlap=TRUE, nthreads=nthreads)
    fc <- fc_SE$counts
  }
  
  if(!is.null(PE_bam_files)){
    fc_PE <- featureCounts(PE_bam_files, annot.ext=gtfFile.flat.fully ,isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="exon_id", useMetaFeatures=FALSE, allowMultiOverlap=TRUE, isPairedEnd=TRUE, nthreads=nthreads)
    fc <- fc_PE$counts
  }
  
  if(!is.null(SE_bam_files) && !is.null(PE_bam_files))
    fc <- cbind(fc_SE$counts,fc_PE$counts)
  
  
  #### count like DEXSeq (overlapping genes from diff strands)
  
  gtf.f <- flattenGtf(gtfFile, gtfFile.flat, ignore.strand=FALSE)
  
  fo <- findOverlaps(gtf.f, gtf.ff)
  x <- queryHits(fo)
  # find exons that should be added 
  duplic.idx <- duplicated(x, fromLast=FALSE) |  duplicated(x, fromLast=TRUE)
  ex <- split(as.data.frame(fo)[duplic.idx, ],  x[duplic.idx] )
  
  
  total.count <- fc
  # add the counts
  cum.counts.tmp <- lapply(ex, function(e){
    # e = ex[[1]]
    colSums(total.count[e[,"subjectHits"], , drop=FALSE])
  } )
  
  cum.counts <- do.call( rbind, cum.counts.tmp)
  new.counts <- data.frame(exon_id=mcols(gtf.f)$exon_id)
  new.counts[, colnames(cum.counts)] <- 0
  
  fo.df <- as.data.frame(fo)
  fo.df.u <- fo.df[!duplic.idx, ]
  
  new.counts[fo.df.u[, 1], -1 ] <- total.count[fo.df.u[, 2], ]
  new.counts[rownames(cum.counts), -1 ] <- cum.counts
  
  write.table(new.counts, paste0(out.dir, "/fc.txt"), quote=F, sep="\t", row.names=F, col.names=T)
  
  annotation <- data.frame(chr=seqnames(gtf.f), start=start(gtf.f), end=end(gtf.f), name=mcols(gtf.f)$exon_id)
  
  write.table(annotation, paste0(out.dir, "/annotation.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  
  
  invisible(list(fc=new.counts, annotation=annotation))
  
}
























