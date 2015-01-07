
########################################################
# Casper junction counting
########################################################

library(casper)

library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)

setwd("/home/gosia/Multinomial_project/SimulationsV2/CasperCounts")

#### find junctions

# gtf must be with chr
gtfFile <- "/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.chr.gtf"


genDB <- makeTranscriptDbFromGFF(file=gtfFile, format="gtf", exonRankAttributeName="exon_number")
genomeDB <- procGenome(genDB=genDB, genome="Dm70", mc.cores=10)

name=""
mc.cores=10


junctionIDs <- function(gtfFile, genomeDB, name="", mc.cores=10){
  
  # gtf
  gtf <- import(gtfFile)
  idx <- mcols(gtf)$type == "exon"
  gtf0 <- gtf[idx]
  
  
  # annotation
  annotation <- unlist(genomeDB@islands)
  
  exonsIDs <- do.call(rbind, (strsplit(names(annotation), ".", fixed=TRUE)))
  annotation$gene <- exonsIDs[,1]
  annotation$exon <- as.numeric(exonsIDs[,2])
  names(annotation) <- paste0("g", annotation$gene, "e", annotation$exon)
  
  write.table( data.frame(chr=seqnames(annotation), start= start(annotation), end=end(annotation), name = names(annotation) ), paste0("Casper_annotation",name,".bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  # find gene number to transcript
  casper.gene.ids <- do.call(rbind, (strsplit(  names(unlist(genomeDB@transcripts, recursive=FALSE)), ".", fixed=TRUE)))
  
  
  # find gaps
  gtf0t <- split(gtf0, mcols(gtf0)$transcript_id)
  
  
  gaps.tmp <- mclapply(gtf0t, function(t){
    # t = gtf0t[[6]]
    
    if(length(t)==1)
      return(NULL)
    
    gps <- gaps(t)
    gps <- gps[start(gps) != 1 ,]
    
    gps$transcript_id <- t$transcript_id[1]
    gps$gene_id <- t$gene_id[1]
    
    return(gps)
    
  }, mc.cores=mc.cores)
  
  
  gaps.tmp <- gaps.tmp[!sapply(gaps.tmp, is.null)]
  
  gaps.tmp <- GRangesList(gaps.tmp)
  
  gtf0t.gaps <- BiocGenerics::unlist(gaps.tmp)
  
  
  # find out to which genes corrresponds the numbers
  gtf0t.gaps$gene <- as.character(merge(gtf0t.gaps$transcript_id, casper.gene.ids, all.x=TRUE, by.x=1, by.y=2)$V1)
  
  names(gtf0t.gaps) <- NULL
  gtf0t.gaps <- sort(unique(gtf0t.gaps))
  
  
  names(gtf0t.gaps) <- paste0("g",gtf0t.gaps$gene, "j" , seq(1:length(gtf0t.gaps)))
  gtf0t.gaps$junction <- seq(1:length(gtf0t.gaps))
  
  start(gtf0t.gaps) <- start(gtf0t.gaps) -1
  end(gtf0t.gaps) <- end(gtf0t.gaps) +1
  
  
  write.table( data.frame(chr=seqnames(gtf0t.gaps), start=start(gtf0t.gaps), end=end(gtf0t.gaps), name = names(gtf0t.gaps)), paste0("Casper_gaps",name,".bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  
  # find junction casperIDs per gene
  
  gtf0t.gaps.splt <- split(gtf0t.gaps, gtf0t.gaps$gene)
  annotation.splt <- split(annotation, annotation$gene)
  
  
  
  junction.ids <- mclapply(names(gtf0t.gaps.splt), function(id){
    # id=2
    gtf0t.gaps <- gtf0t.gaps.splt[[id]]
    annotation <- annotation.splt[[id]]
    
    gaps.overlp <- findOverlaps(gtf0t.gaps, annotation, type="any")
    
    gaps.overlp.spl <- split(as.data.frame(gaps.overlp), queryHits(gaps.overlp))
    
    junction.ids <-  sapply(gaps.overlp.spl, function(j){
      # j=gaps.overlp.spl[[8]]
      
      junction.id <- paste0(min(annotation$exon[j$subjectHits]), ".", max(annotation$exon[j$subjectHits]))
      
      names(junction.id) <- paste0( "-g", id, "j", gtf0t.gaps$junction[j$queryHits[1]])
      
      return(junction.id)    
    })
    
    return(junction.ids)  
  }, mc.cores=mc.cores)
  
  
  names(junction.ids) <- names(gtf0t.gaps.splt)
  
  
  junction.info <- data.frame(junction_id=names(gtf0t.gaps), gene_id=gtf0t.gaps$gene_id, gene=gtf0t.gaps$gene, junction=gtf0t.gaps$junction)
  
  casper_path <- unlist(junction.ids)
  
  casper_path.tmp <- data.frame(junction_id= do.call(rbind, strsplit(names(casper_path), ".-", fixed=TRUE))[,2], casper_path=casper_path, row.names=NULL)
  
  junction.info <- merge(junction.info, casper_path.tmp)
  
  return(invisible(list(junction.ids=junction.ids, junction.info=junction.info)))
  
}


system.time( junctions <- junctionIDs(gtfFile, genomeDB, name="", mc.cores=10) ) 


save(junctions, file=paste0("Casper_junctions", name, ".RData"))

load(paste0("Casper_junctions", name, ".RData"))


#### load bams (sorted & indexed)


library(Rsamtools)
what <- scanBamWhat()
what <- what[!(what %in% c("seq","qual"))]
flag <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE)
param <- ScanBamParam(flag=flag, what=what, tag="XS")


bamFiles <- paste0("/home/Shared_taupo/tmp/Simulation/simulation_drosophila_V2/sam/sim_samp", 1:6, "_s.bam")
names(bamFiles) <- c(paste0("C1", "S", 1:3), paste0("C2", "S", 1:3))



pcs <- lapply(bamFiles, function(bamFile){
  
  bam0 <- scanBam(file=bamFile, param=param)[[1]]
  
  bam0 <- rmShortInserts(bam0, isizeMin=50)
  
  pbam0 <- procBam(bam0)
  
  ## path counting
  pc <- pathCounts(pbam0, DB=genomeDB, mc.cores = mc.cores, verbose=T)
  
  cat(paste("Done for \n", bamFile, "\n"))
  
  return(pc)
  
})


#### count junctions


junctionCounts <- function(pcs, junctions, name="", mc.cores=10){
  
  junction_counts_final <- junctions$junction.info
  samples <- names(pcs)
  
  for(i in 1:length(pcs)){
    # i=1
    pc <- pcs[[i]]
    pcc <- pc@counts[[1]]
    
    junction.ids <- junctions$junction.ids 
    
    junction_counts.tmp <- unlist(mclapply(names(junction.ids), function(id){
      # id=2
      pcc.g <- pcc[[id]]
      junction.ids.g <- junction.ids[[id]]
      
      if(is.null(pcc.g)){
        
        junction_counts <- rep(0, length(junction.ids.g))
        names(junction_counts) <- names(junction.ids.g)
        return(junction_counts)
        
      }
      
      
      paths <- do.call(rbind, strsplit(names(pcc.g), "-", fixed=TRUE))
      paths <- data.frame(paths=c(paths[,1], paths[,2]), counts=c(pcc.g, pcc.g))
      
      
      junction_counts <- sapply(junction.ids.g, function(j){
        # j=junction.ids.g[1]
        
        sum(paths$counts[grepl(j, paths$paths)])
        
      })
      
      return(junction_counts)
      
    }, mc.cores=mc.cores) )
    
    
    junction_counts <- data.frame( junction_id = do.call(rbind, strsplit(names(junction_counts.tmp), ".-", fixed=TRUE))[,2] , row.names=NULL) 
    
    junction_counts[,paste0("counts_", samples[i])] <- junction_counts.tmp
    
    junction_counts_final <- merge(junction_counts_final, junction_counts)
    
    
  }
  
  junction_counts_final$casper_path <- paste0("p", junction_counts_final$casper_path)
  
  junction_counts_final <- junction_counts_final[order(junction_counts_final$junction),]
  
  write.table(junction_counts_final, paste0("Casper_junction_counts",name,".xls"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  write.table(junction_counts_final, paste0("Casper_junction_counts",name,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  
  return(invisible(junction_counts_final))
  
}



system.time(  jc <- junctionCounts(pcs, junctions, name="", mc.cores=10)  )


write.table(jc, paste0("Casper_junction_counts",name,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



########################################################
# 
########################################################





































