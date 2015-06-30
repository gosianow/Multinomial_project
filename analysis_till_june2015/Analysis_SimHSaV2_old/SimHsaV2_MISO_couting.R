#######################################################
# MISO counts

# Created 30 July 2014
# Updated 30 July 2014
#######################################################
# BioC 2.14


setwd("/home/gosia/Multinomial_project/Simulations_hsapiens_V2/")


# create metadata file
metadata <- data.frame(SampleName1 = paste0(1:6), SampleName= paste0(c(rep("C1", 3), rep("C2", 3)), "S",c(1:3, 1:3)), condition=c(rep("C1", 3), rep("C2", 3)))
metadata$condition <- as.factor(metadata$condition)

metadata


# convert gtf to gff3 no chr
# /usr/bin/perl /home/gosia/Programs/gtf2gff3/ensembl_gtf_to_gff.pl Homo_sapiens.GRCh37.71_reduced.gtf > Homo_sapiens.GRCh37.71_reduced.gff3
# add chr
# awk ‘{print “chr”$0}’ Homo_sapiens.GRCh37.71_reduced.gff3 > Homo_sapiens.GRCh37.71_reduced_chr.gff3


miso.path <- "/usr/local/lib/python2.7/dist-packages/misopy-0.4.9-py2.7-linux-x86_64.egg/misopy/"


gff3 <- "/home/Shared/data/annotation/Human_tmp/Ensembl_GRCh37.71/gtf/Homo_sapiens.GRCh37.71_reduced.gff3"
miso.index <- "/home/Shared/data/annotation/Human_tmp/Ensembl_GRCh37.71/miso_0.4.9_index"


# create miso index
cmd <- paste0(miso.path, "index_gff.py --index ", gff3, " ", miso.index)
cat(cmd)
system(cmd)


## PAIRED-END samples


# Get isoform expression estimates
cmd.p <- paste0(miso.path, "run_events_analysis.py --compute-genes-psi ", miso.index, " /home/Shared/tmp/Simulation/simulation_hsapiens_V2/sam/sim_samp", metadata$SampleName1, "_s.bam --output-dir MISOCounts/miso_samp", metadata$SampleName1, " --read-len ", 75," -p 10 --paired-end 200 25 " ,"\n") 
cat(cmd.p)


for(i in 1:length(cmd.p)){
system(cmd.p[i])
}



### create count table 

# sort genes regarding gene ID
all.miso <- dir("MISOCounts/","miso$",recursive=TRUE)
all.genes <- gsub(".miso","",basename(all.miso))

# get the count for each event, gene and replicate
allgenesname <- sapply(all.miso,function(u) strsplit(u,"/")[[1]][1])
allgenesname2 <- sapply(all.miso,function(u) strsplit(u,"/")[[1]][3])
names(all.miso) <- paste(allgenesname,":",allgenesname2,sep="")

# for each .miso file, create a named vector with the event (names) and the count (numeric)
g <- lapply(all.miso, function(u) {
  rt <- read.table(paste0("MISOCounts/", u), comment.char="",nrow=1, sep="\t", stringsAsFactors=FALSE)$V8
  gsub("counts=\\(","",rt)
})

h <- lapply(g, function(u) strsplit(u,",\\(")[[1]])
j <- lapply(h, function(u) {
  z <- unlist(strsplit(u,"\\):"))
  n <- length(z)
  x <- as.numeric(z[seq(2,n,by=2)])
  names(x) <- z[seq(1,n,by=2)]
  x
})

# break into a list for each gene
js <- split(j, all.genes)
jsl <- sapply(js,length)

sum(jsl==nrow(metadata))
sum(!jsl==nrow(metadata))

# only keep genes with data in all samples (could be improved)
js <- js[jsl==nrow(metadata)]  

# function to go from list to table
getTableFromList <- function(u) {
  cn <- sapply(names(u), function(v) strsplit(v,":")[[1]][1],USE.NAMES=FALSE)
  rn <- unique(unlist(lapply(u,names)))
  x <- matrix(0,nrow=length(rn),ncol=length(cn),dimnames=list(rn,cn))
  for(i in 1:length(u))
    x[names(u[[i]]),cn[i]] <- u[[i]]
  x
}

ks <- lapply(js, getTableFromList)

n <- sapply(ks,nrow)
rn <- rep(names(ks),n)

# create 1 big table of everything
MISO_count <- do.call("rbind",ks)
rownames(MISO_count) <- paste0(rn,":",rownames(MISO_count))

write.table(data.frame(isoform=rownames(MISO_count), MISO_count), "MISOCounts/MISO_counts_all_samps.txt", quote=FALSE, sep="\t", row.names=FALSE)


# save counts per sample 

for(i in 1:nrow(metadata)){
  
  write.table(data.frame(rownames(MISO_count), MISO_count[,metadata$SampleName1[i]]), paste0("MISOCounts/MISO_counts_", metadata$SampleName1[i], ".txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  
}










