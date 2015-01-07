############################################################################
# function to simulate data from DM
# Created 19 Nov 2014
############################################################################



simulate_from_DM <- function(sample.size = 5, s = "1" , pi.org = c(1/3, 1/3, 1/3), g0.org = 100, nr.genes = 100, nM = 150, tot="nbinom", nD=3, out.dir, mc.cores=10, save=TRUE){

  
  dir.create(out.dir, recursive = T, showWarnings = FALSE)
  name1 <- paste0("SIM", s , "_")
  # simParams <- as.list(match.call())
  # simParams <- c(mget(names(formals()), envir = as.environment(-1)), sapply(as.list(substitute({ ... })[-1]), deparse))
  simParams <- mget(names(formals()), envir = as.environment(-1))
  sink(file=paste0(out.dir, "/", name1, "simParams",".txt"), split=TRUE)
  print(simParams)
  sink(file=NULL)
  
  if(length(g0.org == 1))
    g0.org <- rep(g0.org, nr.genes)

  
  sim <- mclapply(1:nr.genes, function(i){
    # i=1
    
    g.dir.org <- pi.org * g0.org[i] # gamma
    
    # simulate dirichlets
    g.dir <- rdirichlet( sample.size, g.dir.org)
    
    # simulate total counts
    if(tot=="nbinom")
      t <- rnbinom(sample.size, mu=nM, size=nD)
    if(tot=="norm")
      t <- round(rnorm(sample.size, mean=nM, sd=nD))  
    if(tot=="uni")
      t <- rep(nM, sample.size)
    
    # simulate multinomial
    d <- sapply(1:sample.size, function(u) rmultinom(1, prob=g.dir[u,], size=t[u]))    
    genes <- paste0("g", rep(i, length(g.dir.org)), ":e", 1:length(g.dir.org))   
    rownames(d) <- genes
    
    return(d)
    
  }, mc.cores=mc.cores)
  
  sim <- do.call(rbind, sim)
  ids <- strsplit2(rownames(sim), ":", fixed = TRUE)
  dge <- DGEList( counts=sim, group = rep(1, sample.size), genes=data.frame(gene_id=ids[,1], ete_id = rownames(sim)) )
  
  dge$counts <- split(data.frame(dge$counts), factor(dge$genes$gene_id, levels=unique(dge$genes$gene_id)))
  dge$counts <- lapply(dge$counts, as.matrix)  ## !!! have to conver into matrix, othewise ERROR
  
  
  g.dir.org <- pi.org * g0.org[1] ### Normally do not need this param to return but I keep it by now

  if(save){
    dir.create(out.dir, recursive = T, showWarnings = F)
    save(dge, name1, g.dir.org, file=paste0(out.dir, "/", name1, "dge",".RData"))   
    save(simParams, file = paste0(out.dir, "/", name1, "simParams",".RData"))
  }
  
  
  return(list(dge=dge, name1=name1, g.dir.org=g.dir.org))
  
  
}


