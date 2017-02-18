#Function: ghap.maf
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute marker minor allele frequencies

ghap.maf<-function(phase,only.active.samples=TRUE,only.active.markers=TRUE,ncores=1){
  
  #Check if phase is a GHap object
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  #Check if inactive markers and samples should be reactived
  if(only.active.markers == FALSE){
    phase$marker.in <- rep(TRUE,times=phase$nmarkers)
    phase$nmarkers.in<-length(which(phase$marker.in))
  }
  if(only.active.samples == FALSE){
    phase$id.in <- rep(TRUE,times=2*phase$nsamples)
    phase$nsamples.in<-length(which(phase$id.in))/2
  }
  
  maf.FUN <- function(i){
    p <- sum(phase$phase[i,phase$id.in])/(2*phase$nsamples.in)
    p <- min(c(p,1-p))
    names(p) <- phase$marker[i]
    return(p)
  }
  
  ncores <- min(c(detectCores(),ncores))
  if(Sys.info()["sysname"] == "Windows"){
    cat("\nParallelization not supported yet under Windows (using a single core).\n")
    maf <- lapply(X = which(phase$marker.in), FUN = maf.FUN)
  }else{
    maf <- mclapply(FUN = maf.FUN, X = which(phase$marker.in), mc.cores = ncores)
  }
  maf <- unlist(maf)
  
}