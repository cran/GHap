#Function: ghap.gblup
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: calculate GBLUP solution for each haplotype allele

ghap.blup<-function(
  gebvs, 
  haplo,
  invcov,
  gebvsweights = NULL,
  haploweights = NULL,
  nperm = 1,
  only.active.alleles = TRUE, 
  ncores = 1
){
  
  #General data check
  if (class(haplo) != "GHap.haplo") {
    stop("Argument haplo must be a GHap.haplo object.")
  }
  if (only.active.alleles == FALSE) {
    haplo$allele.in <- rep(TRUE, times = haplo$nalleles)
    haplo$nalleles.in <- length(which(haplo$allele.in))
  }
  activealleles <- which(haplo$allele.in)
  if (is.null(haploweights) == TRUE) {
    haploweights <- rep(1,times=haplo$nalleles.in)
  }
  if (length(haploweights) != haplo$nalleles.in) {
    stop("Vector of haplotype weights must have the same length as the number of haplotype alleles.")
  }
  names(haploweights) <- activealleles
  ids <- rep(NA, times = length(gebvs))
  for (i in 1:length(ids)) {
    ids[i] <- which(haplo$id == names(gebvs)[i])
  }
  if (length(which(names(gebvs) %in% haplo$id)) != length(gebvs)) {
    stop("All ids in the vector of GEBVs must be present in the GHap.haplo object.")
  }
  invcov <- invcov[names(gebvs),names(gebvs)]
  if(identical(colnames(invcov),names(gebvs)) == FALSE){
    stop("Names declared in random effects and covariance matrix do not match!")
  }
  if (is.null(gebvsweights) == TRUE) {
    gebvsweights <- rep(1,times=length(gebvs))
  }
  if (length(gebvsweights) != length(gebvsweights)) {
    stop("Vector of GEBV weights must have the same length as the number of observations.")
  }
  names(gebvsweights) <- names(gebvs)
  k <- invcov%*%Diagonal(x=gebvsweights/mean(gebvsweights))%*%gebvs
  
  #Compute variance and frequencies
  ncores <- min(c(detectCores(),ncores))
  varfun <- function(j) return(haploweights[as.character(j)]*var(haplo$genotypes[j,ids]))
  freqfun <- function(j) return(sum(haplo$genotypes[j,ids])/(2*haplo$nsamples.in))
  if(Sys.info()["sysname"] == "Windows"){
    cat("\nParallelization not supported yet under Windows (using a single core).\n")
    sumvar <- sum(unlist(lapply(X = activealleles, FUN = varfun)))
    freq <- unlist(lapply(X = activealleles, FUN = freqfun))
  }else{
    sumvar <- sum(unlist(mclapply(X=activealleles,FUN = varfun,  mc.cores = ncores)))
    freq <- unlist(mclapply(X=activealleles,FUN = freqfun,  mc.cores = ncores))
  }
  
  #Main BLUP function
  gblup.FUN <- function(j) {
    x <- haplo$genotypes[j, ids]
    cent <- mean(x)
    x <- x - cent
    b <- sum(haploweights[as.character(j)]*x*k)
    b <- b/sumvar
    varxb <- var(x*b)
    return(c(b,varxb,cent))
  }
  
  #Compute effects
  if(Sys.info()["sysname"] == "Windows"){
    cat("\nParallelization not supported yet under Windows (using a single core).\n")
    a <- lapply(FUN = gblup.FUN, X = activealleles)
  }else{
    a <- mclapply(FUN = gblup.FUN, X = activealleles, mc.cores = ncores)
  }
  a <- data.frame(matrix(unlist(a), nrow=haplo$nalleles.in, byrow=TRUE))
  
  #Output data
  hapreg <- NULL
  hapreg$BLOCK <- haplo$block[haplo$allele.in]
  hapreg$CHR <- haplo$chr[haplo$allele.in]
  hapreg$BP1 <- haplo$bp1[haplo$allele.in]
  hapreg$BP2 <- haplo$bp2[haplo$allele.in]
  hapreg$ALLELE <- haplo$allele[haplo$allele.in]
  hapreg$SCORE <- a[,1]
  hapreg$FREQ <- freq
  hapreg$VAR <- a[,2]
  hapreg$pVAR <- hapreg$VAR/sum(hapreg$VAR)
  hapreg$CENTER <- a[,3]
  hapreg$SCALE <- 1
  
  #Permutation test (optional)
  if(nperm > 1){
    cat("A permutation procedure with",nperm,"randomizations will be performed.\n")
    #Re-define BLUP function to randomize effects
    #Main BLUP function
    gblup.FUN <- function(j) {
      x <- haplo$genotypes[j, ids]
      x <- x - mean(x)
      b <- sum(haploweights[as.character(j)]*x*k)
      b <- b/sumvar
      return(b)
    }
    hapreg$P <- 0
    #Permutation iteration
    for(i in 1:nperm){
      cat("Permutation number:",i,"\r")
      k <- sample(k,size=length(k),replace=FALSE)
      if(Sys.info()["sysname"] == "Windows"){
        cat("\nParallelization not supported yet under Windows (using a single core).\n")
        a <- lapply(FUN = gblup.FUN, X = activealleles)
      }else{
        a <- mclapply(FUN = gblup.FUN, X = activealleles, mc.cores = ncores)
      }
      a <- unlist(a)
      a <- as.numeric(abs(max(a)) > abs(hapreg$SCORE))
      hapreg$P <- hapreg$P + a
    }
    hapreg$P <- hapreg$P/nperm
    hapreg$P[hapreg$P == 0] <- 1/nperm
  }
  
  #Return results
  hapreg <- data.frame(hapreg, stringsAsFactors = FALSE)
  return(hapreg)
}