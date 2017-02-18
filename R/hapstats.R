#Function: ghap.hapstats
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Summary statistics for haplotype alleles

ghap.hapstats<-function(
  haplo,
  alpha=c(1,1),
  only.active.samples=TRUE,
  only.active.alleles=TRUE,
  ncores=1
){
  
  
  #Check if haplo is a GHap.haplo object
  if(class(haplo) != "GHap.haplo"){
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Check if inactive alleles and samples should be reactived
  if(only.active.alleles == FALSE){
    haplo$allele.in <- rep(TRUE,times=haplo$nalleles)
    haplo$nalleles.in<-length(which(haplo$allele.in))
  }
  if(only.active.samples == FALSE){
    haplo$id.in <- rep(TRUE,times=haplo$nsamples)
    haplo$nsamples.in<-length(which(haplo$id.in))
  }
  
  #Hapstats iterate function
  hapstats.FUN <- function(j){
    hap.geno <- haplo$genotypes[j,haplo$id.in]
    N<-sum(hap.geno)
    FREQ<-N/(2*haplo$nsamples.in)
    O.HOM <- length(which(hap.geno == 2))
    O.HET <- length(which(hap.geno == 1))
    return(c(N,FREQ,O.HOM,O.HET))
  }
  
  #Compute haplotype allele statistics
  ncores <- min(c(detectCores(),ncores))
  if(Sys.info()["sysname"] == "Windows"){
    cat("\nParallelization not supported yet under Windows (using a single core).\n")
    a <- lapply(X = which(haplo$allele.in), FUN = hapstats.FUN)
  }else{
    a <- mclapply(FUN=hapstats.FUN,X=which(haplo$allele.in),mc.cores = ncores)
  }
  a <- data.frame(matrix(unlist(a), nrow=haplo$nalleles.in, byrow=TRUE))
  hapstats <- NULL
  hapstats$BLOCK <- haplo$block[haplo$allele.in]
  hapstats$CHR <- haplo$chr[haplo$allele.in]
  hapstats$BP1 <- haplo$bp1[haplo$allele.in]
  hapstats$BP2 <- haplo$bp2[haplo$allele.in]
  hapstats$ALLELE<- haplo$allele[haplo$allele.in]
  hapstats$N <- a[,1]
  hapstats$FREQ <- a[,2]
  hapstats$O.HOM <- a[,3]
  hapstats$O.HET <- a[,4]
  hapstats$E.HOM <- (hapstats$FREQ^2)*haplo$nsamples.in
  hapstats$RATIO <- (hapstats$E.HOM+alpha[1])/(hapstats$O.HOM+alpha[2])
  hapstats$BIN.logP <- -1*pbinom(q = hapstats$O.HOM,size = haplo$nsamples.in,prob = hapstats$FREQ^2,lower.tail=TRUE,log.p = TRUE)/log(10)
  hapstats$POI.logP <- -1*ppois(q = hapstats$O.HOM,lambda = hapstats$E.HOM,lower.tail=TRUE,log.p = TRUE)/log(10)
  hapstats <- data.frame(hapstats,stringsAsFactors = FALSE)
  hapstats$TYPE <- NA
  for(i in unique(hapstats$BLOCK)){
    slice <- which(hapstats$BLOCK == i)
    freq <- hapstats$FREQ[slice]
    sumfreq <- sum(freq)
    nalleles <- length(slice)
    type <- rep("REGULAR",times=nalleles)
    type[freq == 0] <- "ABSENT"
    minfreq <- min(freq[freq != 0])
    maxfreq <- max(freq)
    if(sumfreq == 1 & nalleles > 2){
      type[which(freq == minfreq)[1]] <- "MINOR"
      type[which(freq == maxfreq)[1]] <- "MAJOR"
    }else if(sumfreq == 1 & nalleles == 2){
      if(freq[1] == freq[2]){
        type[1] <- "MINOR"
        type[2] <- "MAJOR"
      }else if(freq[1] != freq[2] & pmin(freq[1],freq[2]) != 0){
        type[which(freq == minfreq)] <- "MINOR"
        type[which(freq == maxfreq)] <- "MAJOR"
      }else if(maxfreq == 1){
        type[which(freq == 1)] <- "SINGLETON"
      }
    }else if(sumfreq == 1 & nalleles == 1){
      type <- "SINGLETON"
    # }else if(sumfreq != 1 & nalleles == 1){
    #   type <- "REGULAR"
    # }else{
    #   type[which(freq == maxfreq)[1]] <- "MAJOR"
    # }
    }else{
      type <- "REGULAR"
    }
    hapstats$TYPE[slice] <- type
  }

  #Return object
  return(hapstats)
  
}