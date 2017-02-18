#Function: ghap.assoc
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: fit ordinary least squares for each haplotype allele

ghap.assoc<-function(
  response,
  haplo,
  weights=NULL,
  #type="HapAllele",
  gc=TRUE,
  only.active.alleles=TRUE,
  ncores=1
){
  
  
  #General data check
  if (class(haplo) != "GHap.haplo") {
    stop("Argument haplo must be a GHap.haplo object.")
  }
  if (only.active.alleles == FALSE) {
    haplo$allele.in <- rep(TRUE, times = haplo$nalleles)
    haplo$nalleles.in <- length(which(haplo$allele.in))
  }
  if (length(which(names(response) %in% haplo$id)) != length(response)) {
    stop("All ids in the response must be present in the GHap.haplo object.")
  }
  ids <- rep(NA, times = length(response))
  for (i in 1:length(ids)) {
    ids[i] <- which(haplo$id == names(response)[i])
  }
  
  #Check if response contains missing values
  if(length(na.omit(response)) != length(response)){
    stop("Missing values are not accepted.")
  }
  
  #Include weights
  if(is.null(weights)==FALSE){
    w <- sqrt(weights/mean(weights))
    response <- w*response
  }
  
  #Center response
  response <- response - mean(response)
  
#  if(type == "HapAllele"){
    
    #ols iterate function
    ols.FUN <- function(j){
      x <- haplo$genotypes[j,unique(ids)]
      frq <- sum(x)/(2*length(x))
      x <- haplo$genotypes[j,ids]
      x <- (x-mean(x))/sd(x)
      if(is.null(weights)==F){
        x <- w*x
      }
      xpxi <- 1/sum(x^2)
      xpy <- sum(x*response)
      b <- xpxi*xpy
      se <- sqrt(var(response - x*b)*xpxi)
      return(c(b,se,frq))
    }
    
    #Compute haplotype regression statistics
    ncores <- min(c(detectCores(),ncores))
    if(Sys.info()["sysname"] == "Windows"){
      cat("\nParallelization not supported yet under Windows (using a single core).\n")
      a <- lapply(FUN = ols.FUN, X = which(haplo$allele.in))
    }else{
      a <- mclapply(FUN = ols.FUN, X = which(haplo$allele.in), mc.cores = ncores)
    }
    a <- data.frame(matrix(unlist(a), nrow=haplo$nalleles.in, byrow=TRUE))
    hapreg <- NULL
    hapreg$BLOCK <- haplo$block[haplo$allele.in]
    hapreg$CHR <- haplo$chr[haplo$allele.in]
    hapreg$BP1 <- haplo$bp1[haplo$allele.in]
    hapreg$BP2 <- haplo$bp2[haplo$allele.in]
    hapreg$ALLELE <- haplo$allele[haplo$allele.in]
    hapreg$BETA <- a[,1]
    hapreg$SE <- a[,2]
    hapreg$FREQ <- a[,3]
    hapreg$CHISQ.OBS <- (hapreg$BETA^2)/(hapreg$SE^2)
    hapreg$CHISQ.EXP <- qchisq(p = rank(hapreg$CHISQ.OBS)/(haplo$nalleles.in+1), df=1)
    if(gc == TRUE){
      dev <- sd(hapreg$CHISQ.OBS)
      samp <- which(hapreg$CHISQ.OBS < 3*dev)
      lambda <- lm(hapreg$CHISQ.OBS[samp] ~ hapreg$CHISQ.EXP[samp])
      lambda <- lambda$coefficients[2]
      hapreg$CHISQ.OBS <- hapreg$CHISQ.OBS/lambda
    }
    hapreg$logP <- -1*pchisq(q = hapreg$CHISQ.OBS, df = 1, lower.tail=FALSE, log.p = TRUE)/log(10)
    hapreg <- data.frame(hapreg,stringsAsFactors = FALSE)
  # } else if(type == "HapBlock") {
  #   
  #   #ols iterate function for HapBlock
  #   ols.FUN <- function(j){
  #     hapalleles <- which(haplo$block == j)
  #     nalleles <- length(hapalleles)
  #     if(nalleles == 0){
  #       nu1 <- NA
  #       nu2 <- NA
  #       ftest <- NA
  #     }else if(nalleles == 1){
  #       X <- as.matrix(haplo$genotypes[hapalleles,unique(ids)])
  #       freq <- sum(X)/(2*haplo$nsamples.in)
  #       X <- as.matrix(haplo$genotypes[hapalleles,ids])
  #       if(freq == 1){
  #         nu1 <- NA
  #         nu2 <- NA
  #         ftest <- NA
  #       }else{
  #         X <- X - mean(X)
  #         if(is.null(weights)==F){
  #           X <- w*X
  #         }
  #         b <- sum(X*response)/(sum(X^2))
  #         Xb <- sum(X*b)
  #         RSS <- sum((response - Xb)^2)
  #         ESS <- sum((Xb - mean(response))^2)
  #         nu1 <- 1
  #         nu2 <- length(response)-1
  #         ftest <- (ESS/nu1)/(RSS/nu2)
  #       }
  #     }else if(nalleles == 2){
  #       X <- as.matrix(haplo$genotypes[hapalleles,unique(ids)])
  #       freq <- rowSums(X)/(2*haplo$nsamples.in)
  #       X <- as.matrix(haplo$genotypes[hapalleles,ids])
  #       if(sum(freq) == 1){
  #         haporder <- order(rank(freq))
  #         freq <- freq[haporder]
  #         freq <- freq[-1]
  #         X <- X[haporder,]
  #         X <- X[-1,]
  #         X <- X - mean(X)
  #         if(is.null(weights)==F){
  #           X <- w*X
  #         }
  #         b <- sum(X*response)/(sum(X^2))
  #         Xb <- sum(X*b)
  #         RSS <- sum((response - Xb)^2)
  #         ESS <- sum((Xb - mean(response))^2)
  #         nu1 <- 1
  #         nu2 <- length(response)-1
  #         ftest <- (ESS/nu1)/(RSS/nu2)
  #       }else{
  #         X <- apply(X = X, MARGIN = 1, FUN = function(x){x-mean(x)})
  #         if(is.null(weights)==F){
  #           X <- w*X
  #         }
  #         SQR <- qr(X)
  #         b <- qr.coef(qr=SQR,y=response)
  #         Xb <- X%*%b
  #         RSS <- sum((response - Xb)^2)
  #         ESS <- sum((Xb - mean(response))^2)
  #         nu1 <- nalleles-1
  #         nu2 <- length(response)-nalleles-1
  #         ftest <- (ESS/nu1)/(RSS/nu2)
  #       }
  #     }else if(nalleles > 2){
  #       X <- as.matrix(haplo$genotypes[hapalleles,unique(ids)])
  #       freq <- rowSums(X)/(2*haplo$nsamples.in)
  #       X <- as.matrix(haplo$genotypes[hapalleles,ids])
  #       if(sum(freq) == 1){
  #         haporder <- order(rank(freq))
  #         freq <- freq[haporder]
  #         freq <- freq[-1]
  #         X <- X[haporder,]
  #         X <- X[-1,]
  #         X <- apply(X = X, MARGIN = 1, FUN = function(x){x-mean(x)})
  #         if(is.null(weights)==F){
  #           X <- w*X
  #         }
  #         SQR <- qr(X)
  #         b <- qr.coef(qr=SQR,y=response)
  #         Xb <- X%*%b
  #         RSS <- sum((response - Xb)^2)
  #         ESS <- sum((Xb - mean(response))^2)
  #         nu1 <- nalleles-1
  #         nu2 <- length(response)-nalleles-1
  #         ftest <- (ESS/nu1)/(RSS/nu2)
  #       }else{
  #         X <- apply(X = X, MARGIN = 1, FUN = function(x){x-mean(x)})
  #         if(is.null(weights)==F){
  #           X <- w*X
  #         }
  #         SQR <- qr(X)
  #         b <- qr.coef(qr=SQR,y=response)
  #         Xb <- X%*%b
  #         RSS <- sum((response - Xb)^2)
  #         ESS <- sum((Xb - mean(response))^2)
  #         nu1 <- nalleles-1
  #         nu2 <- length(response)-nalleles-1
  #         ftest <- (ESS/nu1)/(RSS/nu2)
  #         
  #       }
  #     }
  #     return(c(nalleles,nu1,nu2,ftest))
  #   }
  #   
  #   #Get unique blocks
  #   blocks <- unique(data.frame(haplo$block,haplo$chr,haplo$bp1,haplo$bp2,stringsAsFactors = FALSE))
  #   colnames(blocks) <- c("BLOCK","CHR","BP1","BP2")
  #   
  #   #Compute haplotype regression statistics
  #   a <- mclapply(FUN=ols.FUN,X=blocks$BLOCK,mc.cores = ncores)
  #   a <- data.frame(matrix(unlist(a), nrow=nrow(blocks), byrow=TRUE))
  #   
  #   hapreg <- NULL
  #   hapreg$BLOCK <- blocks$BLOCK
  #   hapreg$CHR <- blocks$CHR
  #   hapreg$BP1 <- blocks$BP1
  #   hapreg$BP2 <- blocks$BP1
  #   hapreg$N.ALLELES <- as.integer(a[,1])
  #   hapreg$F.TEST <- a[,4]
  #   hapreg$logP <- -1*pf(q = hapreg$F.TEST, df1 = a[,2], df2 = a[,3], lower.tail=FALSE, log.p = TRUE)/log(10)
  #   hapreg <- data.frame(hapreg,stringsAsFactors = FALSE)
  # } else {
  #   stop("Argument type must be 'HapAllele' or 'HapBlock'")
  # }
  
  
  #Return object
  return(hapreg)
  
}