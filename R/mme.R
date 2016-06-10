#Function: ghap.mme
#License: GPLv3 or later
#Modification date: 6 Apr 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Mixed model equations solver

ghap.mme<-function(
  fixed,             #Formula for the fixed part of the model
  random,            #Column name with the random part of the model
  weights=NULL,      #Weights for records
  env.eff=FALSE,     #Should permanent environmental effects be included
  data,              #Data frame containing model data
  K,                 #Covariance matrix of random effects
  varcomp,           #Vector with variance components
  verbose=TRUE
){
  
  #Log message
  if(verbose==TRUE){
    cat("\nAssuming continuous data and known variance components.\n")
    cat("[This function does not support categorical data]\n")
    cat("\nAssembling design matrices... ")
  }
  
  #Build design matrices
  colnames(data)[colnames(data) == random] <- "MyRanColName"
  random <- "MyRanColName"
  mf <- model.frame(fixed, data = data)
  y <- model.response(mf)
  X <- model.matrix(mf,mf)
  random.fm <- as.formula(paste("~ 0 +",random))
  Z <- model.matrix(random.fm, data = data)
  colnames(Z) <- gsub(pattern = random, replacement = "", x = colnames(Z))
  Z <- Z[,colnames(K)]
  if(is.null(weights) == FALSE){
    if(length(weights) != length(y)){
      stop("The vector of residual weights must have the same length as the number of records.")
    }
  }else{
    weights <- rep(1,times=length(y))
  }
  w <- sqrt(1/weights)
  y <- w*y
  X <- w*X
  Z <- w*Z
  
  #Check variance components
  if(length(varcomp) == 2 & env.eff == FALSE){
    vare <- varcomp[1]
    varu <- varcomp[2]
  }else if(length(varcomp) == 3 & env.eff == TRUE){
    vare <- varcomp[1]
    varu <- varcomp[2]
    varp <- varcomp[3]
  }else{
    stop("Argument varcomp requires 2 elements for env.eff=FALSE and 3 elements for env.eff=TRUE")
  }
  
  #Check data
  if(identical(colnames(K),colnames(Z)) == FALSE){
    stop("Names declared in random effects and covariance matrix do not match!")
  }
  
  #Log message
  if(verbose==TRUE){
    cat("Done.\n")
    cat(length(y),"records will be fitted to an intercept,",ncol(X)-1,"fixed effects and",ncol(Z),"random effects.\n")
    if(env.eff==TRUE){
      cat("Additionally,",ncol(Z),"permanent environmental effects will be included.\n")
    }
    cat("\nComputing fixed and random effects via Henderson's equations... ")
  }
  
  #Inverse relationship matrix
  Kinv <- svd(K)
  Kinv <- Kinv$v%*%diag(1/Kinv$d)%*%t(Kinv$u)

  #Build LHS
  XpX <- crossprod(X)
  XpZ <- crossprod(X,Z)
  ZpX <- crossprod(Z,X)
  ZpZ <- crossprod(Z,Z)
  LHS <- cbind(XpX,XpZ)
  LHS <- rbind(LHS,cbind(ZpX,ZpZ + Kinv*(vare/varu)))
  if(env.eff == T){
    LHS <- cbind(LHS,rbind(XpZ,ZpZ))
    LHS <- rbind(LHS,cbind(ZpX,ZpZ,ZpZ + diag(ncol(ZpZ))*(vare/varp)))
  }
  
  #Build RHS
  Xpy <- crossprod(X,y)
  Zpy <- crossprod(Z,y)
  RHS <- rbind(Xpy,Zpy)
  if(env.eff == TRUE){
    RHS <- rbind(RHS,Zpy)
  }
  
  #Solve for effects
  LHSinv <- svd(LHS)
  LHSinv <- LHSinv$v%*%diag(1/LHSinv$d)%*%t(LHSinv$u)
  eff <- as.vector(LHSinv%*%RHS)
  
  #Log message
  if(verbose==TRUE){
    cat("Done.\n")
    cat("Assembling results... ")
  }
  
  #Prepare output
  results <- NULL
  results$b <- eff[1:ncol(X)]; names(results$b) <- colnames(X)
  results$u <- eff[(ncol(X)+1):(ncol(X)+ncol(Z))]; names(results$u) <- colnames(Z)
  if(env.eff == TRUE){
    results$p <- eff[(ncol(X)+ncol(Z)+1):length(eff)]; names(results$p) <- colnames(Z)
  }
  results$varu <- varu
  if(env.eff==TRUE){
    results$varp <- varp
  }
  results$vare <- vare
  if(env.eff==TRUE){
    results$h2 <- results$varu/(results$varu+results$vare+results$varp)
    results$H2 <- (results$varu+results$varp)/(results$varu+results$vare+results$varp)
  }else{
    results$h2 <- results$varu/(results$varu+results$vare)
  }
  results$k <- as.vector(Kinv%*%results$u); names(results$k) <- colnames(K)
  results$y <- (1/w)*y; names(results$y) <- data[,random]
  results$weights <- weights; names(results$weights) <- data[,random]
  if(env.eff==TRUE){
    results$residuals <- (1/w)*as.vector(results$y - X%*%results$b - Z%*%results$u - Z%*%results$p); names(results$residuals) <- data[,random]
  }else{
    results$residuals <- (1/w)*as.vector(results$y - X%*%results$b - Z%*%results$u); names(results$residuals) <- data[,random]
  }
  results$pdev <- ((-length(results$y)/2)*log(2*pi*results$vare)) + ((sum(results$residuals^2)/(2*results$vare)))
  class(results) <- "GHap.blmm"
  
  #Log message
  if(verbose==TRUE){
    cat("Done.\n\n")
  }
  
  return(results)
}

