#Function: ghap.blockstats
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Calculate block summary statistics

ghap.blockstats <- function(hapstats, ncores = 1){
  
  #Get unique blocks
  blocks <- unique(hapstats[,c("BLOCK","CHR","BP1","BP2")])
  
  #Set of internal functions
  my.fun <- function(block){
    freq <- hapstats$FREQ[hapstats$BLOCK==block & hapstats$TYPE != "ABSENT"]
    if(sum(freq) == 1){
      exp.het <- 1-(sum((freq)^2))
    }else{
      exp.het <- NA
    }
    return(c(exp.het,length(freq)))
  }
    
  #Calculation of expected heterozygosity
  ncores <- min(c(detectCores(),ncores))
  if(Sys.info()["sysname"] == "Windows"){
    cat("\nParallelization not supported yet under Windows (using a single core).\n")
    temp <- lapply(FUN = my.fun, X = blocks$BLOCK)
  }else{
    temp <- mclapply(FUN = my.fun, X = blocks$BLOCK, mc.cores = ncores)
  }
  temp <- unlist(temp)
  blocks$EXP.H <- temp[1:length(temp) %% 2 == 1]
  blocks$N.ALLELES <- temp[1:length(temp) %% 2 == 0]
  
  #Return output
  return(blocks)
  
}

