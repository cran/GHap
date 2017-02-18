#Function: ghap.blockgen
#License: GPLv3 or later
#Modification date: 18 Feb 2017
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Generate blocks based on sliding windows

ghap.blockgen<-function(
  phase,
  windowsize=10,
  slide=5,
  unit="marker",
  nsnp=2
){
  
  #Check if phase is a GHap.phase object
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  if(unit %in% c("marker","kbp") == FALSE){
    stop("Unit must be specified as 'marker' or 'kbp'")
  }
  
  #Initialize vectors
  BP1 <- rep(NA,times=phase$nmarkers.in)
  BP2 <- rep(NA,times=phase$nmarkers.in)
  SIZE <- rep(NA,times=phase$nmarkers.in)
  NSNP <- rep(NA,times=phase$nmarkers.in)
  
  if(unit == "kbp"){
    
    #Transform windowsize to bp
    windowsize <- windowsize*1e+3
    slide <- slide*1e+3
    
    #Get bp from active markers
    markersin <- which(phase$marker.in)
    bp <- phase$bp[markersin]
    
    #Generate indices
    id1 <- seq(1,max(bp),by=slide)
    id2 <- id1 + windowsize
    id1 <- id1[id2 <= max(bp)]
    id2 <- id2[id2 <= max(bp)]
    for(i in 1:length(id1)){
      slice <- which(bp >= id1[i] & bp <= id2[i])
      BP1[i] <- id1[i]
      BP2[i] <- id2[i]
      NSNP[i] <- length(slice)
    }
    
  }else if(unit == "marker"){
    
    id1<-seq(1,phase$nmarkers.in,by=slide);
    id2<-id1+(windowsize-1);
    id1<-id1[id2<=phase$nmarkers.in]
    id2<-id2[id2<=phase$nmarkers.in]
    markersin <- which(phase$marker.in)
    for(i in 1:length(id1)){
      slice <- markersin[id1[i]:id2[i]]
      BP1[i] <- phase$bp[slice[1]]
      BP2[i] <- phase$bp[slice[length(slice)]]
      NSNP[i] <- length(slice)
    }
    
  }
  
  results <- data.frame(BP1,BP2,NSNP,stringsAsFactors = FALSE)
  results <- unique(results)
  results$SIZE <- results$BP2 - results$BP1
  results$SIZE[results$NSNP == 1] <- 1
  results <- results[order(results$BP1,results$BP2),]
  results <- results[results$NSNP >= nsnp,]
  if(nrow(results) == 0){
    stop("No blocks could be generated with the specified options. Try setting the nsnp argument to a smaller value.")
  }
  results <- na.exclude(results)
  results$BLOCK <- paste("CHR",phase$chr,"_B",1:nrow(results),sep="")
  results$CHR <- phase$chr
  results <- results[,c("BLOCK","CHR","BP1","BP2","SIZE","NSNP")]
  return(results);
  
}