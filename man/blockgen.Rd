\name{ghap.blockgen}
\alias{ghap.blockgen}

\title{
 Haplotype block generator
}
\description{
 This function generates HapBlocks based on sliding windows. The window and the step size can be specified in markers or kb. For each window, block coordinates are generated.
}
\usage{
 ghap.blockgen(phase, windowsize = 10, slide = 5, unit = "marker")
}

\arguments{
  \item{phase}{
   A GHap.phase object
  }
  \item{windowsize}{
   A numeric value for the size of the window (default = 10).
  }
  \item{slide}{
   A numeric value for the step size (default = 5).
  }
  \item{unit}{
  A character value for the size unit used for the window and the step. It can be either "marker" or "kb" (default = "marker").
  }
}

\value{
 A data frame with columns:
 \item{BLOCK}{
  Block alias.
 }
 \item{CHR}{
  Chromosome name.
 }
 \item{BP1}{
  Block start position.
 }
 \item{BP2}{
  Block end position.
 }
 \item{SIZE}{
  Haplotype size.
 }
 \item{NSNP}{
  Number of marker.
 }
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy the example data in the current working directory
# ghap.makefile()
# 
# # Load data
# phase <- ghap.loadphase("human.samples", "human.markers", "human.phase")
# 
# # Subset data - randomly select 3000 markers with maf > 0.02
# maf <- ghap.maf(phase, ncores = 2)
# set.seed(1988)
# markers <- sample(phase$marker[maf > 0.02], 3000, replace = FALSE)
# phase <- ghap.subsetphase(phase, unique(phase$id), markers)
# rm(maf,markers)
# 
# 
# ### RUN ###
# 
# # Generate block coordinates based on windows of 10 markers, sliding 5 marker at a time
# blocks <- ghap.blockgen(phase, 10, 5, "marker")
# 
# # Generate block coordinates based on windows of 500 kb, sliding 100 kb at a time
# blocks <- ghap.blockgen(phase, 500, 100, "kb")

}

