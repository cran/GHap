\name{ghap.kinv}
\alias{ghap.kinv}

\title{
Inverse of kinship matrix
}
\description{
Inversion of the haplotype covariance matrix
}
\usage{
ghap.kinv(kinship,  method="nearPD", ncores=1, proven=NULL)
}

\arguments{
  \item{kinship}{
   A HapAllele-based kinship matrix, as supplied by \code{\link{ghap.kinship}}.
  }
  \item{method}{
   Inversion method: common inverse, inverse of the nearest positive definite matrix ("nearPD", default) and approximate inverse based on the Algorithm for Proven and Young animals ("APY").
  }
  \item{ncores}{
  Number of cores to be used for computations (default = 1). Only relevant for method "APY".
 }
  \item{proven}{
   Character vector with the names of proven subjects to be used as reference for method "APY".
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
# # Subset data - markers with maf > 0.05
# maf <- ghap.maf(phase, ncores = 2)
# markers <- phase$marker[maf > 0.05]
# phase <- ghap.subsetphase(phase, unique(phase$id), markers)
# 
# # Generate blocks of 5 markers sliding 5 markers at a time
# blocks.mkr <- ghap.blockgen(phase, windowsize = 5, slide = 5, unit = "marker")
#
# # Generate matrix of haplotype genotypes
# ghap.haplotyping(phase, blocks.mkr, batchsize = 100, ncores = 2, outfile = "human")
#
# # Load haplotype genotypes
# haplo <- ghap.loadhaplo("human.hapsamples", "human.hapalleles", "human.hapgenotypes")
# 
# # Exclude minor alleles and singletons
# hapstats <- ghap.hapstats(haplo, ncores = 2)
# haplo <- ghap.subsethaplo(haplo,ids=haplo$id,alleles = hapstats$TYPE \%in\% c("REGULAR","MAJOR"))
#
# # Compute Kinship matrix
# K <- ghap.kinship(haplo, batchsize = 100)
#
#
# ### RUN ###
# 
# # Common inverse
# Kinv <- ghap.kinv(K, method="common")
#
# # Inverse of nearest positive definite matrix
# Kinv <- ghap.kinv(K, method="nearPD")
#
# # APY inverse
# Kinv <- ghap.kinv(K, method="APY", proven=colnames(K)[1:500], ncores=2)

}