\name{ghap.mergehaplo}
\alias{ghap.mergehaplo}

\title{
  Merging GHap.haplo objects
}
\description{
  This function merges two GHap.haplo objects.
}
\usage{
  ghap.mergehaplo(haplo.1, haplo.2, type, only.active.markers = TRUE,
                  only.active.samples = TRUE, verbose = TRUE)
}

\arguments{
  \item{haplo.1}{
    First GHap.haplo object.
  }
  \item{haplo.2}{
    Second GHap.haplo object.
  }
  \item{type}{
    A character value specifying the merging task to be performed. If type="HapAlleles", the two objects are assumed to have the same individuals but different HapAlleles. If type="individual", the two objects are assumed to have different individuals scored for the same HapAlleles.
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be included in the output (default = TRUE).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
  The returned GHap.haplo object is a list with components specified by \code{\link{ghap.loadhaplo}}. Please notice that no duplicated individuals or HapAlleles (concatenation of \emph{haplo$block} and \emph{haplo$allele}) are allowed.
}

\author{
  Marco Milanesi <marco.milanesi.mm@gmail.com>
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
# 
# ### RUN ###
# 
# # Select CEU and ASW individuals
# CEU.ids <- haplo$id[which(haplo$pop=="CEU")]
# ASW.ids <- haplo$id[which(haplo$pop=="ASW")]
# 
# # Randomly select a set of HapAlleles
# set.seed(1988)
# random.allele <- sample(x=c(TRUE, FALSE), size=haplo$nalleles, replace = TRUE)
# 
# # Subset data 
# haplo.CEU <- ghap.subsethaplo(haplo,ASW.ids,random.allele)
# haplo.ASW <- ghap.subsethaplo(haplo,CEU.ids,random.allele)
# 
# # Merge haplo.CEU and haplo.ASW
# haplo.merge <- ghap.mergehaplo(haplo.1 = haplo.CEU, haplo.2 = haplo.ASW, type = "individual")

  
}

