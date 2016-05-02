\name{ghap.maf}
\alias{ghap.maf}

\title{
  Compute marker minor allele frequencies
}
\description{
  This function takes a GHap.phase object and computes the minor allele frequency for each marker.
}
\usage{
  ghap.maf(phase, only.active.samples = TRUE, only.active.markers = TRUE, ncores = 1)
}

\arguments{
  \item{phase}{
    A GHap.phase object.
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be used for calculations (default = TRUE).
  }
  \item{only.active.markers}{
    A logical value specifying whether only active markers should be included in the output (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
}

\value{
  The function outputs a numeric vector of the same length of active markers containing minor allele frequencies based on the active samples.
}

\author{
  Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
  
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
# ### RUN ###
# 
# # Subset data - randomly select 3000 markers with maf > 0.02
# maf <- ghap.maf(phase, ncores = 2)
# set.seed(1988)
# markers <- sample(phase$marker[maf > 0.02], 3000, replace = FALSE)
# phase <- ghap.subsetphase(phase, unique(phase$id), markers)
# rm(maf,markers)

}
