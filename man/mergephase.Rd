\name{ghap.mergephase}
\alias{ghap.mergephase}

\title{
  Merging GHap.phase objects
}
\description{
  This function merges two GHap.phase objects.  
}
\usage{
  ghap.mergephase(phase.1, phase.2, only.active.markers = TRUE,
                  only.active.samples = TRUE, verbose = TRUE)
}

\arguments{
  \item{phase.1}{
    First GHap.phase object.
  }
  \item{phase.2}{
    Second GHap.phase object.
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
  The returned GHap.phase object is a list with components specified by \code{\link{ghap.loadphase}}. Merging is allowed only if the two original GHap.phase objects have the same markers in the same order. Duplicated individuals are not tolerated. 
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
# 
# ### RUN ###
# 
# # Randomly select 3000 markers with maf > 0.02
# maf <- ghap.maf(phase, ncores = 2)
# set.seed(1988)
# markers <- sample(phase$marker[maf > 0.02], 3000, replace = FALSE)
# 
# # Select ASW and CEU individuals
# ASW.ids <- unique(phase$id[which(phase$pop=="ASW")])
# CEU.ids <- unique(phase$id[which(phase$pop=="CEU")])
# 
# # Subset data 
# phase.ASW <- ghap.subsetphase(phase = phase, ids = ASW.ids, markers = markers)
# phase.CEU <- ghap.subsetphase(phase = phase, ids = CEU.ids, markers = markers)
# 
# # Merge phase.ASW and phase.CEU
# phase.merge <- ghap.mergephase(phase.1 = phase.ASW, phase.2 = phase.CEU)


}

