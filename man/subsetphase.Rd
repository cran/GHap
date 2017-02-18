\name{ghap.subsetphase}
\alias{ghap.subsetphase}

\title{
 Subset GHap.phase object
}

\description{
 This function takes a list of markers and individuals and subsets a GHap.phase object.
}
\usage{
ghap.subsetphase(phase, ids, markers, verbose = TRUE)
}

\arguments{
  \item{phase}{
   A GHap.phase object.
  }
  \item{ids}{
   Character vector of individual names to keep.
  }
  \item{markers}{
   Character vector of marker names to keep.
  }
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
 The returned GHap.phase object (as described in the documentation for the \code{\link{ghap.loadphase}} function) is the same as the one used in the phase argument. However, individuals and markers not included in the provided vectors are set to FALSE (i.e., inactivated) in the phase$samples.in and phase$marker.in vectors, respectively. This procedure avoids expensive subsetting operations by simply flagging which markers and individuals should be used in downstream analyses.
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
# 
# ### RUN ###
# 
# # Subset data - markers with maf > 0.05
# maf <- ghap.maf(phase, ncores = 2)
# markers <- phase$marker[maf > 0.05]
# phase <- ghap.subsetphase(phase, unique(phase$id), markers)

}
