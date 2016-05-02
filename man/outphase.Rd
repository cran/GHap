\name{ghap.outphase}
\alias{ghap.outphase}

\title{
 Output GHap.phase object to a file
}
\description{
 This function saves a GHap.phase object to a text file using the same format specified by \code{\link{ghap.loadphase}}.
}
\usage{
 ghap.outphase(phase, outfile, only.active.markers = TRUE,
  only.active.samples = TRUE, verbose = TRUE)
}

\arguments{
  \item{phase}{
   A GHap.phase object.
  }
  \item{outfile}{
   A character value specifying the name for the output file.
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
  The function outputs files with suffix \emph{.markers}, \emph{.samples} and \emph{.phase}, as specified by \code{\link{ghap.loadphase}}.
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
# # Output data
# ghap.outphase(phase, "example")
# 
# # Reload
# phase <- ghap.loadphase("example.samples", "example.markers", "example.phase")

}

