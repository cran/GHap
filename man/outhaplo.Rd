\name{ghap.outhaplo}
\alias{ghap.outhaplo}

\title{
  Output GHap.haplo object to a file
}
\description{
  This function saves a GHap.haplo object to a text file using the same format specified by \code{\link{ghap.haplotyping}}.
}
\usage{
  ghap.outhaplo(haplo, outfile, only.active.markers = TRUE,
                only.active.samples = TRUE, verbose = TRUE)
}

\arguments{
  \item{haplo}{
    A GHap.haplo object.
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
  The function outputs files with suffix \emph{.hapalleles}, \emph{.hapsamples} and \emph{.hapgenotypes}, as specified by \code{\link{ghap.haplotyping}}.
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
# # Randomly select 500 individuals
# ids <- sample(x = haplo$id, size = 500, replace = FALSE)
#
# #Subset data
# haplo.sub <- ghap.subsethaplo(haplo,ids,haplo$allele.in)
#
#
# ### RUN ###
#
# # Output new GHap.haplo object
# ghap.outhaplo(haplo = haplo.sub, outfile = "humansub")

}

