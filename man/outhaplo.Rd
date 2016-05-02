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
# # Subset data - randomly select 3000 markers with maf > 0.02
# maf <- ghap.maf(phase, ncores = 2)
# set.seed(1988)
# markers <- sample(phase$marker[maf > 0.02], 3000, replace = FALSE)
# phase <- ghap.subsetphase(phase, unique(phase$id), markers)
# rm(maf,markers)
# 
# # Generate block coordinates based on windows of 10 markers, sliding 5 marker at a time
# blocks <- ghap.blockgen(phase, 10, 5, "marker")
# 
# # Generate matrix of haplotype genotypes
# ghap.haplotyping(phase, blocks, batchsize = 100, ncores = 2, freq = 0.05, outfile = "example")
# 
# # Load haplotype genotypes
# haplo <- ghap.loadhaplo("example.hapsamples", "example.hapalleles", "example.hapgenotypes")
# 
# 
# ### RUN ###
# 
# # Select CEU individuals
# CEU.ids <- haplo$id[which(haplo$pop=="CEU")]
# 
# # Select a set of random HapAlleles
# set.seed(1988)
# random.allele <- sample(x=c(TRUE, FALSE), size=haplo$nalleles, replace = TRUE)
# 
# # Subset data
# haplo.CEU <- ghap.subsethaplo(haplo,CEU.ids,random.allele)
# 
# # Output data
# ghap.outhaplo(haplo = haplo.CEU, outfile = "subset")
# 
# # Reload data
# haplo.CEU <- ghap.loadhaplo("CEU.hapsamples", "CEU.hapalleles", "CEU.hapgenotypes")


}

