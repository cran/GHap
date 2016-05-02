\name{ghap.haplotyping}
\alias{ghap.haplotyping}

\title{
Haplotype genotypes
}
\description{
Generate matrix of HapGenotypes for user-defined blocks.
}
\usage{
ghap.haplotyping(phase, blocks, outfile, freq = 0.05,
 batchsize = 500, ncores = 1, verbose = TRUE)

}

\arguments{
  \item{phase}{
   A GHap.phase object.
  }
  \item{blocks}{
   A data frame containing block boundaries, such as supplied by the \code{\link{ghap.blockgen}} function.
  }
  \item{outfile}{
   A character value specifying the name for the output files.
  }
  \item{freq}{
   A numeric value specifying the minimum haplotype allele frequency (default = 0.05).
  }
  \item{batchsize}{
   A numeric value controlling the number of haplotype blocks to be processed and written to output at a time (default = 500).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{

  The function outputs three files with suffix:
  \itemize{
    \item \strong{.hapsamples}: space-delimited file without header containing two columns: Population and Individual ID.
    \item \strong{.hapalleles}: space-delimited file without header containing five columns: Block Name, Chromosome, Start and End Position (in bp), and HapAllele.
    \item \strong{.hapgenotypes}: space-delimited file without header containing the HapGenotype matrix (coded as 0, 1 or 2 copies of the HapAllele). The dimension of the matrix is \emph{m x n}, where \emph{m} is the number of HapAlleles and \emph{n} is the number of individuals.

}
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
# ### RUN ###
# 
# # Generate matrix of haplotype genotypes
# ghap.haplotyping(phase, blocks, batchsize = 100, ncores = 2, freq = 0.05, outfile = "example")

}
