\name{ghap.haplotyping}
\alias{ghap.haplotyping}

\title{
Haplotype genotypes
}
\description{
Generate matrix of HapGenotypes for user-defined blocks.
}
\usage{
ghap.haplotyping(phase, blocks, outfile, freq = c(0,1), drop.minor = FALSE,
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
   A numeric vector of length 2 specifying the range of haplotype allele frequency to be included in the output. Default is c(0,1), which includes all alleles.
  }
  \item{drop.minor}{
   A logical value specfying whether the minor allele should be excluded from the output (default = FALSE).
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
# # Subset data - markers with maf > 0.05
# maf <- ghap.maf(phase, ncores = 2)
# markers <- phase$marker[maf > 0.05]
# phase <- ghap.subsetphase(phase, unique(phase$id), markers)
# 
# # Generate blocks of 5 markers sliding 5 markers at a time
# blocks.mkr <- ghap.blockgen(phase, windowsize = 5, slide = 5, unit = "marker")
#
# 
# ### RUN ###
# 
# # Generate matrix of haplotype genotypes
# ghap.haplotyping(phase, blocks.mkr, batchsize = 100, ncores = 2, outfile = "human")

}
