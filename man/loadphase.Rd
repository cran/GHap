\name{ghap.loadphase}
\alias{ghap.loadphase}

\title{
  Load phased genotype data
}
\description{
  This function loads phased genotype data and converts them into a native GHap.phase object.
}
\usage{
  ghap.loadphase(samples.file, markers.file, phase.file, verbose = TRUE)
}

\arguments{
  \item{samples.file}{
    Individual information.
  }
  \item{markers.file}{
    Variant map information.
  }
  \item{phase.file}{
    Phased genotype matrix.
  }
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}
\value{
  The returned GHap.phase object is a list with components:
    \item{chr}{
      A character value indicating chromosome identity. The current version of the package can handle only one chromosome at time.
    }
  \item{nsamples}{
    An integer value for the sample size.
  }
  \item{nmarkers}{
    An integer value for the number of markers.
  }
  \item{nsamples.in}{
    An integer value for the number of active samples.
  }
  \item{nmarkers.in}{
    An integer value for the number of active markers.
  }
  \item{pop}{
    A character vector relating chromosome alleles to populations. This information is obtained from the first column of the sample file.
  }
  \item{id}{
    A character vector mapping chromosome alleles to samples. This information is obtained from the second column of the sample file.
  }
  \item{id.in}{
    A logical vector indicating active chromosome alleles. By default, all chromosomes are set to TRUE.
  }
  \item{marker}{
    A character vector containing marker names. This information is obtained from the second column of the marker map file.
  }
  \item{marker.in}{
    A logical vector indicating active markers. By default, all markers are set to TRUE.
  }
  \item{bp}{
    A numeric vector with marker positions. This information is obtained from the third column of the marker map file.
  }
  \item{A0}{
    A character vector with reference alleles. This information is obtained from the fourth column of the marker map file.
  }
  \item{A1}{
    A character vector with alternative alleles. This information is obtained from the fifth column of the marker map file.
  }
  \item{phase}{
    A \code{\link{big.matrix}} object containing the phased genotype matrix.
  }
  
  The supported format is composed of three files with suffix:
    \itemize{
      \item \strong{.samples}: space-delimited file without header containing two columns: Population and ID. Please notice that the Population column serves solely for the purpose of grouping samples, so the user can define any arbitrary family/cluster/subgroup and use as a "population" tag. 
      \item \strong{.markers}: space-delimited file without header containing five columns: Chromosome, Marker, Position (in bp), Reference Allele (A0) and Alternative Allele (A1). Markers should be on a single chromosome and sorted by position. Repeated positions are tolerated, but a warning message is given when the data is loaded.
      \item \strong{.phase}: space-delimited file without header containing the phased genotype matrix. The dimension of the matrix is expected to be \emph{m x 2n}, where \emph{m} is the number of markers and \emph{n} is the number of individuals (i.e., two columns per individual, representing the two phased chromosome alleles). Alleles must be coded as 0 and 1. No missing values are allowed, since imputation is assumed to be part of the phasing procedure.
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
# 
# ### RUN ###
# 
# # Load data
# phase <- ghap.loadphase("human.samples", "human.markers", "human.phase")

}

