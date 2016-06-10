\name{ghap.kinship}
\alias{ghap.kinship}

\title{
Kinship matrix from haplotypes
}
\description{
This function computes a HapAllele-based kinship matrix from a GHap.haplo object.
}
\usage{
ghap.kinship(haplo, weights, batchsize = 500, only.active.samples = TRUE,
 only.active.alleles = TRUE, verbose = TRUE)
}

\arguments{
  \item{haplo}{
    A GHap.haplo object.
  }
  \item{weights}{
    A numeric vector providing HapAllele-specific weights.
  }
  \item{batchsize}{
   A numeric value controlling the number of haplotype alleles to be processed at a time (default = 500).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{only.active.alleles}{
    A logical value specifying whether only active haplotype alleles should be included in the output (default = TRUE).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
The function returns a \emph{n} x \emph{n} matrix of HapAllele-based kinships, where \emph{n} is the number of individuals.
}

\details{

Briefly, let \eqn{\mathbf{H}} be the centered \emph{n} x \emph{m} matrix of HapGenotypes, where \emph{n} is the number of individuals and \emph{m} is the number of HapAlleles. The HapAllele covariance among individuals is computed as:

\deqn{\mathbf{K} = q\mathbf{HDH}'}

where \eqn{\mathbf{D} = diag(d_i)}, \eqn{d_i} is the weight of HapAllele \emph{i} (default \eqn{d_i = 1}), and \eqn{q} is the inverse weighted sum of variances in the columns of \eqn{\mathbf{H}}. Notice that this is a generalization of the SNP-based genomic relationship matrix (VanRaden, 2008).

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>

Marco Milanesi <marco.milanesi.mm@gmail.com>
}


\references{
P. M. VanRaden. Efficient methods to compute genomic predictions. J. Dairy. Sci. 2008. 91:4414-4423.
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
# # Compute Kinship matrix
# K <- ghap.kinship(haplo, batchsize = 100)

}
