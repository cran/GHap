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

Let \eqn{\mathbf{M}} be the centered \emph{N} x \emph{H} matrix of HapGenotypes, where \emph{N} is the number of individuals and \emph{H} is the number of HapAlleles. The HapAllele covariance among individuals is computed as:

\deqn{\mathbf{K} = q\mathbf{MDM}'}

where \eqn{\mathbf{D} = diag(d_i)}, \eqn{d_i} is the weight of HapAllele \emph{i} (default \eqn{d_i = 1}), and \eqn{q} is a scaling factor defined as \eqn{tr(\mathbf{MDM}')^{-1}M}. This is a generalization of the SNP-based genomic relationship matrix (VanRaden, 2008).

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
#
# ### RUN ###
# 
# # Exclude minor alleles and singletons
# hapstats <- ghap.hapstats(haplo, ncores = 2)
# haplo <- ghap.subsethaplo(haplo,ids=haplo$id,alleles = hapstats$TYPE \%in\% c("REGULAR","MAJOR"))
#
# # Compute Kinship matrix
# K <- ghap.kinship(haplo, batchsize = 100)

}
