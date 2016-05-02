\name{ghap.simpheno}
\alias{ghap.simpheno}

\title{
Quantitative trait simulation using real HapGenotype data
}
\description{
Simulates phenotypes from a quantitative trait with arbitrary major HapAlleles.
}
\usage{
ghap.simpheno(haplo,K,h2,g2,major=NULL,seed=NULL)
}

\arguments{
  \item{haplo}{
   A GHap.haplo object.
  }
  \item{K}{
    Covariance matrix to be used in polygenic effects simulation.
  }
  \item{h2}{
   A numeric value specifying the heritability.
  }
  \item{g2}{
   A numeric value specifying the proportion of genetic variance explained by major HapAlleles.
  }
  \item{major}{
   A numeric vector specifying the indices of the HapAlleles to be considered as major.
  }
  \item{seed}{
    A numeric value used to set the random number generation state (default = NULL). This is useful for reproducibility of the results.
  }
}

\value{
The function returns a list with items:
  \item{h2}{
   A numeric value specifying the heritability.
  }
  \item{g2}{
   A numeric value specifying the proportion of genetic variance explained by major HapAlleles.
  }
  \item{major}{
   A numeric vector specifying the indices of the HapAlleles to be considered as major.
  }
  \item{major.effect}{
   A numeric vector containing the simulated major HapAllele effects.
  }
  \item{u}{
   A numeric vector containing breeding values.
  }
  \item{varu}{
   A numeric value corresponding to the genetic variance.
  }
  \item{vare}{
   A numeric value corresponding to the residual variance.
  }
  \item{data}{
   A data.frame containing columns: phenotype = a numeric vector containing the simulated phenotypes; individual =
   a character vector containing the IDs of the corresponding phenotypes.
  }

  
}

\details{

The simulation considers the model:

\deqn{\mathbf{y} = \mathbf{u} + \mathbf{e}}

where \eqn{\mathbf{u}} is a vector of breeding values and \eqn{\mathbf{e}} is the vector of residuals. Breeding values are assumed:

\deqn{\mathbf{u} = \mathbf{H}\mathbf{a} + \mathbf{g}}

where \eqn{\mathbf{a}} is the vector of major HapAllele effects, \eqn{\mathbf{H}} is the centered matrix of major HapGenotypes, and \eqn{\mathbf{g}} is a vector of polygenic effects. Vector \eqn{\mathbf{a}} is sampled from \eqn{N(0,1)}, and the variance due to major HapAlleles is computed as \eqn{\sigma_{major}^2 = VAR(\mathbf{H}\mathbf{a})}. Vector \eqn{\mathbf{g}} is sampled from \eqn{N(0,\mathbf{K}\sigma_{major}^2(1-g^2)/g^2)}. Finally, residual effects are sampled from \eqn{N(0,\sigma_{u}^2(1-h^2)/h^2)}.


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
# # Generate block coordinates based on windows of 10 markers, sliding 5 marker at a time
# blocks <- ghap.blockgen(phase, 10, 5, "marker")
# 
# # Generate matrix of haplotype genotypes
# ghap.haplotyping(phase, blocks, batchsize = 100, ncores = 2, freq = 0.05, outfile = "example")
# 
# # Load haplotype genotypes
# haplo <- ghap.loadhaplo("example.hapsamples", "example.hapalleles", "example.hapgenotypes")
# 
# # Compute kinship matrix
# K <- ghap.kinship(haplo, batchsize = 100)
# 
# ### RUN ###
# 
# # Quantitative trait with 50\% heritability
# # One major haplotype accounting for 30\% of the genetic variance
# sim <- ghap.simpheno(haplo = haplo, K = K, h2 = 0.5, g2 = 0.3, major = 1000,seed=1988)
# 
# # Binary trait from the previous example
# # 0 if observation is below the 70\% percentile
# # 1 otherwise
# thr <- quantile(x = sim$data$phenotype, probs = 0.7)
# sim$data$phenotype2 <- sim$data$phenotype
# sim$data$phenotype2[sim$data$phenotype < thr] <- 0
# sim$data$phenotype2[sim$data$phenotype >= thr] <- 1

}
