\name{ghap.hapstats}
\alias{ghap.hapstats}

\title{
Haplotype statistics
}
\description{
HapAllele summary statistics.
}
\usage{
ghap.hapstats(haplo, alpha = c(1,1), only.active.samples = TRUE,
 only.active.alleles = TRUE, ncores = 1)
}

\arguments{
  \item{haplo}{
    A GHap.haplo object.
  }
  \item{alpha}{
   A numeric vector of size 2 specifying the shrinkage parameters for the expected-to-observed homozygotes ratio. Default is c(1,1).
  }
  \item{only.active.samples}{
    A logical value specifying whether only active samples should be included in the output (default = TRUE).
  }
  \item{only.active.alleles}{
    A logical value specifying whether only active haplotype alleles should be included in the output (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
}

\value{
 A data frame with columns:
 \item{BLOCK}{
  Block alias.
 }
 \item{CHR}{
  Chromosome name.
 }
 \item{BP1}{
  Block start position.
 }
 \item{BP2}{
  Block end position.
 }
 \item{ALLELE}{
  Haplotype allele identity.
 }
 \item{N}{
  Number of observations for the haplotype.
 }
 \item{FREQ}{
  Haplotype frequency.
 }
 \item{O.HOM}{
  Observed number of homozygotes.
 }
 \item{O.HET}{
  Observed number of heterozygotes.
 }
 \item{E.HOM}{
  Expected number of homozygotes.
 }
 \item{RATIO}{
  Shrinkage expected-to-observed ratio for the number of homozygotes.
 }
 \item{BIN.logP}{
  log10(1/P) or -log10(P) for Hardy-Weinberg equilibrium assuming number of homozygotes follows a Binomial distribution.
 }
 \item{POI.logP}{
  log10(1/P) or -log10(P) for Hardy-Weinberg equilibrium assuming number of homozygotes follows a Poisson distribution.
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
# # Generate matrix of haplotype genotypes
# ghap.haplotyping(phase, blocks, batchsize = 100, ncores = 2, freq = 0.05, outfile = "example")
# 
# # Load haplotype genotypes
# haplo <- ghap.loadhaplo("example.hapsamples", "example.hapalleles", "example.hapgenotypes")
# 
# 
# ### RUN ###
# 
# # Compute haplotype allele statistics
# hapstats <- ghap.hapstats(haplo, ncores = 2)

}
