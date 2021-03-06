\name{ghap.blockstats}
\alias{ghap.blockstats}

\title{
HapBlock statistics
}
\description{
Generate HapBlock summary statistics from pre-computed HapAlleles statistics.
}
\usage{
ghap.blockstats(hapstats, ncores=1, verbose=TRUE)
}

\arguments{
  \item{hapstats}{
   A data.frame containing HapAllele statistics, as generated by the \code{\link{ghap.hapstats}} function.
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
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
 \item{EXP.H}{
  Block expected heterozygosity.
 }
 \item{N.ALLELES}{
  Number of HapAlleles per block.
 }
 
}

\details{
 For each HapBlock, the function counts the number of unique HapAlleles and computes the expected heterozygosity \eqn{ 1 - \sum {p_{i}}^2 }, where \eqn{p_{i}} is the frequency of HapAllele \emph{i}. Please notice that when HapAlleles are prunned out by frequency the block statistics can retrieve high expected heterozygosity for blocks with small number of HapAlleles.
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\examples{

# #### DO NOT RUN IF NOT NECESSARY ###
# 
# # Copy the example data in the current working directory
# exfiles <- ghap.makefile()
# file.copy(from = exfiles, to = "./")
# 
# # Compress phase data
# ghap.compress(input.file = "human", out.file = "human")
# 
# # Load phase data
# phase <- ghap.loadphase(input.file = "human")
# 
# # Generate blocks of 5 markers sliding 5 markers at a time
# blocks <- ghap.blockgen(phase, windowsize = 5, slide = 5, unit = "marker")
# 
# # Haplotyping
# ghap.haplotyping(phase = phase, blocks = blocks, outfile = "human",
#                  binary = T, ncores = 1)
#
# # Load haplotype genotypes
# haplo <- ghap.loadhaplo(input.file = "human")
#
# Subset
# ids <- haplo$id[which(haplo$pop == "CEU")]
# haplo <- ghap.subsethaplo(haplo = haplo, ids = ids, alleles = haplo$allele.in)
#
# #Compute haplotype statistics
# hapstats <- ghap.hapstats(haplo, ncores = 1)
#
#
# ### RUN ###
#
# #Compute block statistics
# blockstats <- ghap.blockstats(hapstats, ncores = 1)

}

