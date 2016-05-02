\name{ghap.profile}
\alias{ghap.profile}

\title{
Haplotype allele profile
}
\description{
Given a data.frame of user-defined haplotype allele scores, compute individual profiles.
}
\usage{
ghap.profile(score, haplo, only.active.samples = TRUE, ncores = 1)
}

\arguments{
  \item{score}{
   A data.frame containing six columns: BLOCK, CHR, BP1, BP2, ALLELE, SCORE.
  }
  \item{haplo}{
   A GHap.haplo object.
  }
  \item{only.active.samples}{
    A logical value specifying whether calculations should be reported only for active samples (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
}

\value{
The function returns a data.frame with the following columns:
 \item{POP}{
 Population ID.
 }
 \item{ID}{
 Individual name.
 }
 \item{PROFILE}{
 Individual profile.
 }
}

\details{
The profile for each individual is calculated as sum(genotypes*score), where genotypes is a vector of number of copies of each haplotype allele and score is a vector of user-defined scores for each haplotype allele. By default, if scores are provided for only a subset of the haplotype alleles, the missing alleles scores will be set to zero. This function has the same spirit as the profiling routine implemented in the \emph{score} option in PLINK (Purcell et al., 2007; Chang et al., 2015).
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
Marco Milanesi <marco.milanesi.mm@gmail.com>
}

\references{
C. C. Chang et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience. 2015. 4, 7.

S. Purcell et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am. J. Hum. Genet. 2007. 81, 559-575.
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
# # Create a score data.frame
# score <- NULL
# score$BLOCK <- haplo$block
# score$CHR <- haplo$chr
# score$BP1 <- haplo$bp1
# score$BP2 <- haplo$bp2
# score$ALLELE <- haplo$allele
# set.seed(1988)
# score$SCORE <- rnorm(length(score$ALLELE))
# score <- data.frame(score,stringsAsFactors = FALSE)
# 
# # Compute profiles
# profile <- ghap.profile(score, haplo, ncores = 2)

}
