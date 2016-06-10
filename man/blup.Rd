\name{ghap.blup}
\alias{ghap.blup}

\title{
Convert breeding values into BLUP solutions of HapAllele effects
}
\description{
Given solutions for a mixed linear model, compute Best Linear Unbiased Predictor (BLUP) solutions for HapAllele effects.
}
\usage{
ghap.blup(blmm,haplo,weights=NULL,nperm=1,only.active.alleles = TRUE,ncores = 1)
}

\arguments{
  \item{blmm}{
   A GHap.blmm object, such as supplied by the \code{\link{ghap.blmm}} function.
  }
  \item{haplo}{
   A GHap.haplo object.
  }
  \item{weights}{
    A numeric vector providing HapAllele-specific weights.
  }
  \item{nperm}{
   Number of permutations to be performed for significance assessment (default = 1).
  }
  \item{only.active.alleles}{
    A logical value specifying whether only active haplotype alleles should be included in the calculations (default = TRUE).
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computing (default = 1).
  }
}

\value{
The function returns a dataframe with columns:
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
 \item{SCORE}{
  BLUP for the random effect of the haplotype allele.
 }
 \item{FREQ}{
  Frequency of the haplotype allele.
 }
 \item{VAR}{
  Variance in allele-specific breeding values.
 }
  \item{pVAR}{
  Proportion of variance explained by the haplotype allele.
 }
 \item{CENTER}{
  Average genotype (meaningful only for predictions with \code{\link{ghap.profile}}).
 }
 \item{SCALE}{
  A constant set to 1 (meaningful only for predictions with \code{\link{ghap.profile}}).
 }
 \item{P}{
  P-value for the permutation test. This column is suppressed if nperm < 1.
 }

}

\details{

The function uses the equation:

\deqn{\mathbf{\hat{a}} = q\mathbf{DH}'\mathbf{K}^{-1}\mathbf{\hat{u}}}

where \eqn{\mathbf{H}} is the \emph{n} x \emph{m} centered matrix of HapGenotypes observed for \emph{n} individuals and \emph{m} HapAlleles, \eqn{\mathbf{D} = diag(d_i)}, \eqn{d_i} is the weight of HapAllele \emph{i} (default \eqn{d_i = 1}), \eqn{q} is the inverse weighted sum of variances in the columns of \eqn{\mathbf{H}}, \eqn{\mathbf{K}} is the haplotype-based kinship matrix and \eqn{\hat{u}} is the vector of estimated breeding values. The permutation procedure consists in randomizing the vector \eqn{\hat{\mathbf{u}}} and computing the null statistic \eqn{max(\mathbf{\hat{a}})}. The permutation p-value is computed as the number of times the HapAllele effect was smaller than the null statistic.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
I. Stranden and D.J. Garrick. Technical note: derivation of equivalent computing algorithms for genomic predictions and reliabilities of animal merit. J Dairy Sci. 2009. 92:2971-2975.
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
# # Quantitative trait with 50\% heritability
# # One major haplotype accounting for 30\% of the genetic variance
# sim <- ghap.simpheno(haplo = haplo, K = K, h2 = 0.5, g2 = 0.3, major = 1000,seed=1988)
# 
# # Continuous model
# model <- ghap.blmm(fixed = phenotype ~ 1, random = "individual", data = sim$data, K = K)
# 
# ### RUN ###
# 
# 
# # BLUP solutions for HapAllele effects
# blup <- ghap.blup(blmm = model, haplo = haplo, ncores=2)
# blup$POS <- (blup$BP1+blup$BP2)/2e+6
# plot(blup$POS,100*blup$pVAR, xlab="Chromosome 2 position (in Mb)",
#      ylab="Haplotype variance explained (\%)", pch=20, col="#471FAA99")
# abline(v=(haplo$bp1[1000]+haplo$bp2[1000])/2e+6,lty=3)
# 
# # Permutation test
# blup <- ghap.blup(blmm = model, haplo = haplo, ncores=2, nperm = 1000)
# blup$POS <- (blup$BP1+blup$BP2)/2e+6
# plot(blup$POS,-log10(blup$P), xlab="Chromosome 2 position (in Mb)",
#      ylab=expression(-log[10](p)), pch=20, col="#471FAA99")
# abline(v=(haplo$bp1[1000]+haplo$bp2[1000])/2e+6,lty=3)
# abline(h=-log10(0.05))
# 
# # Re-weighted solutions
# w <- nrow(blup)*blup$pVAR
# K <- ghap.kinship(haplo, weights = w, batchsize = 100)
# model2 <- ghap.mme(fixed = phenotype ~ 1, random = "individual", data = sim$data, K = K, 
#                    varcomp=c(model$vare,model$varu))
# blup <- ghap.blup(blmm = model2, haplo = haplo, weights = w, ncores=2)
# blup$POS <- (blup$BP1+blup$BP2)/2e+6
# plot(blup$POS,100*blup$pVAR, xlab="Chromosome 2 position (in Mb)",
#      ylab="Haplotype variance explained (\%)", pch=20, col="#471FAA99")
# abline(v=(haplo$bp1[1000]+haplo$bp2[1000])/2e+6,lty=3)
# sim$major.effect
# blup$SCORE[1000]

}
