\name{ghap.assoc}
\alias{ghap.assoc}

\title{
Association analysis for HapAlleles and HapBlocks
}
\description{
Given a GHap.blmm object (as supplied by the  \code{\link{ghap.blmm}} function), this function returns association statistics for HapAllele or HapBlocks.
}
\usage{
ghap.assoc(blmm, haplo, type = "HapAllele", gc = TRUE, 
           only.active.alleles = TRUE, ncores = 1)
}

\arguments{
  \item{blmm}{
   A GHap.blmm object, such as supplied by the \code{\link{ghap.blmm}} function.
  }
  \item{haplo}{
   A GHap.haplo object.
  }
  \item{type}{
    A character value specifying whether association analysis should be carried on alleles or blocks. If type="HapAlleles" (default), association results are based on chi-squared test computed from least squares regression. If type="HapBlock", HapAlleles within each block are fitted together and an F test is reported.
  }
  \item{only.active.alleles}{
    A logical value specifying whether calculations should be reported only for active haplotype alleles (default = TRUE).
  }
  \item{gc}{
   A logical value specifying whether genomic control should be performed (default = TRUE). Currently, this option does not take effect on HapBlocks.
  }
  \item{ncores}{
   A numeric value specifying the number of cores to be used in parallel computations (default = 1).
  }
}

\value{
If type="HapAllele", the function returns a data.frame with the following columns:
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
 \item{BETA}{
  BLUE for the fixed effect of the haplotype allele.
 }
 \item{SE}{
  Standard error for the fixed effect.
 }
 \item{FREQ}{
  Frequency of the haplotype allele.
 }
 \item{CHISQ.OBS}{
  Observed value for the test statistics. If gc = TRUE (default), these values are scaled by the inflation factor. Inflation is computed through regression of observed quantiles onto expected quantiles. In order to avoid overestimation, only HapAlleles with test statistics within three standard deviations from the mean are used to compute the inflation factor.
 }
 \item{CHISQ.EXP}{
  Expected values for the test statistics.
 }
 \item{logP}{
  log10(1/P) or -log10(P) for the fixed effect.
 }
 
 If type="HapBlock", the function returns a data.frame with the following columns:
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
 \item{N.ALLELES}{
  Number of HapAlleles per block.
 }
 \item{F.TEST}{
  F statistic for HapBlock association.
 }
 \item{logP}{
  log10(1/P) or -log10(P) for the F statistic.
 }

}

\details{
Consider the residuals from the linear mixed model analysis:

\deqn{\mathbf{\hat{e}} = \mathbf{y} - \mathbf{X\hat{b}} - \mathbf{Z\hat{u}}}


These residuals can be viewed as adjusted records accounting for covariates and polygenic effects, i.e. \eqn{\mathbf{y}_{adj} = \mathbf{\hat{e}}}. We can then conduct least squares regression to test each HapAllele at a time for association with phenotypes. The fixed effect, error variance and test statistic of a given HapAllele are estimated as:

\deqn{\hat{a}_i = (\mathbf{h}_i'\mathbf{h}_i)^{-1}\mathbf{h}_i'\mathbf{y}_{adj}}
\deqn{VAR(\hat{a}_i) = (\mathbf{h}_i'\mathbf{h}_i)^{-1}\hat{\sigma}_e^2}
\deqn{t_i^2 = \frac{\hat{a}_i^2}{VAR(\hat{a}_i)}}

Under the null hypothesis that the regression coefficient is zero \eqn{t_i^2 \sim \chi^2(\nu = 1)}. Given a GHap.blmm object, the function regresses adjusted records on each HapAllele. As the linear mixed model admits repeated measures through \eqn{\mathbf{Z}}, the adjusted records are dully mapped against the vector of HapGenotypes.

The user must be aware of two known caveats associated with this approach. First, by pre-adjusting records instead of estimating HapAllele effects based on generalized least squares equations we ignore covariance structure and therefore bias the estimates downwards (Svishcheva et al., 2012). Second, each HapAllele being tested is also included in the kinship matrix, such that the HapAllele is included twice in the model: as fixed and random effect. This problem is known as proximal contamination (Listgarten et al., 2012). In the first case, we can use genomic control to recover p-values to an unbiased scale (Devlin and Roeder, 1999; Amin et al., 2007). However, not much can be done regarding the estimates of the effects. As a general recommendation, if the user is only interested in the p-values, the ghap.assoc analysis should be sufficient. When effect estimates are of interest, the user can include the candidate HapAllele as a fixed effect in the full model in ghap.blmm. For the second case, a leave-one-chromosome-out (LOCO analysis) procedure can mitigate proximal contamination (Yang et al., 2014).

As an alternative to the HapAllele test, the function also offers a HapBlock test. If the argument type is supplied with the option "HapBlock" instead of the default "HapAllele", all \emph{m} HapAlleles within a block are fitted simultaneously:

\deqn{\hat{\mathbf{a}} = (\mathbf{H}'\mathbf{H})^{-1}\mathbf{H}'\mathbf{y}_{adj}}

where \eqn{\mathbf{H}} is the \emph{n} x \emph{m} sub-HapGenotypes matrix for the HapAlleles within the block. Then, given the explained (\eqn{ESS}) and the residual (\eqn{RSS}) sum of squares, an F test is computed instead of the chi-squared test:

\deqn{F = \frac{ESS/\nu_1}{RSS/\nu_2}}
\deqn{ESS = (\mathbf{H}\hat{\mathbf{a}} - \hat{\mu}_y)'(\mathbf{H}\hat{\mathbf{a}} - \hat{\mu}_y)}
\deqn{RSS = (\mathbf{y}_{adj} - \mathbf{H}\hat{\mathbf{a}})'(\mathbf{y}_{adj} - \mathbf{H}\hat{\mathbf{a}})}
\deqn{\hat{\mu}_y = \frac{1}{N} \sum_{i=1}^{n} y_i}
\deqn{\nu_1 = m}
\deqn{\nu_2 = n - m}

Under the null hypothesis that the variance due to the HapBlock is zero \eqn{F \sim F(\nu_1,\nu_2)}.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
N. Amin et al. A Genomic Background Based Method for Association Analysis in Related Individuals. PLoS ONE. 2007. 2:e1274. 

B. Devlin and K. Roeder. Genomic control for association studies. Biometrics. 1999. 55:997-1004.

J. Listgarten et al. Improved linear mixed models for genome-wide association studies. Nat. Methods. 2012. 9:525-526.

G. R. Svishcheva et al. Rapid variance components-based method for whole-genome association analysis. Nat Genet. 2012. 44:1166-1170.

J. Yang et al. Advantages and pitfalls in the application of mixed-model association methods. Nat. Genet. 2014. 46: 100-106.
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
# #Continuous model
# model <- ghap.blmm(fixed = phenotype ~ 1, random = "individual", data = sim$data, K = K)
# 
# ### RUN ###
# 
# #HapAllele GWAS
# gwas.allele <- ghap.assoc(blmm = model, haplo = haplo, type = "HapAllele", gc = TRUE, ncores=2)
# gwas.allele$POS <- (gwas.allele$BP1+gwas.allele$BP2)/2e+6
# plot(gwas.allele$POS,gwas.allele$logP, xlab="Chromosome 2 position (in Mb)", 
#      ylab=expression(-log[10](p)), pch=20, col="#471FAA99")
# abline(v=(haplo$bp1[1000]+haplo$bp2[1000])/2e+6,lty=3)
# 
# #HapBlock GWAS
# gwas.block <- ghap.assoc(blmm = model, haplo = haplo, type = "HapBlock", ncores=2)
# gwas.block$POS <- (gwas.block$BP1+gwas.block$BP2)/2e+6
# plot(gwas.block$POS,gwas.block$logP, xlab="Chromosome 2 position (in Mb)",
#      ylab=expression(-log[10](p)), pch=20, col="#471FAA99")
# abline(v=(haplo$bp1[1000]+haplo$bp2[1000])/2e+6,lty=3)

}
