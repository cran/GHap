\name{ghap.lmm}
\alias{ghap.lmm}

\title{
Linear mixed model
}
\description{
Linear mixed model fitting for fixed effects, random effects and variance components.
}
\usage{
ghap.lmm(fixed, random, covmat = NULL, data, weights = NULL, family = "gaussian", 
         REML = TRUE, verbose = TRUE)
}

\arguments{
  \item{fixed}{
   Formula describing the fixed effects part of the model, e.g. y ~ a + b + c ... If the model does not include any covariate simply state the response variable with an intercept, i.e. y ~ 1.
  }
  \item{random}{
   Formula describing the random effects part of the model, e.g., ~ x + w + z.
  }
  \item{covmat}{
   A list of covariance matrices for each group of random effects. If a matrix is not defined for a given group, an identity matrix will be used.
 }
  \item{data}{
   A dataframe containing the data.
  }
  \item{weights}{
   A numeric vector with weights for observations. These weights are treated as diagonal elements of the inverse weight matrix. If not supplied, the analysis is carried out assuming all observations are equally important.
 }
 \item{family}{
  A GLM family, see \code{\link{glm}} and \code{\link{family}}. Default is "gaussian".
 }
 \item{REML}{
  A logical value specifying whether the likelihood should be restricted regarding fixed effects (default = TRUE). Only relevant for the gaussian family with identity link function.
 }
 \item{verbose}{
   A logical value specifying whether log messages should be printed (default = TRUE).
 }
}

\value{
The returned GHap.lmm object is a list with the following items:
 \item{fixed}{
  A numeric vector containing the fixed effects.
 }
  \item{random}{
  A numeric vector containing the random effects.
 }
  \item{vcp}{
  A numeric vector with variance components.
 }
  \item{residuals}{
  A numeric vector containing residuals.
 }
  \item{lme4}{
  An object of class \code{\link{merMod}}.
 }


}

\details{

The function fits mixed models with correlated random effects using Cholesky factorization of covariance matrices. The default behaviour is to use the REstricted Maximum Likelihood (REML) algorithm implemented in \code{\link{lmer}} assuming a Gaussian family and an identity link function. However, regular Maximum Likelihood (ML) fit can be specified by setting the REML argument to FALSE. Additionally, generalized linear mixed models (GLMM) can be fit by specifying a different family and link function.

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
# # Subset common haplotypes in Europeans
# EUR.ids <- haplo$id[haplo$pop \%in\% c("TSI","CEU")]
# haplo <- ghap.subsethaplo(haplo,EUR.ids,rep(TRUE,times=haplo$nalleles))
# hapstats <- ghap.hapstats(haplo, ncores = 2)
# common <- hapstats$TYPE \%in\% c("REGULAR","MAJOR") &
#  hapstats$FREQ > 0.05 &
#  hapstats$FREQ < 0.95
# haplo <- ghap.subsethaplo(haplo,EUR.ids,common)
#
# #Compute relationship matrix
# K <- ghap.kinship(haplo, batchsize = 100)
# 
# # Quantitative trait with 50\% heritability
# # Unbalanced repeated measurements (0 to 30)
# # Two major haplotypes accounting for 50\% of the genetic variance
# myseed <- 123456789
# set.seed(myseed)
# major <- sample(which(haplo$allele.in == TRUE),size = 2)
# g2 <- runif(n = 2, min = 0, max = 1)
# g2 <- (g2/sum(g2))*0.5
# sim <- ghap.simpheno(haplo, kinship = K, h2 = 0.5, g2 = g2, nrep = 30,
#                      balanced = FALSE, major = major, seed = myseed)
# 
# #Fit model using REML
# model <- ghap.lmm(fixed = phenotype ~ 1, random = ~ individual,
#                   covmat = list(individual = K), data = sim$data)
# 
# #Estimated heritability and repeatability
# model$vcp/sum(model$vcp)
# 
# #True versus estimated breeding values
# plot(model$random$individual,sim$u,xlab="Estimated BV",ylab="True BV"); abline(0,1)
# summary(lm(sim$u ~ as.numeric(model$random$individual)))

}

\references{
D. Bates et al. Fitting Linear Mixed-Effects Models Using lme4. J. Stat. Soft., 67:1-48.

A. I. Vazquez. Technical note: An R package for fitting generalized linear mixed models in animal breeding. J. Anim. Sci. 2010. 88, 497-504.
}