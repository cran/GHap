\name{ghap.mme}
\alias{ghap.mme}

\title{
Mixed linear model solver
}
\description{
Given known variance components, Best Linear Unbiased Estimators (BLUE) of fixed effects and Best Linear Unbiased Predictors (BLUP) of random effects are obtained using Henderson's mixed model equations.
}
\usage{
ghap.mme(fixed,random,weights=NULL,env.eff=FALSE,data,K,varcomp,verbose=TRUE)
}

\arguments{
  \item{fixed}{
   Formula describing the fixed effects part of the model, e.g. y ~ a + b + c ... If the model does not include any covariate simply state the response variable with an intercept, i.e. y ~ 1.
  }
  \item{random}{
   A character value with name of the column containing labels for the random effects.
  }
  \item{weights}{
   A numeric vector with residual weights. If not supplied, all residual weights are set to 1.
  }
  \item{env.eff}{
   A logical value indicating if permanent environmental effects should be included (default=FALSE).
  }
  \item{data}{
   A dataframe containing the data.
  }
  \item{K}{
    A covariance matrix for random effects.
  }
  \item{varcomp}{
   A numeric vector specifying the known variance components, in the following order: residual, correlated (e.g., additive genetic) and uncorrelated (e.g., permanent environmental) random effects variances, respectively. If env.eff=FALSE, the third variance component is ignored.
  }
  \item{verbose}{
   A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
The returned GHap.blmm object is a list with the following items:
 \item{b}{
  A numeric vector containing the BLUE of fixed effects.
 }
 \item{u}{
  A numeric vector containing BLUP of correlated random effects.
 }
  \item{p}{
  A numeric vector containing the BLUP of permanent environmental effects. This vector is suppressed if env.eff=FALSE.
 }
  \item{varu}{
  A numeric value for the known variance of correlated random effects.
 }
  \item{varp}{
  A numeric value for the known variance of permanent environmental effects. This value is suppressed if env.eff=FALSE.
 }
 \item{vare}{
  A numeric value for the known residual variance.
 }
 \item{h2}{
  A numeric value for the known variance explained by correlated random effects only.
 }
  \item{H2}{
  A numeric value for the known variance explained by random effects. This value is suppressed if env.eff=FALSE.
 }
 \item{k}{
  A numeric vector containing the solutions for \eqn{\mathbf{K}^{-1}\mathbf{\hat{u}}}. This vector is used by the \code{\link{ghap.blup}} function.
 }
 \item{y}{
  A numeric vector containing the records used to fit the model.
 }
  \item{weights}{
  A numeric vector containing the residual weights used to fit the model.
 }
  \item{residuals}{
  A numeric vector containing residuals computed based on the BLUE and BLUP solutions.
 }
 \item{pdev}{
  Deviance evaluated at the BLUE and BLUP solutions.
 }

}

\details{
Consider the model described in \code{\link{ghap.blmm}} and \code{\link{ghap.lmm}}. In the special case where variance components are known, this function uses Henderson's mixed model equations to solve for fixed and random effects.
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
C. R. Henderson. Sire evaluation and genetic trends. Journal of Animal Science. 1973. Symposium: 10-41.
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
# ### RUN ###
# 
# # Solving for effects
# model <- ghap.mme(fixed = phenotype ~ 1, random = "individual", data = sim$data, K = K
#                   varcomp = c(sim$vare,sim$varu))
# plot(model$u,sim$u, ylab="True Breeding Value", xlab="Estimated Breeding Value")
# cor(model$u,sim$u)


}
