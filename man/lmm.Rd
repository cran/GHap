\name{ghap.lmm}
\alias{ghap.lmm}

\title{
Frequentist linear mixed model
}
\description{
Maximum likelihood estimates for the parameters of a linear mixed model.
}
\usage{
ghap.lmm(fixed,random,weights=NULL,env.eff=FALSE,data,K,verbose=TRUE)
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
  A numeric value for the maximum likelihood estimate of the variance of correlated random effects.
 }
  \item{varp}{
  A numeric value for the maximum likelihood estimate of the variance of permanent environmental effects. This value is suppressed if env.eff=FALSE.
 }
 \item{vare}{
  A numeric value for the maximum likelihood estimate of the variance of residual variance.
 }
 \item{h2}{
  A numeric value for the maximum likelihood estimate of the variance explained by correlated random effects only.
 }
  \item{H2}{
  A numeric value for the maximum likelihood estimate of the variance explained by random effects. This value is suppressed if env.eff=FALSE.
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

The function uses a frequentist framework to fit the following linear mixed model:

\deqn{\mathbf{y} = \mathbf{Xb} + \mathbf{Zu} + \mathbf{Zp} + \mathbf{e}}

where \eqn{\mathbf{X}} is a matrix relating \eqn{\mathbf{y}} to the vector of fixed effects \eqn{\mathbf{b}}, \eqn{\mathbf{Z}} is an incidence matrix relating \eqn{\mathbf{y}} to random effects \eqn{\mathbf{u}} and \eqn{\mathbf{p}}, and \eqn{\mathbf{e}} is the vector of residuals. The likelihood of the data is assumed:

\deqn{\mathbf{y} \mid \mathbf{b},\mathbf{u},\mathbf{p},\sigma_{u}^{2},\sigma_{p}^{2},\sigma_{e}^2 \sim N(\mathbf{Xb},\mathbf{V})}

where \eqn{\mathbf{V} = \mathbf{ZKZ}'\sigma_{u}^2 + \mathbf{ZZ}'\sigma_{p}^2 + \mathbf{W}\sigma_{e}^2}, \eqn{\mathbf{K}} is a covariance matrix for \eqn{\mathbf{u}}, \eqn{\sigma_{u}^{2}} and \eqn{\sigma_{p}^{2}} are the variances of \eqn{\mathbf{u}} and \eqn{\mathbf{p}}, respectively, \eqn{\mathbf{W}} is a residual covariance matrix and \eqn{\sigma_{e}^{2}} is the residual variance. The current implementation assumes \eqn{\mathbf{W} = diag(w_i)}. More details about the maximization algorithm can be found in our vignette.

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
# # Quantitative trait with 50\% heritability
# # One major haplotype accounting for 30\% of the genetic variance
# sim <- ghap.simpheno(haplo = haplo, K = K, h2 = 0.5, g2 = 0.3, major = 1000,seed=1988)
# 
# 
# ### RUN ###
# 
# #Continuous model
# model <- ghap.lmm(fixed = phenotype ~ 1, random = "individual", data = sim$data, K = K)
# model$h2
# plot(model$u,sim$u, ylab="True Breeding Value", xlab="Estimated Breeding Value")
# cor(model$u,sim$u)


}
