\name{ghap.blmm}
\alias{ghap.blmm}

\title{
Bayesian linear mixed model
}
\description{
A Gibbs sampler is used to find posterior mean estimates for the parameters of a linear mixed model.
}
\usage{
ghap.blmm(fixed,random,weights=NULL,ordinal=FALSE,env.eff=FALSE,data,K,vu=4,vp=4,ve=4,
          R2=0.5,nchain=1,nsim=500,burnin=0,thin=1,ncores=1,verbose=TRUE)
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
  \item{ordinal}{
   A logical value indicating if responses are ordered categorical values (default=FALSE).
  }
  \item{env.eff}{
   A logical value indicating if permenant environmental effects should be included (default=FALSE).
  }
  \item{data}{
   A dataframe containing the data.
  }
  \item{K}{
    A covariance matrix for random effects.
  }
  \item{vu}{
    A numeric value specifying the prior degrees of freedom for the variance of random effects (default = 4).
  }
  \item{vp}{
    A numeric value specifying the prior degrees of freedom for the variance of permanent environmental effects (default = 4).
  }
  \item{ve}{
    A numeric value specifying the prior degrees of freedom for the residual variance (default = 4).
  }
  \item{R2}{
   A numeric value specifying the prior variance explained by random effects. In case the covariance matrix is a genomic relationship matrix, this corresponds to the prior for the narrow sense heritability (default = 0.5). If permanent environmental effects are included, this becomes the prior for the broad sense heritability.
  }
  \item{nchain}{
   A numeric value indicating the number of independent Markov chains (default = 1).
  }
  \item{nsim}{
   A numeric value specifying the number of simulations to be performed by the main MCMC algorithm (default = 500).
  }
  \item{burnin}{
   A numeric value specifying the number of simulations to be performed prior to the main MCMC algorithm (default = 0).
  }
  \item{thin}{
   A numeric value specifying the thinning interval (default = 1).
  }
  \item{ncores}{
   A numeric value specifying the number of processors to be used in parallelization of independent Markov chains (default = 1).
  }
  \item{verbose}{
   A logical value specifying whether log messages should be printed (default = TRUE).
  }
}

\value{
The returned GHap.blmm object is a list with the following items:
  \item{nchain}{
   Number of independent Markov chains.
  }
  \item{nsim}{
   Number of simulations performed by the main MCMC algorithm.
  }
  \item{thin}{
   Thinning interval used by the MCMC algorithm.
  }
  \item{eff.nsim}{
   The effective number of samples (nsim/thin) used for posterior computations.
  }
 \item{b}{
  A numeric vector containing the posterior means of the fixed effects.
 }
 \item{u}{
  A numeric vector containing the posterior means for the correlated random effects.
 }
  \item{p}{
  A numeric vector containing the posterior means for the permanent environmental effects. This vector is suppressed if env.eff=FALSE.
 }
  \item{varu}{
  A numeric value for the posterior mean of the variance of correlated random effects.
 }
  \item{varp}{
  A numeric value for the posterior mean of the variance of permanent environmental effects. This value is suppressed if env.eff=FALSE.
 }
 \item{vare}{
  A numeric value for the posterior mean of the residual variance.
 }
 \item{h2}{
  A numeric value for the posterior mean of the variance explained by correlated random effects.
 }
  \item{H2}{
  A numeric value for the posterior mean of the variance explained by permanent environmental effects. This value is suppressed if env.eff=FALSE.
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
  A numeric vector containing residuals computed based on the posterior mean of the model parameters.
 }
  \item{dev}{
  Posterior mean of the deviance (-2*log-likelihood).
 }
 \item{pdev}{
  Deviance evaluated at the posterior mean of model parameters.
 }

If the model is fitted to ordered categorical data, the following items are added to the object:

  \item{liability}{
  Posterior means of liabilities.
 }
  \item{thresholds}{
  A numeric vector containing the posterior means of the thresholds in the liability scale.
 }

Additionally, if multiple independent chains are ran, the following items are included under inter-chain:

 \item{b.sd}{
  A numeric vector containing standard deviations for inter-chain estimates of the fixed effects.
 }
 \item{u.sd}{
  A numeric vector containing standard deviations for inter-chain estimates of the correlated random effects.
 }
  \item{p.sd}{
  A numeric vector containing standard deviations for inter-chain estimates of the permanent environmental effects. This vector is suppressed if env.eff=FALSE.
 }
  \item{varu.sd}{
  A numeric value containing standard deviations for inter-chain estimates of the variance of correlated random effects.
 }
  \item{varp.sd}{
  A numeric value containing standard deviations for inter-chain estimates of the variance of permanent environmental effects. This value is suppressed if env.eff=FALSE.
 }
 \item{vare.sd}{
  A numeric value containing standard deviations for inter-chain estimates of the residual variance.
 }
 \item{h2.sd}{
  A numeric value containing standard deviations for inter-chain estimates of the variance explained by the correlated random effects.
 }
 \item{H2.sd}{
  A numeric value containing standard deviations for inter-chain estimates of the variance explained by the correlated random effects and permanent environmental effects. This value is suppressed if env.eff=FALSE.
 }
 \item{dev.sd}{
  A numeric value containing standard deviations for inter-chain estimates of the deviance.
 }

}

\details{

The function uses a Bayesian framework to fit the following linear mixed model:

\deqn{\mathbf{y} = \mathbf{Xb} + \mathbf{Zu} + \mathbf{Zp} + \mathbf{e}}

where \eqn{\mathbf{X}} is a matrix relating \eqn{\mathbf{y}} to the vector of fixed effects \eqn{\mathbf{b}}, \eqn{\mathbf{Z}} is an incidence matrix relating \eqn{\mathbf{y}} to random effects \eqn{\mathbf{u}} and \eqn{\mathbf{p}}, and \eqn{\mathbf{e}} is the vector of residuals. The likelihood of the data and the prior distribution of the parameters are assumed:

\deqn{\mathbf{y} \mid \mathbf{b},\mathbf{u},\mathbf{p},\sigma_{u}^{2},\sigma_{p}^{2},\sigma_{e}^2 \sim N(\mathbf{Xb}+\mathbf{Zu},\mathbf{W}\sigma_{e}^2)}
\deqn{\mathbf{b} \propto constant}
\deqn{\mathbf{u} \mid \sigma_{u}^2 \sim N(0,\mathbf{K}\sigma_{u}^2)}
\deqn{\mathbf{p} \mid \sigma_{p}^2 \sim N(0,\mathbf{I}\sigma_{u}^2)}
\deqn{\sigma_{u}^2 \sim \chi^{-2}(\nu_{u},S_u^2)}
\deqn{\sigma_{p}^2 \sim \chi^{-2}(\nu_{p},S_p^2)}
\deqn{\sigma_{e}^2 \sim \chi^{-2}(\nu_{e},S_e^2)}

where \eqn{\mathbf{K}} is a covariance matrix for \eqn{\mathbf{u}}, \eqn{\sigma_{u}^{2}} and \eqn{\sigma_{p}^{2}} are the variances of \eqn{\mathbf{u}} and \eqn{\mathbf{p}}, respectively, \eqn{\mathbf{W}} is a residual covariance matrix and \eqn{\sigma_{e}^{2}} is the residual variance. The current implementation assumes \eqn{\mathbf{W} = diag(w_i)}. The hyper-parameters \eqn{\nu_{u}}, \eqn{\nu_{p}}, \eqn{\nu_{e}}, \eqn{S_{u}^2}, \eqn{S_{p}^2} and \eqn{S_{e}^2} are the random effects and residual variance degrees of freedom and scale parameters, respectively. 

In the case of ordered categorical data, categories are assumed to emerge from thresholds of a latent normal variable (i.e., "liability"). The function uses data augmentation to sample thresholds and observations from the underlying latent variable, which are then treated as responses in the main algorithm. More details about the MCMC algorithm can be found in our vignette.

}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\references{
P. Perez and G. de Los Campos. Genome-Wide Regression and Prediction with the BGLR Statistical Package. Genetics. 2014. 198:483-495.

D. A. Sorensen et al. Bayesian inference in threshold models using Gibbs sampling. Genet Sel Evol. 1995. 27:229-249.

C. S. Wang et al. Bayesian analysis of mixed linear models via Gibbs sampling with an application to litter size in Iberian pigs. Genet Sel Evol. 1994. 26:91-115.
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
# # Binary trait from the previous example
# # 0 if observation is below the 70\% percentile
# # 1 otherwise
# thr <- quantile(x = sim$data$phenotype, probs = 0.7)
# sim$data$phenotype2 <- sim$data$phenotype
# sim$data$phenotype2[sim$data$phenotype < thr] <- 0
# sim$data$phenotype2[sim$data$phenotype >= thr] <- 1
# 
# ### RUN ###
# 
# #Continuous model
# model <- ghap.blmm(fixed = phenotype ~ 1, random = "individual", data = sim$data, K = K)
# model$h2
# plot(model$u,sim$u, ylab="True Breeding Value", xlab="Estimated Breeding Value")
# cor(model$u,sim$u)
# 
# #Threshold model
# model <- ghap.blmm(fixed = phenotype2 ~ 1, random = "individual",
#                    ordinal = TRUE, data = sim$data, K = K)
# model$h2
# plot(model$u,sim$u, ylab="True Breeding Value", xlab="Estimated Breeding Value")
# cor(model$u,sim$u)
# model$thresholds[2]
# quantile(x = scale(sim$data$phenotype), probs = 0.7)

}
