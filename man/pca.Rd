\name{ghap.pca}
\alias{ghap.pca}

\title{
Principal Components Analysis
}
\description{
PCA from a HapAllele-based kinship matrix.
}
\usage{
ghap.pca(haplo, K, npc = 2)
}

\arguments{
  \item{haplo}{
   A GHap.haplo object.
  }
  \item{K}{
    A HapAllele-based kinship matrix, as supplied by \code{\link{ghap.kinship}}.
  }
  \item{npc}{
    Number of principal components to be retrieved (default = 2).
  }
}

\value{
The returned object is a list with items:
 \item{eigenvec}{
  A data.frame containing the principal components of the kinship matrix.
 }
 \item{eigenval}{
  Vector with eigenvalues of the kinship matrix.
 }
  \item{propvar}{
  Vector with the proportion of variance explained by each principal component.
 }
}

\details{
Principal components are computed using singular value decomposition:

\deqn{\mathbf{K} = \mathbf{USV}'}

Since \eqn{\mathbf{K}} is a \emph{n} x \emph{n} matrix, \eqn{\mathbf{U} = \mathbf{V}}. Also, \eqn{\mathbf{S} = diag(s_i)}. Matrix \eqn{\mathbf{U}} contains the eigenvectors, whereas \eqn{s_{i}} are the eigenvalues.
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
# # Compute Kinship matrix
# K <- ghap.kinship(haplo, batchsize = 100)
# 
# ### RUN ###
# 
# # PCA analysis
# pca <- ghap.pca(haplo,K)
# 
# # Plot
# plot(x=pca$eigenvec$PC1, y=pca$eigenvec$PC2, xlab="PC1", ylab="PC2", pch="")
# pop <- pca$eigenvec$POP
# pop.col <- as.numeric(as.factor(pop))
# pop <- sort(unique(pop))
# legend("bottomright", legend = pop, col = 1:length(pop), pch = 1:length(pop), ncol = 3)
# points(x=pca$eigenvec$PC1, y=pca$eigenvec$PC2,, pch = pop.col, col = pop.col, cex = 1.2)

}
