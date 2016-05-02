\name{ghap.makefile}
\alias{ghap.makefile}

\title{
 Create example input files
}
\description{
 Create example files to test the package.
}
\usage{
 ghap.makefile() 
}

\details{
This function copies the following example files to the current working directory:

\emph{human.phase}

\emph{human.markers}

\emph{human.samples}

For details about the input format, see \code{\link{ghap.loadphase}}. The dataset was extracted from the reference phased data available at the IMPUTE2 (Howie et al., 2009) software website (\url{https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference}).

The genotypes derive from the International HapMap Project Phase 3 (The International HapMap 3 Consortium, 2010), and comprise 1,011 subjects (from 11 populations) and 20,000 SNPs (randomly sampled from chromosome 2) mapped to the NCBI build 36 (hg18) assembly.

}

\author{
Marco Milanesi <marco.milanesi.mm@gmail.com>
}

\examples{
# Copy the example data in the current working directory
ghap.makefile() 
}

\references{
B. N. Howie, P. Donnelly, and J. Marchini. A flexible and accurate genotype imputation method for the next generation of genome-wide association studies. PLOS Genet. 2009. 5, e1000529.

The International HapMap 3 Consortium. Integrating common and rare genetic variation in diverse human populations. Nature. 2010. 467, 52-58.
}
