\name{ghap.subsethaplo}
\alias{ghap.subsethaplo}

\title{
 Subset GHap.hapo object
}

\description{
 This function takes a list of alleles and individuals and subsets a GHap.haplo object.
}
\usage{
ghap.subsethaplo(haplo, ids, alleles, verbose = TRUE)
}

\arguments{
  \item{haplo}{
   A GHap.haplo object.
  }
  \item{ids}{
   Character vector of individual names to keep.
  }
  \item{alleles}{
   Logical vector indicating alleles to be set to active (TRUE) or inactive (FALSE).
  }
  \item{verbose}{
    A logical value specfying whether log messages should be printed (default = TRUE).
  }
}

\value{
 The returned GHap.haplo object (as described in the documentation for the \code{\link{ghap.loadhaplo}} function) is the same as the one used in the haplo argument. However, individuals not included in the ids vector are set to FALSE (i.e., inactivated) in the haplo$samples.in vector. The vector provided in the alleles argument replaces the haplo$allele.in vector in the new GHap.haplo object. This procedure avoids expensive subsetting operations by simply flagging which haplotype alleles and individuals should be used in downstream analyses.
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
# # Randomly select 500 individuals
# ids <- sample(x = haplo$id, size = 500, replace = FALSE)
#
# #Subset data
# haplo.sub <- ghap.subsethaplo(haplo,ids,haplo$allele.in)

}

