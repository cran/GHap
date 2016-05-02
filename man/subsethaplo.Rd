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
# # Subset CEU
# CEU.ids <- haplo$id[which(haplo$pop=="CEU")]
# haplo <- ghap.subsethaplo(haplo,CEU.ids,haplo$allele.in)

}

