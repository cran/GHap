\name{ghap.fst}
\alias{ghap.fst}

\title{
 Haplotype Fst
}

\description{
Multi-alleleic Fst computed using block summary statistics generated from \code{\link{ghap.blockstats}}.
}

\usage{
 ghap.fst(blockstats.pop1, blockstats.pop2, blockstats.tot)
}

\arguments{
  \item{blockstats.pop1}{
   A data.frame containing block statistics computed on population 1.
  }
  \item{blockstats.pop2}{
   A data.frame containing block statistics computed on population 2.
  }
  \item{blockstats.tot}{
   A data.frame containing block statistics computed on population 1 + population 2.
  }
}

\value{
The function returns a data.frame with the following columns:
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
  \item{EXP.H.pop1}{
  Expected heterozygosity in population 1.
 }
  \item{EXP.H.pop2}{
  Expected heterozygosity in population 2.
 }
  \item{EXP.H.tot}{
  Expected heterozygosity in the total population.
 }
 \item{FST}{
  Fst value.
 }

}

\details{
 This function calculates Fst (Nei, 1973) based on the formula for multi-allelic markers:
 
 \deqn{Fst = (Ht - Hs) / Ht}
 
 where \emph{Ht} is the total gene diversity (i.e., expected heterozygosity in the population) and \emph{Hs} is the subpopulation gene diversity (i.e., the average expected heterozygosity in the subpopulations).
  
}

\author{
Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>

Marco Milanesi <marco.milanesi.mm@gmail.com>
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
# # Compute haplotype allele statistics for each group
# CHB.ids <- haplo$id[which(haplo$pop=="CHB")]
# CEU.ids <- haplo$id[which(haplo$pop=="CEU")]
# haplo <- ghap.subsethaplo(haplo,CHB.ids,haplo$allele.in)
# CHB.hapstats <- ghap.hapstats(haplo,ncores = 2)
# haplo <- ghap.subsethaplo(haplo,CEU.ids,haplo$allele.in)
# CEU.hapstats <- ghap.hapstats(haplo,ncores = 2)
# haplo <- ghap.subsethaplo(haplo,c(CHB.ids,CEU.ids),haplo$allele.in)
# TOT.hapstats <- ghap.hapstats(haplo,ncores = 2)
# 
# # Compute haplotype block statistics for each group
# CHB.blockstats <- ghap.blockstats(CHB.hapstats, ncores = 2)
# CEU.blockstats <- ghap.blockstats(CEU.hapstats, ncores = 2)
# TOT.blockstats <- ghap.blockstats(TOT.hapstats, ncores = 2)
# 
# # Calculate Fst
# fst<-ghap.fst(CHB.blockstats, CEU.blockstats, TOT.blockstats)
# 
# # Plot results
# top.fst <- fst[fst$FST == max(fst$FST, na.rm=TRUE),]
# plot(
#   x = (fst$BP1+fst$BP2)/2e+6,
#   y = fst$FST, pch = "",
#   ylab = expression(paste("Haplotype ", F[ST])),
#   xlab = "Chromosome 2 (in Mb)",
#   ylim=c(0,1)
# )
# abline(v=108.7, col="gray")
# points(x = (fst$BP1+fst$BP2)/2e+6, y = fst$FST, pch = 20, col="#471FAA99")
# points(x = (top.fst$BP1+top.fst$BP2)/2e+6, y = top.fst$FST, pch = 20, col="red")
# text(x = 125, y = max(fst$FST, na.rm=TRUE), "EDAR", col="red")
# CEU.hapstats[CEU.hapstats$BLOCK == top.fst$BLOCK & CEU.hapstats$FREQ > 0,1:9]
# CHB.hapstats[CHB.hapstats$BLOCK == top.fst$BLOCK & CHB.hapstats$FREQ > 0,1:9]

}

\references{
M. Nei. Analysis of Gene Diversity in Subdivided Populations. PNAS. 1973. 70, 3321-3323.
}
