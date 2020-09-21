\name{ghap.fast2phase}
\alias{ghap.fast2phase}

\title{
  Convert fastPHASE data into the GHap phase format
}
\description{
  This function takes phased genotype data in fastPHASE format and converts them into a GHap plan and compressed binary format.
}
\usage{
  ghap.fast2phase(input.files = NULL, switchout.files = NULL,
                  map.files = NULL, out.file, verbose = TRUE)
}

\arguments{

  If all input files share the same prefix, the user can use the following shortcut options:

  \item{input.files}{
    Character vector with the list of prefixes for input files.
  }
  \item{out.file}{
    Character value for the output file name.
  }

  The user can also opt to point to input files separately:
  
  \item{switchout.files}{
    Character vector containing the list of fastPHASE files
  }
  \item{map.files}{
    Character vector containing the list of map files.
  }
  
  To turn conversion progress-tracking on or off please use:

  \item{verbose}{
  	A logical value specfying whether log messages should be printed (default = TRUE).
  }

}

\details{
Currently this function handles only _switch.out files from fastPHASE v1.4.0. The map files should contain the following 5 space-delimited columns: chromosome, marker, position, allele 0 and allele 1.
}

\references{
Scheet, P., Stephens, M., 2006. A Fast and Flexible Statistical Model for Large-Scale Population Genotype Data: Applications to Inferring Missing Genotypes and Haplotypic Phase. Am. J. Hum. Genet. 78, 629–644. https://doi.org/10.1086/502802.
}

\author{
Mario Barbato <mario.barbato@unicatt.it>, Yuri Tani Utsunomiya <ytutsunomiya@gmail.com>
}

\seealso{
  \code{\link{ghap.compress}}, \code{\link{ghap.loadphase}}, \code{\link{ghap.oxford2phase}}, \code{\link{ghap.vcf2phase}}
}