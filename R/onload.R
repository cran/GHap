.onAttach <- function(libname, pkgname) {
  mymsg <- "Loading required package: GHap\n\n\n"
  mymsg <- paste(mymsg,"Thanks for using GHap v1.2.1!\n")
  mymsg <- paste(mymsg,"For more information use: browseVignettes('GHap')\n\n")
  mymsg <- paste(mymsg,"Citation:\n")
  mymsg <- paste(mymsg,"  Authors: Utsunomiya YT, Milanesi M, Utsunomiya ATH, Ajmone-Marsan P and Garcia JF.\n")
  mymsg <- paste(mymsg,"  Title: GHap: An R package for Genome-wide Haplotyping.\n")
  mymsg <- paste(mymsg,"  Journal: Bioinformatics.\n")
  mymsg <- paste(mymsg,"  DOI: 10.1093/bioinformatics/btw356.\n\n")
  packageStartupMessage(mymsg)
}