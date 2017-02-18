.onAttach <- function(libname, pkgname) {
  mymsg <- "Loading required package: GHap\n\n\n"
  mymsg <- paste(mymsg,"Thanks for using GHap v1.2.2!\n")
  mymsg <- paste(mymsg,"For more information use: help(package = 'GHap')\n")
  mymsg <- paste(mymsg,"                          browseVignettes(package = 'GHap')\n\n")
  mymsg <- paste(mymsg,"Citation:\n")
  mymsg <- paste(mymsg,"  Authors: Utsunomiya YT, Milanesi M, Utsunomiya ATH, Ajmone-Marsan P and Garcia JF.\n")
  mymsg <- paste(mymsg,"  Title: GHap: An R package for Genome-wide Haplotyping.\n")
  mymsg <- paste(mymsg,"  Journal: Bioinformatics 2016, 32(18):2861-2862.\n")
  mymsg <- paste(mymsg,"  DOI: 10.1093/bioinformatics/btw356.\n\n")
  mymsg <- paste(mymsg,"Last update:\n")
  mymsg <- paste(mymsg,"  18 Feb 2017\n\n")
  mymsg <- paste(mymsg,"Help us improve the package by reporting bugs and other issues to:\n")
  mymsg <- paste(mymsg,"  ytutsunomiya@gmail.com\n")
  mymsg <- paste(mymsg,"  marco.milanesi.mm@gmail.com\n\n")
  packageStartupMessage(mymsg)
}
