.onAttach <- function(libname, pkgname) {
  mymsg <- "Loading required package: GHap\n\n\n"
  mymsg <- paste(mymsg,"Thanks for using GHap v2.0.0!\n")
  mymsg <- paste(mymsg,"For more information use: help(package = 'GHap')\n")
  mymsg <- paste(mymsg,"                          citation(package = 'GHap')\n")
  mymsg <- paste(mymsg,"                          browseVignettes(package = 'GHap')\n\n")
  mymsg <- paste(mymsg,"Version date: 11 Sep 2020\n\n")
  packageStartupMessage(mymsg)
}
