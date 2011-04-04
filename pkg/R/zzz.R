.onAttach <- function(libname, pkgname){
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("This is ",paste(pkgname, version))
  packageStartupMessage(pkgname, " is BETA software\nPlease report any bugs to marco.maier@wu.ac.at")
}
