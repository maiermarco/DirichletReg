.onAttach <- function(libname, pkgname){
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("This is ",paste(pkgname, version))
  packageStartupMessage(pkgname, " is BETA software\nPlease report any bugs to marco.maier@wu.ac.at")
}

swrap <- function(text, type=c("stop","warning","message"), xdent){
  if(missing(xdent)) xdent <- ifelse(match.arg(type) == "message", 0, 4)
  width <- ifelse(getOption("width") < 40, 40, getOption("width"))

  paste(
    ifelse(match.arg(type) == "stop", "\n", ""),
    paste(strwrap(text, width=width, indent=xdent, exdent=xdent), sep="", collapse="\n"),
  sep="", collapse="")
}

na.delete <- function(x){
  return(x[which(!is.na(x), arr.ind=TRUE)])
}
