.onUnload <- function(libpath){
  library.dynam.unload("DirichletReg", libpath) # unload the package's DLL
}



blank.trim <- function(x){
  .Deprecated(msg = "This function is unused and will be removed in future versions.")
  
  x_split <- unlist(strsplit(x, "^\\s+|\\s+$")) # trim leading/trailing space
  paste(x_split[x_split != ""], collapse = " ") # combine w/o empty char. elements # nolint
}



deparse_nocutoff <- function(expr, width.cutoff = 500L, ...){ # similar to deparse1() but removes duplicated whitespace
  gsub("[[:space:]]{2,}", "\\ ", paste(deparse(expr, width.cutoff = width.cutoff, ...), collapse = " ")) # nolint
}



inv.logit <- function(x){
  .Deprecated(new = "qlogis()")
  log(x) - log(1.0 - x)
}



swrap <- function(text, type = c("stop", "warning", "message"), xdent){
  .Deprecated(msg = "This function is unused and will be removed in future versions.")
  
  if(missing(xdent)) xdent <- ifelse(match.arg(type) == "message", 0L, 4L)
  width <- ifelse(getOption("width") < 40L, 40L, getOption("width")) # nolint
  
  paste(
    ifelse(match.arg(type) == "stop", "\n", ""),
    paste(strwrap(text, width = width, indent = xdent, exdent = xdent), sep = "", collapse = "\n"),
  sep = "", collapse = "")
}



na.delete <- function(x){
  if(is.null(dim(x))){
    return( x[!is.na(x)] )
  } else {
    return( x[rowSums(is.na(x)) == 0L, ] )
  }
}



make.symmetric <- function(x){
  .Deprecated(msg = "This function is unused and will be removed in future versions.")
  
  if(nrow(x) != ncol(x)) stop("x must be a square matrix")
  cell.ind <- which(is.na(x), arr.ind = TRUE)
  x[cell.ind] <- x[cell.ind[, 2:1]] # nolint
  return(x)
}
