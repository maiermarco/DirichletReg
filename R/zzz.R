blank.trim <- function(x){
  x_split <- unlist(strsplit(x, "^\\s+|\\s+$")) # trim leading/trailing space
  paste(x_split[x_split != ""], collapse = " ") # combine w/o empty char. elements
}



deparse_nocutoff <- function(expr){ # deparse, paste, remove multiple blanks
  gsub("[[:space:]]{2,}", " ", paste(deparse(expr), collapse = " "))
}



inv.logit <- function(x){
  log(x) - log(1.0 - x)
}



swrap <- function(text, type = c("stop", "warning", "message"), xdent){
  if(missing(xdent)) xdent <- ifelse(match.arg(type) == "message", 0, 4)
  width <- ifelse(getOption("width") < 40L, 40L, getOption("width"))

  paste(
    ifelse(match.arg(type) == "stop", "\n", ""),
    paste(strwrap(text, width = width, indent = xdent, exdent = xdent), sep = "", collapse = "\n"),
  sep = "", collapse = "")
}



na.delete <- function(x){
  if(is.null(dim(x))){
    return( x[!is.na(x)] )
  } else {
    return( x[ suppressWarnings( rowSums(is.na(x)) ) == 0L , ] )
  }
}



make.symmetric <- function(x){
  if(nrow(x) != ncol(x)) stop("x must be a square matrix")

  cell.ind <- which(is.na(x), arr.ind = TRUE)

  x[cell.ind] <- x[cell.ind[, 2:1]]

  return(x)
}
