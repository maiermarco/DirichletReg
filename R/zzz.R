.onUnload <- function(libpath){
  library.dynam.unload(chname = "DirichletReg", libpath = libpath) # unload DirichletReg's DLL
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
  
  if(missing(xdent)) xdent <- if(match.arg(type) == "message") 0L else 4L
  width <- if(getOption("width") < 40L) 40L else getOption("width")
  
  paste0(
    if(match.arg(type) == "stop") "\n" else "",
    paste0(strwrap(text, width = width, indent = xdent, exdent = xdent), collapse = "\n"),
  collapse = "")
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



pretty_list <- function(x, quotes = FALSE, oxford = TRUE, ampersand = FALSE, remove_empty_NA = TRUE){
  et <- sprintf(" %s ", if(ampersand) "&" else gettext("and")) # supports translations!
  
  if(remove_empty_NA) x <- x[!is.na(x) & nzchar(x)]
  if(quotes) x <- dQuote(x) # automatically inserts regular or fancy quotation marks according to getOption("useFancyQuotes")
  
  if(length(x) < 2L) return(x) else
  if(length(x) == 2L) return(paste(x, collapse = et)) else {
    if(oxford) et <- sprintf(",%s", et)
    return(paste(paste(head(x, -2L), collapse = ", "), paste(tail(x, 2L), collapse = et), sep = ", "))
  }
}



maxLik_code_msg <- function(maxLik_code){
  msg <- switch(as.character(maxLik_code),
    "1" = "gradient close to zero.",
    "2" = "successive function values within tolerance limit.",
    "3" = "last step could not find higher value - try changing \"steptol\".",
    "4" = "iteration limit exceeded.",
    "5" = "Infinite value.",
    "6" = "Infinite gradient.",
    "7" = "Infinite Hessian.",
    "8" = "Successive function values within relative tolerance limit.",
    "9" = "Hessian approximation cannot be improved because of gradient did not change.",
   "10" = "Lost patience: The optimizer has hit an inferior value too many times.",
  "100" = "Initial value out of range."
  )
  
  if(maxLik_code %notin% c(1L, 2L, 8L)){
    if(maxLik_code %in% c(seq_len(10L), 100L)){
      warning(gettextf(
        "Possible convergence problems!\n  Code %d: %s\n  help(\"returnCode\", package = \"maxLik\") # See \"Details\"",
        maxLik_code,
        dQuote(msg)
      ), call. = FALSE, domain = NA)
    } else {
      warning(gettextf("Unknown maxLik-Code %d. Please report this at\n  https://github.com/maiermarco/DirichletReg/issues", maxLik_code), call. = FALSE, domain = NA)
    }
  }
}
