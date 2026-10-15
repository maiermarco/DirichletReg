print.DirichletRegData <- function(x, type = c("processed", "original"), ...){
  allowed_dddots <- c("digits", "na.print", "print.gap", "max", "width") # "quote", "useSource", "right"
  chkDots(..., allowed = allowed_dddots)
  
  type <- match.arg(type)
  
  if(type == "processed"){
    out <- unclass(x)
    attr(out, "Y.original" ) <- NULL
    attr(out, "dims"       ) <- NULL
    attr(out, "dim.names"  ) <- NULL
    attr(out, "obs"        ) <- NULL
    attr(out, "valid_obs"  ) <- NULL
    attr(out, "normalized" ) <- NULL
    attr(out, "transformed") <- NULL
    attr(out, "base"       ) <- NULL
    attr(out, "class"      ) <- NULL
  } else {
    out <- attr(x, "Y.original")
  }
  
  do.call("print", c(list("x" = out), list(...)[...names() %in% allowed_dddots]))
  
  invisible(x)
}



`[<-.DirichletRegData` <- function(x, i, j, value){
  stop("objects of class \"DirichletRegData\" cannot be edited manually.", "Manipulate your data in another object and use DR_data() instead.") # TODO
}



`[.DirichletRegData` <- function(x, i, j, drop){
  has_i <- !missing(i)
  has_j <- !missing(j)
  
  if(any(names(syscall <- sys.call()) %notin% c("", "drop"))){ # && !isS4(x) # should be if(any(... %notin% ...)) # TODO
    warning("named arguments other than 'drop' are discouraged") # TODO
  }
  syscall[which(names(syscall) != "")] <- NULL # remove unknown arguments and drop
  
  len_syscall <- length(syscall)
  
  if(len_syscall == 3L){                # call with 2 arguments (x, i)
    if(!has_i) return(x) else           # covers: x[]            --> as is
      return(unclass(x)[i])             # covers: x[vec], x[mat] --> vector (unclassed)
  } else if(len_syscall == 4L){         # call with 3 arguments (x, i, j)
    if(!has_i && !has_j) return(x) else # covers: x[,]         --> as is
    if((has_i && is.matrix(i)) || (has_j && is.matrix(j))) stop("attempt to index rows/columns with a matrix.") else { # covers: x[m,], x[,m] --> error! # TODO
      if(!has_i) i <- seq_len(nrow(x)) # if i is not present, set TRUE for all rows
      if(!has_j) j <- seq_len(ncol(x)) # if j is not present, set TRUE for all columns
      
      # check logical indices
      if((is.logical(i) && (length(i) != nrow(x))) || (is.logical(j) && (length(j) != ncol(x)))) stop("if logical, lengths of the indices must exactly match the number of rows or columns of x.") # TODO
      # convert logical indices to index numbers
      if(is.logical(i)) i <- which(i)
      if(is.logical(j)) j <- which(j)
      
      # if row indices are names, check and convert to index numbers
      if(is.character(i)){
        x_rownames <- rownames(x)
        if(is.null(x_rownames)) stop("x does not have row names. Indexing by names failed.") else if(any(i %notin% x_rownames) || anyDuplicated(i) > 0L) warning("invalid or duplicated names in row index were removed.") # TODO
        i <- which(x_rownames %in% i)
      }
      # if column indices are names, check and convert to index numbers
      if(is.character(j)){
        x_colnames <- colnames(x)
        if(is.null(x_colnames)) stop("x does not have column names. Indexing by names failed.") else if(any(j %notin% x_colnames) || anyDuplicated(i) > 0L) warning("invalid or duplicated names in column index were removed.") # TODO
        j <- which(x_colnames %in% j)
      }
      
      # if drop was not specified and either index is exactly one, set to TRUE, else FALSE
      if(missing(drop)) drop <- if(length(i) == 1L || length(j) == 1L) TRUE else FALSE
      
      # yet another warning because discarding columns is a bad idea in this case
      # if(nay(seq_len(ncol(x)) %notin% j)) warning("extracting a subset of columns. Check if this was intended! If not, change the data elsewhere and use DR_data().") # TODO
      
      new_obj <- as.matrix(unclass(x))                    [i = i, j = j, drop = drop]
      new_o_o <- as.matrix(unclass(attr(x, "Y.original")))[i = i, j = j, drop = drop]
    }
  } else { stop("unexpected error!") }
  
  # restore attributes from x
  attr(new_obj, "Y.original" ) <- new_o_o
  attr(new_obj, "dims"       ) <- ncol(new_obj)
  attr(new_obj, "dim.names"  ) <- colnames(new_obj)
  attr(new_obj, "obs"        ) <- nrow(new_obj)
  attr(new_obj, "valid_obs"  ) <- sum(complete.cases(new_obj))
  attr(new_obj, "normalized" ) <- attr(x, "normalized" )
  attr(new_obj, "transformed") <- attr(x, "transformed")
  attr(new_obj, "base"       ) <- attr(x, "base"       )
  
  class(new_obj) <- "DirichletRegData"
  
  new_obj
}



str.DirichletRegData <- function(object, ...){
  NextMethod("str", ...)
}
