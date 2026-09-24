DR_data <- function(
  Y,                                       # response (compositional variable)
  trafo       = sqrt(.Machine$double.eps), # transform (compress) the data?
  base        = 1L,                        # base variable for the reparametrized (alternative) model
  norm_tol    = sqrt(.Machine$double.eps), # tolerance for normalization [0, ?]
  no_guessing                              # disables "name guessing" when Y is a vector
){
  
  # initialization
  force.norm     <- isTRUE(norm_tol) # was normalization forced?
  force.norm.gt1 <- FALSE            # normalization because of the data
  state.tran     <- FALSE            # was Y transformed?
  force.tran     <- isTRUE(trafo)    # was transformation forced?
  
  if(length(trafo) != 1L || is.na(trafo) || !(is.logical(trafo) || (trafo > 0))){
    stop('"trafo" must be a small value > 0 or TRUE/FALSE. See ?DR_data') # error if trafo is not specified correctly
  }
  
  # set up beta-distributed matrix if a variables with values in [0, 1] is supplied
  if(is.null(dim(Y)) || (ncol(Y) == 1L)){ # if Y is a vector
    if(!is.null(dim(Y)) && ncol(Y) == 1L){ # handle Y with only 1 column
      if(!is.null(colnames(Y)) && nchar(colnames(Y)) > 1L) .name <- col # if that single column has a non-empty name, store it in .name
      Y <- Y[, 1L, drop = TRUE] # store the column in a vector format
    }    
    
    if((length(na.delete(Y)) < 1L) || any((na.delete(Y) < 0) | (na.delete(Y) > 1))){ # error if no non-missing or just values outside [0, 1] supplied
      stop('only one variable with values outside [0, 1] supplied.\nbeta distribution cannot safely be assumed.\ncheck and prepare your data first.')
    }
    
    Y <- cbind(1.0 - Y, Y) # create a two-column matrix
    
    # get the deparsed name from the call unless Y was a single column matrix with a column name
    if(!exists(".name")) .name <- deparse_nocutoff(match.call()$Y)
    
    # loop to debug regex
    # for(.name in c("ob[,1L]", "ob[ , 1L ]", "ob[[ 1L ]]", "ob$var", "ob$`var`", "ob[var]", "ob[[var]]", "ob[\"var\"]", "ob[[\"var\"]]", "ob['var']", "ob[['var']]")){
    .name <- local({ # here, we try to guess the name from the call   o_O
      .original.name <- .name # keep a backup in case anything goes wrong
      .name <- gsub("\\s", "", .name) # remove whitespace
      if(grepl("^.+\\[?\\[\\,?[0-9]+L?\\]\\]?$", .name)) return(.name) # shortening is too risky, returning as is
      .name <- gsub("^.+\\$", "", .name) # eliminate any references to the object the variable comes from (e.g., object$variable)
      if(grepl("^.+\\[\\[?.+\\]\\]?$", .name)){ # if .name looks something like object[["variables"]], ob['name'] etc.
        .name <- paste(strsplit(.name, split = "^.+\\[\\[?[\"']?|[\"']?\\]\\]?$")[[1L]], collapse = "") # split the string and try to extract the variable
      }
      if(grepl("^`(.+)`$", .name)) .name <- gsub("^`|`$", "", .name) # if Y was supplied as object$`variable` remove the backticks
      if(grepl("[[:alpha:]]+", .name)) return(.name) else return(.original.name) # if .name looks good, return it, else return the original
    })
    #writeLines(.name)} # debug loop end
    
    if(nchar(.name) < 1L) .name <- "Y" # fallback if, for some reason, .name is an empty string
    colnames(Y) <- c(paste0("(1 - ", .name, ")"), .name) # name the columns accordingly
    
    # print a message how the vector was processed
    message('only one variable in [0, 1] supplied - beta-distribution assumed.\ncheck this assumption.')
  }
  
  # set all rows containing NAs to NA
  if(anyNA(Y)){ Y[which(rowSums(is.na(Y)) > 0L), ] <- NA }
  
  # check if remaining matrix has at least 1 row
  if(nrow(na.delete(Y)) < 1L) stop('"Y" has no valid rows.')
  
  # check for negative values in Y
  if(any(na.delete(Y) < 0)) stop('"Y" contains values < 0.')
  
  # save the original data for reference
  Y.original <- Y
  
  
  
  # more checks
  if(is.null(dim(Y))) stop('"Y" must be either a matrix or a data.frame.') # this should not be possible
  if(ncol(Y) < 2L) stop('"Y" must at least have two columns.') # neither should this
  if(!is.integer(base) || (base < 1L) || (base > ncol(Y))) stop('"base" must be an integer in the range of variables.') # check base category
  if(length(norm_tol) != 1L || is.na(norm_tol) || (norm_tol <= 0)) stop('"norm_tol" must be a small number > 0. See ?DR_data')
  if(is.null(colnames(Y))) colnames(Y) <- paste0("v", seq_len(ncol(Y))) # if Y has no column names, assign a sequence v1, v2, v3, ...
  
  
  
  # Normalization - either by user-request or forced if row sums != 1 (with tolerance = norm_tol)
  row.sums <- rowSums(Y) # na.rm is irrelevant, because rows containing NAs have been set to NA above
  
  if(
    force.norm || # either forced by the user or
    isTRUE(all.equal(na.delete(row.sums), rep(1.0, length(na.delete(row.sums))), tolerance = norm_tol, check.attributes = FALSE))
  ){
    Y <- Y / row.sums # normalize rows
    force.norm.gt1 <- any(Y > 1, na.rm = TRUE) # was normalization necessary because of values over 1
  }
  
  
  
  # Transformation
  if(
    force.tran || # if either transformation is forced by the user or
    (is.numeric(trafo) && (any(Y < trafo, na.rm = TRUE) || any(Y > (1 - trafo), na.rm = TRUE))) # values are too close to 0 or 1 -- should be state.tran?
  ){
    n.obs      <- length(na.delete(row.sums))             # number of valid observations
    Y          <- (Y * (n.obs - 1) + 1 / ncol(Y)) / n.obs # Smithson, M. & Verkuilen, J. (2006)
    state.tran <- TRUE                                    # was Y transformed?
  }
  
  if(any(Y <= 0, na.rm = TRUE) || any(Y >= 1, na.rm = TRUE)){ # this should not happen, checking anyways
    stop('"trafo" was suppressed, yet values on the boundary of the support are present (0 and 1).\nConsider setting "trafo" = TRUE or to a threshold.\nSee ?DR_data')
  }
  
  
  
  # Object definition
  res <- structure(
    ".Data"       = as.matrix(Y),                 # the final, possible normalized/transformed data
    "Y.original"  = as.data.frame(Y.original),    # the original data
    "dims"        = ncol(Y),                      # the number of dimensions/components
    "dim.names"   = colnames(Y),                  # names of dimensions/components
    "obs"         = nrow(Y),                      # number of observations (including NAs)
    "valid_obs"   = length(na.delete(row.sums)),  # number of valid observations
    "normalized"  = force.norm || force.norm.gt1, # normalizations?
    "transformed" = force.tran || state.tran,     # transformation?
    "base"        = base,                         # index of the base category
    "class"       = "DirichletRegData"            # class definition
  )
  
  
    
  # Issue warnings
  if((force.norm || force.norm.gt1) && (force.tran || state.tran)){
    warning("not all rows sum up to 1 => normalization forced\n  some entries are 0 or 1 => transformation forced")
  } else if(force.norm || force.norm.gt1){
    warning("not all rows sum up to 1 => normalization forced")
  } else if(force.tran || state.tran){
    warning("some entries are 0 or 1 => transformation forced")
  }
  
  
  
  return(res)
}
