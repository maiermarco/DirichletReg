### 3 components - ternary
toTernary <- function(abc){
  cbind(
    x = (abc[, 1L] + 2.0 * abc[, 3L]) / sqrt(3.0),
    y = abc[, 1L]
  )
}


toTernaryVectors <- function(c1, c2, c3){ toTernary(cbind(c1, c2, c3)) }


### 4 components - quaternary
toQuaternary <- function(abcd){
  cbind(
    x = (abcd[, 1L] + 2.0 * abcd[, 3L] + abcd[, 4L]) / sqrt(3.0),
    y = abcd[, 1L] + abcd[, 4L] / 3.0,
    z = abcd[, 4L]
  )
}


toQuaternaryVectors <- function(c1, c2, c3, c4){ toQuaternary(cbind(c1, c2, c3, c4)) }


### "generic" function
toSimplex <- function(x){
  # checks
  if(is.null(dim(x))) stop("\"x\" must be a matrix-like object.")
  if((ncol(x) < 3L) || (ncol(x) > 4L)) stop("\"x\" must have 3 or 4 columns.")
  if(any((x < 0.0) | (x > 1.0))) stop("all values in \"x\" must be in [0, 1].")
  if(!isTRUE(all.equal(rowSums(x), rep(1.0, nrow(x)), check.attributes = FALSE))) stop("all row sums of \"x\" must be 1.")
  
  # transformations
  if(ncol(x) == 3L){
    return(toTernary(x))
  } else if(ncol(x) == 4L){
    return(toQuaternary(x))
  } else {
    stop("unexpected error")
  }
}
