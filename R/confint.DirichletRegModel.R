confint.DirichletRegModel <- function(
  object
, parm
, level = 0.95
, ...
, type  = c("all", "beta", "gamma")
, exp   = FALSE
){
  e    <- exp
  type <- match.arg(type)

  if(any(level <= 0.0) || any(level >= 1.0)) stop("level must be in (0, 1)")
  level <- sort(level)

  repar <- if(object$parametrization == "alternative") TRUE else FALSE
  if(!repar && (type == "gamma")) warning("type = gamma not meaningful for the common parametrization - ignored")

  co <- object$coefficients
  se <- object$se
  names(co) <- object$coefnames
  names(se) <- object$coefnames

  res <- list(
      level        = level
    , type         = type
    , coefficients = coef(object)
    , se           = se
    , e            = e
    , repar        = repar
    )

  rci <- lapply(level, function(L){
    cbind(qnorm(    (1.0 - L) / 2.0, co, se),
          qnorm(L + (1.0 - L) / 2.0, co, se))
  })

  res$ci <- lapply(seq_along(rci), function(i) list())

  if(repar){
    Xc <- rev(rev(cumsum(c(1L, object$n.vars)))[-1L])
    Zc <- c(1L, ncol(object$Z)) + rev(Xc)[1L] - 1L

    for(ll in seq_along(rci)){
      inti <- 0L
      for(i in seq_len(object$dims)){
        res$ci[[ll]][i] <- if(i == object$base){
          list(NULL)
        } else {
          inti <- inti + 1L
          list(rci[[ll]][Xc[inti]:(Xc[inti + 1L] - 1L), , drop = FALSE])
        }
      }
      res$ci[[ll]][[length(res$ci[[ll]]) + 1L]] <- rci[[ll]][Zc[1L]:Zc[2L] , , drop = FALSE]
      names(res$ci[[ll]]) <- c(object$varnames, "gamma")
    }
  } else {
    Xc <- cumsum(c(1L, object$n.vars))
    for(ll in seq_along(rci)){
      inti <- 0L
      for(i in 1L:object$dims){
        inti <- inti + 1L
        res$ci[[ll]][[i]] <- rci[[ll]][Xc[inti]:(Xc[inti + 1L] - 1L) , , drop = FALSE]
      }
    }
  }

  if(e){
    for(ll in seq_along(rci)){
      for(these in which(!unlist(lapply(res$ci[[1L]], is.null)))){
        res$ci[[ll]][[these]] <- exp(res$ci[[ll]][[these]])
      }
    }
  }

  class(res) <- "DirichletRegConfint"

  return(res)

}
